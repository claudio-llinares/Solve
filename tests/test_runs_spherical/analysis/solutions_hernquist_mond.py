#!/usr/bin/env python
# coding: utf-8

import numpy as np
np.set_printoptions(linewidth=210)  # 210

import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec

from scipy.integrate import odeint
from scipy.integrate import solve_ivp


# Redshift is needed because $a$ factors appear in the equations
# (although, they are not included here, so all this is at redshift 0).

class hernquist_solver_mond_with_mass:

    def __init__(self):

        self.cm2mpc  = 1.0/3.08568025e24
        self.km2mpc  = 1.0/3.08568025e19
        self.m2mpc   = 1.0/(3.08568025e22)
        self.sec2gyr = 1.0 / (1.0e+9 * 365.25 * 23.56 * 3600)
        self.kg2msun = 1.0/(1.9198e+30)

    def init_hernquist_solver_mond_with_mass(self):

        # MOND parameters:
        self.C       = 3.0*self.omegam*(self.H0*self.km2mpc/self.sec2gyr)**2/(2.0-self.KB)  # constant in front of the density in Poisson's equation
        self.A       = (3.0*self.omegam*(self.H0*self.km2mpc/self.sec2gyr)**2*self.aH**2*self.ratio_rho)/4.0  # some constant that appears in the solutions

        self.G_mks    = 6.67428e-11   # m^3 kg^-1 s^-2
        self.rho_crit = (3.0*(self.H0*self.km2mpc)**2)/(8.0*np.pi*self.G_mks*self.m2mpc**3/self.kg2msun)
        self.rho_mean = self.omegam*self.rho_crit

        print("a0       = ", self.a0, "mpc/gyr^2")
        print("H0       = ", self.H0, "1/gyr")
        print("C        = ", self.C, "1/gyr^2")
        print("rho_mean = %e" % self.rho_mean, "msun/mpc**3")

        # Output file:
        self.output_file = "output_z_%f" % self.z

        # Boundaries of the domain and initial conditions:
        self.domain  = [self.r_left, self.r_right]
        self.ics = [self.phi_bar(self.r_left), self.grad_phi_bar(self.r_left), self.phi_bar(self.r_left)]  # initial conditions for the solution

        # Some more initialization:
        self.t_eval = np.logspace(np.log10(self.domain[0]), np.log10(self.domain[1]),1000)  # points where the numerical solution will be clculated
        self.r_plot = np.logspace(np.log10(self.r_left), np.log10(self.r_right), 1000)  # radius for plotting

    # Density:
    def delta(self, r):

        return self.ratio_rho/((r/self.aH)*(1.0+r/self.aH)**3) - 1.0

    # Analytic solutions with no mass term:
    #======================================
    def phi_newt(self, r):

        return -self.A/(1.0+r/self.aH)

    def grad_phi_bar(self, r):

        return self.A/(self.aH*(1.0+r/self.aH)**2)

    def grad_chi(self, r):

        #return (self.A*self.aH + np.sqrt(self.A)*np.sqrt(self.aH)*np.sqrt(self.A*self.aH + 4*self.a0*(self.aH + r)**2))/(2.*(self.aH + r)**2)
        return (self.A*self.aH + np.sqrt(self.A)*np.sqrt(self.aH)*np.sqrt(self.A*self.aH))/(2.*(self.aH + r)**2)

    def grad_phi(self, r):

        return self.grad_phi_bar(r) + self.grad_chi(r)

    def phi_bar(self, r):
        # This one is the solution of Poisson's equation for the 
        # Hernquist profile which can be used to determine initial conditions
        # assuming that the profile approaches the Newtonian regimen in the center
        # (which may not necessarily happen)
        return -self.A/(1.0+r/self.aH)

    def integrand(self, r, psi):

        # Integrand for the model with mass term.
        # psi[0] = phi_bar
        # psi[1] = u_phi
        # psi[2] = chi

        integrand_0 = psi[1]
        integrand_1 = self.C*self.delta(r) - 2.0/r*psi[1] - self.M2*(psi[0]+psi[2])
        integrand_2 = 0.5*(psi[1]+np.sqrt(psi[1]**2+4.0*self.a0*psi[1]))

        return [integrand_0, integrand_1, integrand_2]

    def solve_hernquist(self):

        # Density:
        delta_plot    = self.delta(self.r_plot)
        phi_newt_plot = self.phi_newt(self.r_plot)

        # Solutions for M=0 and M2:
        tmp = self.M2
        self.M2 = 0.0
        psi_zero = solve_ivp(self.integrand, self.domain,  self.ics, method='RK45', rtol=1.0e-10, atol=1.0e-10, t_eval=self.t_eval, dense_output="True")
        self.M2 = tmp
        psi      = solve_ivp(self.integrand, self.domain,  self.ics, method='RK45', rtol=1.0e-10, atol=1.0e-10, t_eval=self.t_eval, dense_output="True")

        return self.r_plot, delta_plot, phi_newt_plot, psi_zero, psi

    def store_solution_hernquist(self, psi_zero):
        # Store solution:
        # It may be necessary to add $a$ factors to compare this solution with Ramses.
        # (Ramses fields contain $a$ factors in their definitions)
        np.savetxt(self.output_file, np.column_stack([psi_zero.t, psi_zero.y[0], psi_zero.y[2], psi_zero.y[0]+psi_zero.y[2]]))


    # Plot solutions:
    #===============
    def plot_solutions(self, psi_zero, psi):

        plt.figure(figsize=(15,10))
        plt.rcParams.update({'font.size': 15})  # 15
        #gs1 = gridspec.GridSpec(nrows=1, ncols=1, wspace=0.3, hspace=0.3) # The plots

        plt.subplot(221)
        plt.plot(self.r_plot, self.delta(self.r_plot), lw=3)
        plt.axvline(1.0/self.M2, c="k")
        plt.grid()
        plt.xlabel("r (Mpc)")
        plt.ylabel("Overdensity")
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()

        plt.subplot(222)
        sol1 = [ self.grad_phi_bar(r) for r in self.r_plot]
        sol2 = [ self.grad_chi(r) for r in self.r_plot]
        sol3 = [ self.grad_phi(r) for r in self.r_plot]
        plt.plot(self.r_plot, np.multiply(self.r_plot,sol1),                                 label="d$\\tilde\Phi$/d$r$ (no mass)", lw=6  )  
        plt.plot(self.r_plot, np.multiply(self.r_plot,sol2),                                 label="d$\chi$/d$r$ (no mass)",        lw=6  ) 
        plt.plot(self.r_plot, np.multiply(self.r_plot,sol3),                                 label="d$\Phi$/d$r$ (no mass)",        lw=6  ) 
        plt.plot(psi.t,  np.multiply(psi.t,psi.y[1]),                              label="d$\\tilde\Phi$/d$r$ (mass)",    lw=3  )
        plt.plot(psi.t,  np.multiply(psi.t,self.integrand(psi.t, psi.y)[2]),            label="d$\chi$/d$r$ (mass)",    lw=3  ) 
        plt.plot(psi.t,  np.multiply(psi.t,psi.y[1] + self.integrand(psi.t, psi.y)[2]), label="d$\Phi$/d$r$ (mass)",    lw=3  ) 
        plt.axvline(1.0/self.M2, c="k")
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r (Mpc)')
        plt.ylabel("Circular velocity (Mpc/Gyr)")
        plt.grid(True)
        plt.legend()   

        plt.subplot(223)
        plt.plot(psi_zero.t, psi_zero.y[0],               label="$\\tilde{\Phi}$ (mass)", lw=6  )
        plt.plot(psi_zero.t, psi_zero.y[2],               label="$\chi$ (mass)", lw=6  )
        plt.plot(psi_zero.t, psi_zero.y[0]+psi_zero.y[2], label="$\Phi$ (mass)", lw=6  )
        plt.plot(psi.t, psi.y[0],                         label="$\\tilde{\Phi}$ (no mass)", lw=3  )
        plt.plot(psi.t, psi.y[2],                         label="$\chi$ (no mass)", lw=3  )
        plt.plot(psi.t, psi.y[0]+psi.y[2],                label="$\Phi$ (no mass)", lw=3  )
        #plt.ylim(-10.0, 25.0)
        plt.axvline(1.0/self.M2, c="k")
        plt.xscale('log')
        #plt.yscale('log')
        plt.xlabel('r (Mpc)')
        plt.ylabel("Solutions (Mpc$^2$/Gyr$^2$)")
        plt.grid(True)
        plt.legend()   

        plt.subplot(224)
        sol1 = [ self.grad_phi_bar(r) for r in self.r_plot]
        sol2 = [ self.grad_chi(r) for r in self.r_plot]
        sol3 = [ self.grad_phi(r) for r in self.r_plot]
        plt.plot(self.r_plot, sol1, label="d$\\tilde\Phi$/d$r$ (no mass)", lw=6  )  
        plt.plot(self.r_plot, sol2, label="d$\chi$/d$r$ (no mass)", lw=6  ) 
        plt.plot(self.r_plot, sol3, label="d$\Phi$/d$r$ (no mass)", lw=6  ) 

        plt.plot(psi_zero.t, psi_zero.y[1], label="d$\\tilde\Phi$/d$r$ (mass)", lw=3  )
        plt.plot(psi_zero.t, self.integrand(psi_zero.t, psi_zero.y)[2], label="d$\\tilde\Phi$/d$r$ (mass)", lw=3  ) 
        plt.plot(psi_zero.t, psi_zero.y[1] + self.integrand(psi_zero.t, psi_zero.y)[2], label="d$\\tilde\Phi$/d$r$ (mass)", lw=3  ) 

        plt.plot(psi.t, psi.y[1], label="d$\\tilde\Phi$/d$r$ (mass)", lw=3  )
        plt.plot(psi.t, self.integrand(psi.t, psi.y)[2], label="d$\\tilde\Phi$/d$r$ (mass)", lw=3  ) 
        plt.plot(psi.t, psi.y[1] + self.integrand(psi.t, psi.y)[2], label="d$\\tilde\Phi$/d$r$ (mass)", lw=3  ) 
        plt.plot(self.r_plot, [self.a0]*len(self.r_plot), c="k")
        plt.plot(self.r_plot, [-self.a0]*len(self.r_plot), c="k", label="")
        plt.axvline(1.0/self.M2, c="k")
        #plt.ylim(-10.0, 25.0)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r (Mpc)')
        plt.ylabel("Derivatives (Mpc/Gyr$^2$)")
        plt.grid(True)
        plt.legend()   

        plt.tight_layout()

        plt.savefig("solutions_hernquist.pdf", dpi=300, format="pdf")
        plt.show()




