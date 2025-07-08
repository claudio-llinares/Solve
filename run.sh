#!/bin/bash

# Bring the source:
cd sphinx-c-apidoc
./run.sh
cd ../

# Build the documentation:
make html
make html



