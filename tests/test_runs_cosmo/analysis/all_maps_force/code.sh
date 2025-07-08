#!/bin/bash

# To make movies:
#mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc copy -lavcopts vcodec=mpeg4
mencoder "mf://*.png" -mf fps=25 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4






