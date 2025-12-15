#!/bin/bash

# cs7
mkdir -p cs7/data
rig run -f cs7/r/simulation.r
rig run -f cs7/r/jags.r
rig run -f cs7/r/stan.r
rig run -f cs7/r/analysis.r