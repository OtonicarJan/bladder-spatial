#!/bin/bash

source ~/.bashrc

micromamba activate cell2location

python3 cell2location_visium.py -a 20 -n 20 -c Gouin_muscle_merged