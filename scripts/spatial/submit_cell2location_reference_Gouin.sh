#!/bin/bash

source ~/.bashrc

micromamba activate cell2location

python cell2location_reference.py -s gouin
