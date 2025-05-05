#!/bin/bash

# Install tools
pip install pyani

# Calculate ANI distance between strain pairs
average_nucleotide_identity.py -i ./genomes/ -o ./ani_output -m ANIm -g  # Note：This step is super time-consuming!
