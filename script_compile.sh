#!/bin/bash

# our comment is here
squeue -u rameshrp &&

nvcc sourcecode.cu -o TestRun &&

./TestRun

