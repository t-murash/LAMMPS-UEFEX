#!/bin/bash

mpirun ./lmp -in in.uefex
pid=$!
wait $pid
python smooth.py
pid=$!
wait $pid
gnuplot visc.plt
