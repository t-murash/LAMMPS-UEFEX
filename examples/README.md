# UEFEX/examples

## Contents
- N100M100.data : Configuration data of Kremer-Grest model with N=100 and M=100
- README.md : This document
- anim.py : Python script for processing a PNG file
- anim.sh : Bash script for processing a gif file
- in.uefex : Example for UEFEX(nve/uefex, temp/uefex, pressure/uefex)
- in.uefex.data : Input file to produce the above gif file
- in.uefex.dump : Example for rotation/uefex
- in.uefex.eng : Example for engrate option
- in.uefex.stress.atom : Example for stress/atom/uefex
- in.uefex.uni : Example for fix_modify u option
- in.uefex.uni.integrate : Example for fix_modify ui option

## Elongational viscosity
First, install `numpy`, `scipy`, `pandas`, and `gnuplot`.

Then, run the following commands.
```
mpirun ./lmp -in in.uefex
python smooth.py
gnuplot visc.plt
```
You will get the following image (PNG file).

<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/visc.png" title="Elongational viscosity" width=300/>

## GIF animation
First, install OVITO python module.
- [OVITO](https://www.ovito.org/python-downloads/)

Then, run the following commands.
```
mpirun ./lmp -in in.uefex.data
bash anim.sh
```

You will get data files, png files, and a gif file (movie.gif) as shown in the top page.

<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/movie-2022-04-18.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

**movie.gif**

<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/figure.000.png" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

**figure.000.png**

`in.uefex.data` summarizes how to use the UEFEX package.
```
fix 2 all nve/uefex erate 0.0 0.0          # set nve with UEF deformation (elongation rate is 0)
compute mytemp all temp/uefex              # compute temperature under elongational flow
compute mypress all pressure/uefex mytemp  # compute pressure under elongational flow (pxx,pyy,pzz,pxy,pxz,pyz order)
fix_modify 1 temp mytemp                   # velocity correction under elongational flow
fix_modify 2 u 0.001                       # set uniaxial elongational rate as 0.001 (u:uniaxial, b:biaxial, p:planar)
compute rmatrix all rotation/uefex         # compute rotation matrix
```

Particle positions and the unit cell saved in a LAMMPS data file are not consistent with the elongational directions under UEF.
To correct this, a rotation matrix is multiplied to the particle positions and the unit cell in `anim.py`.