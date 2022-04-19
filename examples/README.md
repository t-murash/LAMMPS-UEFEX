# UEFEX/examples

## Contents

| file name | explanation |
| ---- | ----------- |
| N100M100.data | Configuration data of Kremer-Grest model with N=100 and M=100 |
| README.md | This document |
| anim.py   | Python script for processing a snapshot of polymers (PNG) |
| anim.sh   | Bash script for processing an animation of polymers (gif) |
| in.uefex  | Example for UEFEX(nve/uefex, temp/uefex, pressure/uefex) |
| in.uefex.data          | Input file to produce the above gif file |
| in.uefex.dump          | Example for rotation/uefex |
| in.uefex.eng           | Example for engrate option |
| in.uefex.stress.atom   | Example for stress/atom/uefex |
| in.uefex.uni           | Example for fix_modify u option |
| in.uefex.uni.integrate | Example for fix_modify ui option |
| smooth.py              | Python script for smoothing data |
| visc.plt               | Gnuplot script for processing a graph of viscosity growth |



## GIF animation
First, install OVITO python module.
- [OVITO](https://www.ovito.org/python-downloads/)

Then, execute the following commands.
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


## Elongational viscosity
First, install `numpy`, `scipy`, `pandas`, and `gnuplot`.

Then, execute the following commands.
```
mpirun ./lmp -in in.uefex
python smooth.py
gnuplot visc.plt
```
You will get the following image (PNG file).

<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/visc.png" title="Elongational viscosity" width=300/>

`mpirun ./lmp -in in.uefex` produces `press.txt`. The first 4 lines of this text file (`press.txt`) will start as follows.
```
# Time-averaged data for fix a1
# TimeStep v_pxx v_pyy v_pzz v_pxy v_pxz v_pyz
0 4.81638 4.70122 5.01902 0.0253983 0.0392435 0.0244829
100 4.63578 4.60468 4.80104 0.120378 -0.029315 -0.0383236
```
The first line is the comment line.
Although this line says "Time-averaged data", the data (`press.txt`) are not time-averaged because of `${freq} 1 ${freq}` at `fix ave/time` in `in.uefex`.

The first column is "Time Step". To obtain "Time", you need to multiply <img src="https://render.githubusercontent.com/render/math?math=\Delta t (=0.01)"> to this column values.
The second to seventh columns present the components of "Pressure tensor". "Stress tensor" is the negative value of "Pressure tensor" <img src="https://render.githubusercontent.com/render/math?math=\sigma=-P">.
Uniaxial elongational viscosity is calculated by
<img src="https://render.githubusercontent.com/render/math?math=\eta_{\rm uni}=\{\sigma_{zz}-(\sigma_{xx} %2B \sigma_{yy})/2 \} / \dot{\varepsilon}">.
Here, <img src="https://render.githubusercontent.com/render/math?math=\dot{\varepsilon}={\rm d}\varepsilon / {\rm d}t (=0.001)"> is the elongational rate.

Since the raw data of the elongational viscosity is noisy,
we apply Savitzky-Golay filter to smooth out the high frequency noise through `python smooth.py`.
Finally, we will get the above graph by `gnuplot visc.plt`.
