# USER-UEFEX
The extensional package for the USER-UEF package in LAMMPS. This package is useful for Langevin dynamics and dissipative particle dynamics to apply uniform extensional flow.

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/gif/movie.gif" title="N=100, M=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>



Authored by:
[Takahiro Murashima](https://github.com/t-murash)<br>
Tohoku University, Japan<br>
Initial commit: Feb 22, 2018<br>
Support provided via [issues](https://github.com/t-murash/USER-UEFEX/issues) and/or [email](mailto:murasima@cmpt.phys.tohoku.ac.jp).

## Installation
This package is compiled within LAMMPS and depends on the original USER-UEF package.
Download and install LAMMPS and the USER-UEF package according to the following site.

* [RutledgeGroupMIT/UEF](https://github.com/RutledgeGroupMIT/UEF)

```
wget http://lammps.sandia.gov/tars/lammps-stable.tar.gz
tar xvf lammps-stable.tar.gz
git clone https://github.com/RutledgeGroupMIT/UEF.git
cp -r UEF-master/USER-UEF/ lammps-*/src
cd lammps-*/src/
make yes-user-uef
```

Then, get this package and install.

```
git clone https://github.com/t-murash/USER-UEFEX.git
cp -r USER-UEFEX-master/USER-UEFEX/ lammps-*/src
cd lammps-*/src/
make yes-user-uefex
```

Finally, compile LAMMPS.

```
cd lammps-*/src/
make mpi
```

## Usage
You can find example files in `USER-UEF-master/examples`.
```
mpirun ./lmp_mpi -in in.example
```

## Citing the USER-UEFEX package

Users of this package are encouraged to cite the following articles in scientific publications:

* D. A. Nicholson, G. C. Rutledge, "Molecular simulation of flow-enhanced nucleation in *n*-eicosane melts under steady shear and uniaxial extension", *J. Chem Phys.*, 2016, **145** (24), http://aip.scitation.org/doi/full/10.1063/1.4972894.

* T. Murashima, K. Hagita, T. Kawakatsu, "Elongational viscosity of weakly entangled polymer melt via coarse-grained molecular dynamics simulation", in preparation.
