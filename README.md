# UEFEX
The extensional package for the UEF package in LAMMPS. This package is useful for Langevin dynamics and dissipative particle dynamics to apply uniform extensional flows.

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/unwrap.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

Authored by:
[Takahiro Murashima](https://github.com/t-murash)<br>
Tohoku University, Japan<br>
Initial commit: Feb 22, 2018<br>
Last updated: Oct 20, 2021<br>
Support provided via [issues](https://github.com/t-murash/USER-UEFEX/issues) and/or [email](mailto:murasima@cmpt.phys.tohoku.ac.jp).


<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/movie.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/original-view.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/cubic-view.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>



## Installation
This package is compiled within LAMMPS and depends on the original UEF package.
Download and install LAMMPS and the UEF package according to the following sites.
* [LAMMPS](https://lammps.sandia.gov/)
* [RutledgeGroupMIT/UEF](https://github.com/RutledgeGroupMIT/UEF)
(The UEF package has already been included in LAMMPS. You do not need to get the UEF package from this site.)

```
wget http://lammps.sandia.gov/tars/lammps-stable.tar.gz
tar xvf lammps-stable.tar.gz
cd lammps-*/src/
make yes-molecule
make yes-uef
```
(`yes-molecule` is necessary for examples.)

Then, get this package and install.

```
git clone https://github.com/t-murash/USER-UEFEX.git
cp -r USER-UEFEX/UEFEX lammps-*/src/.
cd lammps-*/src/
make yes-uefex
```

Finally, compile LAMMPS.

```
cd lammps-*/src/
make mpi mode=static
```

## For cmake users

```
wget http://lammps.sandia.gov/tars/lammps-stable.tar.gz
tar xvf lammps-stable.tar.gz
git clone https://github.com/t-murash/USER-UEFEX.git
cp -r USER-UEFEX/UEFEX lammps-*/src/.
mv lammps-*/src/UEFEX/domain.cpp lammps-*/src/.
```

Edit the following files to include `UEFEX`.

```
lammps-*/cmake/CMakeLists.txt
lammps-*/cmake/presets/all_off.cmake
lammps-*/cmake/presets/all_on.cmake
```
The easiest way is, find `UEF` in the above files and place `UEFEX` below `UEF`.


Then, build using `cmake` with `-D PKG_UEFEX=yes`

```
cd lammps-*
mkdir build
cd build
cmake ../cmake -D BUILD_MPI=yes　-D PKG_MOLECULE=yes -D PKG_UEF=yes -D PKG_UEFEX=yes
make
```

## Usage
You can find several example files in `USER-UEFEX/examples`.
```
mpirun ./lmp_mpi -in in.example
```
or
```
mpirun ./lmp -in in.example
```


## Citing the UEFEX package

Users of this package are encouraged to cite the following articles in scientific publications:

* D. A. Nicholson, G. C. Rutledge, "Molecular simulation of flow-enhanced nucleation in *n*-eicosane melts under steady shear and uniaxial extension", *J. Chem Phys.*, **145** (24), 244903 (2016), http://aip.scitation.org/doi/full/10.1063/1.4972894.

* T. Murashima, K. Hagita, T. Kawakatsu, "Elongational Viscosity of Weakly Entangled Polymer Melt via Coarse-Grained Molecular Dynamics Simulation", *Nihon Reoroji Gakkaishi (J. Soc. Rheol. Jpn.)* , **46** (5), 207-220 (2018), https://doi.org/10.1678/rheology.46.207.
