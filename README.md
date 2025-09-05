# UEFEX
The extensional package for the UEF package in LAMMPS. This package is useful for Langevin dynamics and dissipative particle dynamics to apply uniform extensional flows.

<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/movie-2022-04-18.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

Authored by:
[Takahiro Murashima](https://github.com/t-murash)<br>
Tohoku University, Japan<br>
Initial commit: Feb 22, 2018<br>
Last updated: Sep 5, 2025<br>
Support provided via [issues](https://github.com/t-murash/LAMMPS-UEFEX/issues) and/or [email](mailto:murasima@cmpt.phys.tohoku.ac.jp).

<!--
<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/movie.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/original-view.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

<img src="https://github.com/t-murash/LAMMPS-UEFEX/blob/master/img/cubic-view.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>
-->



## Installation
This package is compiled within LAMMPS and depends on the original UEF package.
Download and install LAMMPS and the UEF package according to the following sites.
* [LAMMPS](https://lammps.org/)
* [RutledgeGroupMIT/UEF](https://github.com/RutledgeGroupMIT/UEF)

**(Note)** The UEF package has been included in LAMMPS several years ago. You **DO NOT** need to get the UEF package from Rutledge Group's site.


### Build LAMMPS with cmake

```
wget https://download.lammps.org/tars/lammps-22Jul2025.tar.gz
tar zxvf lammps-22Jul2025.tar.gz
git clone https://github.com/t-murash/LAMMPS-UEFEX
cp -r LAMMPS-UEFEX/UEFEX lammps-*/src/.
mv lammps-*/src/UEFEX/domain.cpp lammps-*/src/.
cp LAMMPS-UEFEX/cmake/CMakeLists.txt lammps-*/cmake/.
cp LAMMPS-UEFEX/cmake/presets/uefex.cmake lammps-*/cmake/presets/.

```
Then, build
```
cd lammps-*
mkdir build
cd build
cmake ../cmake -DBUILD_MPI=yes -C ../cmake/presets/uefex.cmake
make
```


### Build LAMMPS with make
```
wget https://download.lammps.org/tars/lammps-22Jul2025.tar.gz
tar xvf lammps-22Jul2025.tar.gz
cd lammps-*/src/
make yes-molecule
make yes-uef
```
(`yes-molecule` is necessary for examples.)

Then, get this package and install.

```
git clone https://github.com/t-murash/LAMMPS-UEFEX.git
cp -r LAMMPS-UEFEX/UEFEX lammps-*/src/.
mv lammps-*/src/UEFEX/domain.cpp lammps-*/src/.
cd lammps-*/src/
make yes-uefex
```

Finally, compile LAMMPS.

```
cd lammps-*/src/
make mpi mode=static
```



### For old version users (23Jun2022, 29Sep2021, 29Oct20, 3Mar20)
You can find `2021-Sep-29`, `2020-Oct-29`, and `2020-Mar-3` directories in `LAMMPS-UEFEX/UEFEX/old`.
`2022-Jun-23` contains source files compatible with `23Jun2022` of LAMMPS.
`2021-Sep-29` contains source files compatible with `29Sep2021` of LAMMPS.
`2020-Oct-29` contains source files compatible with `29Oct2020` of LAMMPS.
The source files in `2020-Mar-3` are compatible with `3Mar20` of LAMMPS and the more previous versions of LAMMPS.
Copy the source files to the `lammps-*/src/UEFEX` directory according to your versions.
Then, make or cmake in the same way as described above.

**(Note)** These old versions are not up-to-date.

## Usage
You can find several example files in `LAMMPS-UEFEX/examples`.
```
mpirun ./lmp_mpi -in in.example
```
or
```
mpirun ./lmp -in in.example
```


## Citing the UEFEX package

Users of this package are encouraged to cite the following articles in scientific publications:

* D. A. Nicholson, G. C. Rutledge, "Molecular simulation of flow-enhanced nucleation in *n*-eicosane melts under steady shear and uniaxial extension", *J. Chem Phys.* (2016) **145** (24) 244903, https://doi.org/10.1063/1.4972894.

* T. Murashima, K. Hagita, T. Kawakatsu, "Elongational Viscosity of Weakly Entangled Polymer Melt via Coarse-Grained Molecular Dynamics Simulation", *Nihon Reoroji Gakkaishi (J. Soc. Rheol. Jpn.)*  (2018) **46** (5) 207-220, https://doi.org/10.1678/rheology.46.207.
