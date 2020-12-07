# USER-UEFEX
The extensional package for the USER-UEF package in LAMMPS. This package is useful for Langevin dynamics and dissipative particle dynamics to apply uniform extensional flows.

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/movie.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/original-view.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/cubic-view.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>

<img src="https://github.com/t-murash/USER-UEFEX/blob/master/img/unwrap.gif" title="M=100, N=100 Kremer-Grest chains in a uniaxial elongational flow" width=300/>



Authored by:
[Takahiro Murashima](https://github.com/t-murash)<br>
Tohoku University, Japan<br>
Initial commit: Feb 22, 2018<br>
Last updated: Dec 07, 2020<br>
Support provided via [issues](https://github.com/t-murash/USER-UEFEX/issues) and/or [email](mailto:murasima@cmpt.phys.tohoku.ac.jp).

## Installation
This package is compiled within LAMMPS and depends on the original USER-UEF package.
Download and install LAMMPS and the USER-UEF package according to the following sites.
* [LAMMPS](https://lammps.sandia.gov/)
* [RutledgeGroupMIT/UEF](https://github.com/RutledgeGroupMIT/UEF)
(The USER-UEF package has already been included in LAMMPS.)

```
wget http://lammps.sandia.gov/tars/lammps-stable.tar.gz
tar xvf lammps-stable.tar.gz
cd lammps-*/src/
make yes-molecule
make yes-user-uef
```
(`yes-molecule` is necessary for examples.)

Then, get this package and install.

```
git clone https://github.com/t-murash/USER-UEFEX.git
cp -r USER-UEFEX/USER-UEFEX lammps-*/src/.
cd lammps-*/src/
make yes-user-uefex
```

Finally, compile LAMMPS.

```
cd lammps-*/src/
make mpi
```

## Usage
You can find several example files in `USER-UEFEX/examples`.
```
mpirun ./lmp_mpi -in in.example
```

## Citing the USER-UEFEX package

Users of this package are encouraged to cite the following articles in scientific publications:

* D. A. Nicholson, G. C. Rutledge, "Molecular simulation of flow-enhanced nucleation in *n*-eicosane melts under steady shear and uniaxial extension", *J. Chem Phys.*, **145** (24), 244903 (2016), http://aip.scitation.org/doi/full/10.1063/1.4972894.

* T. Murashima, K. Hagita, T. Kawakatsu, "Elongational Viscosity of Weakly Entangled Polymer Melt via Coarse-Grained Molecular Dynamics Simulation", *Nihon Reoroji Gakkaishi (J. Soc. Rheol. Jpn.)* , **46** (5), 207-220 (2018), https://doi.org/10.1678/rheology.46.207.
