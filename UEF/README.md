# compute_pressure_uef.cpp

`compute pressure/uef` failed to calculate `bond`, `angle`, `dihedral`, etc. This bug is fixed in `compute_pressure_uef.cpp`.

## Instalation

Copy `compute_pressure_uef.cpp` here into `${lammps source directory}/UEF`.

Then, `cmake` or `make`.

## Usage

```
compute uef_temp all temp/uef
compute uef_press all pressure/uef uef_temp
compute uef_press_bond all pressure/uef uef_temp bond
thermo_style custome step c_uef_press_bond[*]
```
