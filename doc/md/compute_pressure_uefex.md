<div class="index">

compute pressure/uefex

</div>

# compute pressure/uefex command

## Syntax

``` LAMMPS
compute ID group-ID pressure/uefex temp-ID keyword ...
```

- ID, group-ID are documented in [compute](compute.md) command
- pressure/uefex = style name of this compute command
- temp-ID = ID of a temp/uefex compute, or NULL
- zero or more keywords may be appended
- keyword = *ke* or *pair* or *bond* or *angle* or *dihedral* or
  *improper* or *kspace* or *fix* or *virial*

## Examples

``` LAMMPS
compute mypress all pressure/uefex uefex_temp
compute mypress all pressure/uefex uefex_temp virial
```

## Description

This compute calculates the pressure tensor in the laboratory (LAB)
frame during UEFEX extensional flow simulations. The virial
contributions and kinetic energy tensor are rotated from the
upper-triangular (UT) LAMMPS frame to the LAB frame using the rotation
matrix Q from [fix nve/uefex](fix_nve_uefex.md).

The 6-component pressure tensor is returned in Voigt notation: Pxx, Pyy,
Pzz, Pxy, Pxz, Pyz, where the directions correspond to the LAB frame.

For uniaxial extension (`erate -e/2 -e/2` where e is the strain rate), the extensional stress is:

$$N_1 = \sigma_{zz} - \frac{\sigma_{xx} + \sigma_{yy}}{2}
    = -P_{zz} + \frac{P_{xx} + P_{yy}}{2}$$

This compute is automatically created by [fix
nve/uefex](fix_nve_uefex.md) with the ID `{fix-ID}_press`.

The keywords and general output information are the same as for [compute
pressure](compute_pressure.md).

## Restrictions

This compute is part of the UEFEX package. It is only enabled if LAMMPS
was built with that package.

This compute can only be used when [fix nve/uefex](fix_nve_uefex.md) is
active.

The kinetic contribution to the pressure tensor will be accurate only
when the compute specified by *temp-ID* is a [compute
temp/uefex](compute_temp_uefex.md).

## Related commands

[compute pressure](compute_pressure.md), [compute
pressure/uef](compute_pressure_uef.md), [fix
nve/uefex](fix_nve_uefex.md), [compute
temp/uefex](compute_temp_uefex.md)

## Default

none
