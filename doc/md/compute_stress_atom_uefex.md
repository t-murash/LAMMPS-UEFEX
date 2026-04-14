<div class="index">

compute stress/atom/uefex

</div>

# compute stress/atom/uefex command

## Syntax

``` LAMMPS
compute ID group-ID stress/atom/uefex temp-ID keyword ...
```

- ID, group-ID are documented in [compute](compute.md) command
- stress/atom/uefex = style name of this compute command
- temp-ID = ID of a temperature compute for kinetic energy contribution,
  or NULL
- zero or more keywords may be appended
- keyword = *ke* or *pair* or *bond* or *angle* or *dihedral* or
  *improper* or *kspace* or *fix* or *virial*

## Examples

``` LAMMPS
compute mystress all stress/atom/uefex uefex_temp
compute mystress all stress/atom/uefex NULL virial
```

## Description

This compute calculates the per-atom stress tensor in the laboratory
(LAB) frame during UEFEX extensional flow simulations. It is analogous
to [compute stress/atom](compute_stress_atom.md), but the stress tensor
is rotated from the upper-triangular (UT) LAMMPS frame to the LAB frame
using the rotation matrix Q from [fix nve/uefex](fix_nve_uefex.md).

The per-atom stress tensor has 6 components in Voigt notation: xx, yy,
zz, xy, xz, yz, where the directions correspond to the LAB frame.

The keywords are the same as for [compute
stress/atom](compute_stress_atom.md).

## Output info

This compute calculates a per-atom array with 6 columns. The values are
in stress\*volume units. See [compute
stress/atom](compute_stress_atom.md) for details on units and usage.

## Restrictions

This compute is part of the UEFEX package. It is only enabled if LAMMPS
was built with that package.

This compute can only be used when [fix nve/uefex](fix_nve_uefex.md) is
active.

## Related commands

[compute stress/atom](compute_stress_atom.md), [fix
nve/uefex](fix_nve_uefex.md), [compute
pressure/uefex](compute_pressure_uefex.md)

## Default

none
