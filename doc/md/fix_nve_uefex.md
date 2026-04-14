<div class="index">

fix nve/uefex

</div>

# fix nve/uefex command

## Syntax

``` LAMMPS
fix ID group-ID nve/uefex erate edot_x edot_y keyword value ...
```

- ID, group-ID are documented in [fix](fix.md) command

- nve/uefex = style name of this fix command

- erate = required keyword

- edot_x, edot_y = strain rates in x and y directions (1/time units)

- one or more keyword/value pairs may be appended

  <div class="parsed-literal">

  keyword = *erate* or *engrate* or *strain*  
  *erate* values = edot_x edot_y = true strain rates (required unless
  engrate is used) *engrate* value = edot_eng = engineering strain rate
  for uniaxial extension *strain* values = e_x e_y = initial strain (for
  resuming from data files)

  </div>

## Examples

``` LAMMPS
fix uefex all nve/uefex erate -0.0005 -0.0005
fix uefex all nve/uefex erate 0.0005 -0.001
fix uefex all nve/uefex engrate 0.001
fix uefex all nve/uefex erate -0.0005 -0.0005 strain 0.5 0.5
```

## Description

This fix is an extension of the UEF (Uniform Extensional Flow) package
for LAMMPS. It implements NVE time integration under uniform extensional
flow fields using the SLLOD equations of motion with Kraynik-Reinelt
compatible boundary conditions.

Unlike the original [fix nvt/uef](fix_nh_uef.md) which uses Nose-Hoover
thermostatting, this fix performs NVE integration only. A separate
thermostat (e.g., [fix langevin](fix_langevin.md)) can be applied on
top. This makes it suitable for Langevin dynamics and dissipative
particle dynamics (DPD) simulations under extensional flow.

The applied flow field is set by the *erate* keyword. The values
*edot_x* and *edot_y* correspond to the strain rates in the xx and yy
directions. The strain rate in the zz direction is implicitly
$\dot{\varepsilon}_z = -(\dot{\varepsilon}_x + \dot{\varepsilon}_y)$ to
maintain a traceless flow field.

Common flow types:

- Uniaxial extension: `erate -e -e` (z-direction stretching)
- Biaxial extension: `erate e e` (xy-plane stretching)
- Planar extension: `erate e -e` (x stretching, y compression)

The *engrate* keyword provides an alternative way to specify uniaxial
extension using the engineering strain rate
$\dot{\varepsilon}_{\mathrm{eng}}$. The true strain rates are computed
as:

$$\dot{\varepsilon}_x = \dot{\varepsilon}_y
= -\frac{\dot{\varepsilon}_{\mathrm{eng}}}{2(1 + \varepsilon_{\mathrm{eng}})}$$

where $\varepsilon_{\mathrm{eng}}$ is the cumulative engineering strain.

The boundary conditions use the same lattice reduction approach as the
UEF package, based on `(Nicholson and Rutledge) <Nicholson_uefex>`.

------------------------------------------------------------------------

This fix automatically creates the following computes when it is
defined:

``` LAMMPS
compute {fix-ID}_temp group-ID temp/uefex
compute {fix-ID}_press group-ID pressure/uefex {fix-ID}_temp
```

See [compute temp/uefex](compute_temp_uefex.md) and [compute
pressure/uefex](compute_pressure_uefex.md) for details.

## Restart, fix_modify, output, run start/stop, minimize info

This fix writes the cumulative strain to `binary restart files
<restart>`.

<div class="note">

<div class="title">

Note

</div>

It is not necessary to set the *strain* keyword when resuming from a
restart file. The *strain* keyword is only needed when resuming from a
data file.

</div>

This fix supports the following [fix_modify](fix_modify.md) options to
change the strain rate during a simulation:

``` LAMMPS
fix_modify uefex u 0.001
fix_modify uefex b 0.001
fix_modify uefex p 0.001
fix_modify uefex ui 0.0005 0.001 100000
```

- *u* value = edot : switch to uniaxial extension with strain rate edot.
  Sets erate = (-edot/2, -edot/2), so the z-direction stretches at rate
  edot.
- *b* value = edot : switch to biaxial extension with strain rate edot.
  Sets erate = (edot, edot), so the xy-plane stretches at rate edot.
- *p* value = edot : switch to planar extension with strain rate edot.
  Sets erate = (-edot, 0), so x compresses and z stretches.
- *ui* values = irate frate nsteps : uniaxial extension with linearly
  ramped strain rate. The rate changes from *irate* to *frate* over
  *nsteps* steps.

These options allow changing the flow type or strain rate without
restarting the simulation.

This fix computes a global scalar (the kinetic energy) and a global
vector.

This fix is not invoked during [energy minimization](minimize.md).

## Restrictions

This fix is part of the UEFEX package. It is only enabled if LAMMPS was
built with that package.

The simulation box must be triclinic. If the box is initially
orthogonal, convert it before invoking the fix:

``` LAMMPS
change_box all triclinic
```

When the *strain* keyword is set to zero (or unset), the initial
simulation box must be cubic.

This fix requires the modified `domain.cpp` included in the UEFEX
package.

## Related commands

[fix nvt/uef](fix_nh_uef.md), [fix langevin](fix_langevin.md), [fix
barostat/uefex](fix_barostat_uefex.md), [compute
temp/uefex](compute_temp_uefex.md), [compute
pressure/uefex](compute_pressure_uefex.md), [compute
rotation/uefex](compute_rotation_uefex.md)

## Default

The default for *strain* is 0 0.

------------------------------------------------------------------------

<div id="Nicholson_uefex">

**(Nicholson and Rutledge)** Nicholson and Rutledge, J Chem Phys, 145,
244903 (2016).

</div>

<div id="MurashimaHagitaKawakatsu18_uefex">

**(Murashima, Hagita, and Kawakatsu)** Murashima, Hagita, and Kawakatsu,
Nihon Reoroji Gakkaishi, 46, 207-220 (2018).

</div>
