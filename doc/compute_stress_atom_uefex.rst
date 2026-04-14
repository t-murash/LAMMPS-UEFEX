.. index:: compute stress/atom/uefex

compute stress/atom/uefex command
=================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID stress/atom/uefex temp-ID keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* stress/atom/uefex = style name of this compute command
* temp-ID = ID of a temperature compute for kinetic energy contribution, or NULL
* zero or more keywords may be appended
* keyword = *ke* or *pair* or *bond* or *angle* or *dihedral* or *improper* or *kspace* or *fix* or *virial*

Examples
""""""""

.. code-block:: LAMMPS

   compute mystress all stress/atom/uefex uefex_temp
   compute mystress all stress/atom/uefex NULL virial

Description
"""""""""""

This compute calculates the per-atom stress tensor in the laboratory
(LAB) frame during UEFEX extensional flow simulations.  It is analogous
to :doc:`compute stress/atom <compute_stress_atom>`, but the stress
tensor is rotated from the upper-triangular (UT) LAMMPS frame to the
LAB frame using the rotation matrix Q from
:doc:`fix nve/uefex <fix_nve_uefex>`.

The per-atom stress tensor has 6 components in Voigt notation:
xx, yy, zz, xy, xz, yz, where the directions correspond to the
LAB frame.

The keywords are the same as for
:doc:`compute stress/atom <compute_stress_atom>`.

Output info
"""""""""""

This compute calculates a per-atom array with 6 columns.  The values
are in stress*volume units.  See :doc:`compute stress/atom <compute_stress_atom>`
for details on units and usage.

Restrictions
""""""""""""

This compute is part of the UEFEX package.  It is only enabled if LAMMPS
was built with that package.

This compute can only be used when :doc:`fix nve/uefex <fix_nve_uefex>`
is active.

Related commands
""""""""""""""""

:doc:`compute stress/atom <compute_stress_atom>`,
:doc:`fix nve/uefex <fix_nve_uefex>`,
:doc:`compute pressure/uefex <compute_pressure_uefex>`

Default
"""""""

none
