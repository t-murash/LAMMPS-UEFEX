.. index:: compute pressure/uefex

compute pressure/uefex command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID pressure/uefex temp-ID keyword ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* pressure/uefex = style name of this compute command
* temp-ID = ID of a temp/uefex compute, or NULL
* zero or more keywords may be appended
* keyword = *ke* or *pair* or *bond* or *angle* or *dihedral* or *improper* or *kspace* or *fix* or *virial*

Examples
""""""""

.. code-block:: LAMMPS

   compute mypress all pressure/uefex uefex_temp
   compute mypress all pressure/uefex uefex_temp virial

Description
"""""""""""

This compute calculates the pressure tensor in the laboratory (LAB)
frame during UEFEX extensional flow simulations.  The virial
contributions and kinetic energy tensor are rotated from the
upper-triangular (UT) LAMMPS frame to the LAB frame using the rotation
matrix Q from :doc:`fix nve/uefex <fix_nve_uefex>`.

The 6-component pressure tensor is returned in Voigt notation:
Pxx, Pyy, Pzz, Pxy, Pxz, Pyz, where the directions correspond to the
LAB frame.

For uniaxial extension (``erate -e -e``), the extensional stress is:

.. math::

   N_1 = \sigma_{zz} - \frac{\sigma_{xx} + \sigma_{yy}}{2}
       = -P_{zz} + \frac{P_{xx} + P_{yy}}{2}

This compute is automatically created by
:doc:`fix nve/uefex <fix_nve_uefex>` with the ID ``{fix-ID}_press``.

The keywords and general output information are the same as for
:doc:`compute pressure <compute_pressure>`.

Restrictions
""""""""""""

This compute is part of the UEFEX package.  It is only enabled if LAMMPS
was built with that package.

This compute can only be used when :doc:`fix nve/uefex <fix_nve_uefex>`
is active.

The kinetic contribution to the pressure tensor will be accurate only
when the compute specified by *temp-ID* is a
:doc:`compute temp/uefex <compute_temp_uefex>`.

Related commands
""""""""""""""""

:doc:`compute pressure <compute_pressure>`,
:doc:`compute pressure/uef <compute_pressure_uef>`,
:doc:`fix nve/uefex <fix_nve_uefex>`,
:doc:`compute temp/uefex <compute_temp_uefex>`

Default
"""""""

none
