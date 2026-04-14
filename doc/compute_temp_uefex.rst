.. index:: compute temp/uefex

compute temp/uefex command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/uefex

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/uefex = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute mytemp all temp/uefex

Description
"""""""""""

This compute calculates the temperature by subtracting the streaming
velocity from each atom's velocity before computing the kinetic energy.
The streaming velocity is determined from the current deformation rate
of the simulation box (``h_rate``) in the upper-triangular (UT) LAMMPS frame.

The kinetic energy tensor (6-component vector output) is rotated from
the UT frame to the laboratory (LAB) frame using the rotation matrix Q
obtained from :doc:`fix nve/uefex <fix_nve_uefex>`.  This ensures that
the stress tensor components correspond to the physical flow directions.

The scalar temperature is frame-invariant and can also be obtained from
a standard :doc:`compute temp <compute_temp>`.

This compute is automatically created by
:doc:`fix nve/uefex <fix_nve_uefex>` with the ID ``{fix-ID}_temp``.

This compute can be used as a temperature bias for thermostats:

.. code-block:: LAMMPS

   fix lang all langevin 1.0 1.0 2.0 12345
   fix_modify lang temp uefex_temp

Output info
"""""""""""

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (the kinetic energy tensor in Voigt notation:
xx, yy, zz, xy, xz, yz) in the LAB frame.

Restrictions
""""""""""""

This compute is part of the UEFEX package.  It is only enabled if LAMMPS
was built with that package.

This compute can only be used when :doc:`fix nve/uefex <fix_nve_uefex>`
is active.

Related commands
""""""""""""""""

:doc:`compute temp <compute_temp>`,
:doc:`compute temp/uef <compute_temp_uef>`,
:doc:`fix nve/uefex <fix_nve_uefex>`,
:doc:`compute pressure/uefex <compute_pressure_uefex>`

Default
"""""""

none
