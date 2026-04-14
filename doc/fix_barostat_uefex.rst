.. index:: fix barostat/uefex

fix barostat/uefex command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID barostat/uefex keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* barostat/uefex = style name of this fix command
* one or more keyword/value pairs must be appended

  .. parsed-literal::

     keyword = *iso* or *const_xy* or *const_z* or *modulus*
       *iso* values = Pstart Pstop Pdamp
         Pstart, Pstop = desired pressure at start/end of run (pressure units)
         Pdamp = damping parameter (time units)
       *const_xy* values = Pstart Pstop Pdamp
         control pressure in xy-plane only
       *const_z* values = Pstart Pstop Pdamp
         control pressure in z-direction only
       *modulus* value = B
         B = bulk modulus for Berendsen scaling (pressure units, default 10.0)

Examples
""""""""

.. code-block:: LAMMPS

   fix baro all barostat/uefex iso 0.0 0.0 1.0
   fix baro all barostat/uefex const_xy 0.0 0.0 1.0
   fix baro all barostat/uefex iso 0.0 0.0 1.0 modulus 20.0

Description
"""""""""""

This fix implements a Berendsen-style barostat for use with
:doc:`fix nve/uefex <fix_nve_uefex>` under extensional flow.  It rescales
the simulation box isotropically (or anisotropically) to achieve the
target pressure.

The *iso* keyword controls the hydrostatic pressure (Pxx+Pyy+Pzz)/3
by scaling all three box dimensions equally.

The *const_xy* keyword controls the pressure in the xy-plane only
(Pxx+Pyy)/2, leaving the z-dimension uncontrolled.  This is useful
for uniaxial extension where the transverse directions should be at
a target pressure (e.g., free surface conditions).

The *const_z* keyword controls the pressure in the z-direction only,
leaving x and y uncontrolled.

The *modulus* keyword sets the bulk modulus used in the Berendsen
scaling formula.  The default value is 10.0 in pressure units.

This fix automatically creates the following computes:

.. code-block:: LAMMPS

   compute {fix-ID}_temp all temp/uefex
   compute {fix-ID}_press all pressure/uefex {fix-ID}_temp

Restrictions
""""""""""""

This fix is part of the UEFEX package.  It is only enabled if LAMMPS
was built with that package.

This fix must be used in conjunction with
:doc:`fix nve/uefex <fix_nve_uefex>`.

Related commands
""""""""""""""""

:doc:`fix nve/uefex <fix_nve_uefex>`,
:doc:`fix press/berendsen <fix_press_berendsen>`

Default
"""""""

The default for *modulus* is 10.0.
