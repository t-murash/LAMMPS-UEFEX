.. index:: compute rotation/uefex

compute rotation/uefex command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID rotation/uefex

* ID, group-ID are documented in :doc:`compute <compute>` command
* rotation/uefex = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute myrot all rotation/uefex
   fix rotout all ave/time 5000 1 5000 c_myrot[1] c_myrot[2] c_myrot[3] &
       c_myrot[4] c_myrot[5] c_myrot[6] c_myrot[7] c_myrot[8] c_myrot[9] &
       file rotation_uefex.txt

Description
"""""""""""

This compute returns the 3x3 rotation matrix Q that relates the
upper-triangular (UT) LAMMPS frame to the laboratory (LAB) frame
during UEFEX extensional flow simulations.

The rotation matrix satisfies:

.. math::

   Q \cdot L_{\mathrm{LAB}} = L_{\mathrm{UT}}

where :math:`L_{\mathrm{LAB}}` is the lattice in the LAB frame and
:math:`L_{\mathrm{UT}}` is the upper-triangular lattice used by LAMMPS.
To convert from UT to LAB frame: :math:`L_{\mathrm{LAB}} = Q^T \cdot L_{\mathrm{UT}}`.

This compute is primarily used for post-processing and visualization.

Output info
"""""""""""

This compute calculates a global vector of length 9 containing the
rotation matrix elements in row-major order:
Q[0][0], Q[0][1], Q[0][2], Q[1][0], Q[1][1], Q[1][2], Q[2][0], Q[2][1], Q[2][2].

Restrictions
""""""""""""

This compute is part of the UEFEX package.  It is only enabled if LAMMPS
was built with that package.

This compute can only be used when :doc:`fix nve/uefex <fix_nve_uefex>`
is active.

Related commands
""""""""""""""""

:doc:`fix nve/uefex <fix_nve_uefex>`,
:doc:`fix ave/time <fix_ave_time>`

Default
"""""""

none
