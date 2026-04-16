/*
  ----------------------------------------------------------------------
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
  https://www.lammps.org
  Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

  Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.  This software is distributed under
  the GNU General Public License.

  See the README file in the top-level LAMMPS directory.

  Contributing author: T. Murashima (Tohoku Univ, JPN)
  -------------------------------------------------------------------------
  atom->v is treated as peculiar velocity. No streaming subtraction.
  The KE tensor is rotated to the LAB frame using the UEF rotation matrix.
*/
#include <cstring>
#include "compute_temp_uefex.h"
#include "fix_nve_uefex.h"
#include "modify.h"
#include "fix.h"
#include "error.h"

using namespace LAMMPS_NS;

ComputeTempUefex::ComputeTempUefex(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempUef(lmp, narg, arg)
{
  rot_flag = true;
}

void ComputeTempUefex::init()
{
  ComputeTemp::init();
  int i = 0;
  for (i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "nve/uefex") == 0) break;
  }
  if (i == modify->nfix)
    error->all(FLERR,
               "Can't use compute temp/uefex without defining a fix nve/uefex");
  ifix_uef = i;
}

void ComputeTempUefex::compute_vector()
{
  ComputeTemp::compute_vector();
  if (rot_flag) {
    double rot[3][3];
    (dynamic_cast<FixNVEUefex*>(modify->fix[ifix_uef]))->get_rot(rot);
    virial_rot(vector, rot);
  }
}
