/*
  ----------------------------------------------------------------------
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
  https://www.lammps.org/
  Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

  Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.  This software is distributed under
  the GNU General Public License.

  See the README file in the top-level LAMMPS directory.

  Contributing author: T. Murashima (Tohoku Univ, JPN)
  -------------------------------------------------------------------------
*/

#ifdef COMPUTE_CLASS

ComputeStyle(temp/uefex,ComputeTempUefex);

#else

#ifndef LMP_COMPUTE_TEMP_UEFEX_H
#define LMP_COMPUTE_TEMP_UEFEX_H

#include "compute_temp_uef.h"

namespace LAMMPS_NS {

class ComputeTempUefex : public ComputeTempUef {
 public:
  ComputeTempUefex(class LAMMPS *, int, char **);
  void init() override;
  void compute_vector() override;
};

}

#endif
#endif
