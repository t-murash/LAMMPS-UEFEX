/* -*- c++ -*- ----------------------------------------------------------
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
  http://lammps.sandia.gov, Sandia National Laboratories
  Steve Plimpton, sjplimp@sandia.gov

  Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.  This software is distributed under
  the GNU General Public License.

  See the README file in the top-level LAMMPS directory.

  Contributing author: T. Murashima (Tohoku Univ, JPN)
  -------------------------------------------------------------------------
*/

#ifdef COMPUTE_CLASS

ComputeStyle(rotation/uefex,ComputeRotationUefex)

#else

#ifndef LMP_COMPUTE_ROTATION_UEFEX_H
#define LMP_COMPUTE_ROTATION_UEFEX_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRotationUefex : public Compute {
 public:
  ComputeRotationUefex(class LAMMPS *, int, char **);
  virtual ~ComputeRotationUefex();
  virtual void init();
  virtual void compute_vector();

 protected:
  int ifix_uef;
  double rot[3][3];
  int uef_flag;// 0:nve/uefex, 1:nvt/uef, 2:npt/uef
};

}

#endif
#endif

/*
   ERROR/WARNING messages:

   This class inherits most of the warnings from ComputePressure. The
   only addition is:

   E: Can't use compute rotation/uefex without defining fix nve/uefex or fix nvt/uef or fix npt/uef

   Self-explanatory.  

   W: The temperature used in compute rotation/uefex is not of style temp/uefex

   Self-explanatory.

*/
