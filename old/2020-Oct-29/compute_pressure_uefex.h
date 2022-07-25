/*
  ----------------------------------------------------------------------
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

ComputeStyle(pressure/uefex,ComputePressureUefex)

#else

#ifndef LMP_COMPUTE_PRESSURE_UEFEX_H
#define LMP_COMPUTE_PRESSURE_UEFEX_H

#include "compute_pressure.h"
#include "compute_pressure_uef.h"

namespace LAMMPS_NS {

  class ComputePressureUefex : public ComputePressureUef {
  public:
    ComputePressureUefex(class LAMMPS *, int, char **);
    virtual ~ComputePressureUefex(){}
    virtual void init();
    virtual double compute_scalar();
    virtual void compute_vector();
  };


}

#endif
#endif

/*
   ERROR/WARNING messages:

   This class inherits most of the warnings from ComputePressure. The
   only addition is:

   E: Can't use compute pressure/uefex without defining a fix nve/uefex

   Self-explanatory.  

   W: The temperature used in compute pressure/uefex is not of style temp/uefex

   Self-explanatory.

*/
