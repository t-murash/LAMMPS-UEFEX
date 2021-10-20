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

ComputeStyle(temp/uefex,ComputeTempUefex)

#else

#ifndef LMP_COMPUTE_TEMP_UEFEX_H
#define LMP_COMPUTE_TEMP_UEFEX_H

#include "compute_temp.h"
#include "compute_temp_uef.h"

namespace LAMMPS_NS {

  class ComputeTempUefex : public ComputeTempUef {
  public:
    ComputeTempUefex(class LAMMPS *, int, char **);
    virtual ~ComputeTempUefex(){}
    virtual void init();
    virtual double compute_scalar(); // Murashima 2018/12/25
    virtual void compute_vector(); // Murashima 2018/12/25
    void remove_bias(int i, double *v); // Murashima 2019/01/02
    void restore_bias(int i, double *v); // Murashima 2019/01/02
    
  };


}

#endif
#endif

/* ERROR/WARNING messages:

   This class inherits most of the warnings from ComputePressure. The
   only addition is:

   E: Can't use compute temp/uefex without defining a fix nve/uefex

   Self-explanatory.

*/
