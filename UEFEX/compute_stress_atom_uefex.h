/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: T. Murashima (Tohoku Univ, JPN)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(stress/atom/uefex,ComputeStressAtomUefex);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRESS_ATOM_UEFEX_H
#define LMP_COMPUTE_STRESS_ATOM_UEFEX_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressAtomUefex : public Compute {
 public:
  ComputeStressAtomUefex(class LAMMPS *, int, char **);
  ~ComputeStressAtomUefex();
  void init();
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void virial_rot(double *x, const double r[3][3]);

 private:
  int keflag, pairflag, bondflag, angleflag, dihedralflag, improperflag;
  int kspaceflag, fixflag, biasflag;
  Compute *temperature;
  char *id_temp;

  int nmax;
  double **stress;

  // 2022/04/15 by TM
  int ifix_uef;
  int uef_flag;
  double rot[3][3];
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute stress/atom/uefex temperature ID

Self-explanatory.

E: Compute stress/atom/uefex temperature ID does not compute temperature

The specified compute must compute temperature.

E: Per-atom virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to have
tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
