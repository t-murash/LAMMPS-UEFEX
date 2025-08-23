/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: T. Murashima (Tohoku Univ, JPN)
------------------------------------------------------------------------- */


#ifdef FIX_CLASS
// clang-format off
FixStyle(barostat/uefex,FixBarostatUefex);
// clang-format on
#else

#ifndef LMP_FIX_BAROSTAT_UEFEX_H
#define LMP_FIX_BAROSTAT_UEFEX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBarostatUefex : public Fix {
 public:
  FixBarostatUefex(class LAMMPS *, int, char **);
  ~FixBarostatUefex() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  int modify_param(int, char **) override;

 protected:
  int dimension, which;
  double bulkmodulus;

  int pstyle, pcouple, allremap;
  int p_flag[3];    // 1 if control P on this dim, 0 if not
  double p_start[3], p_stop[3];
  double p_period[3], p_target[3];
  double p_current[3], dilation;
  double factor[3];
  int kspace_flag;    // 1 if KSpace invoked, 0 if not
  int nrigid;         // number of rigid fixes
  int *rfix;          // indices of rigid fixes

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tflag, pflag;

  void couple();
  void remap();
};

}    // namespace LAMMPS_NS

#endif
#endif
