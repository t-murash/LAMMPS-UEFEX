/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: T. Murashima (Tohoku Univ, JPN)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(cfg/uefex,DumpCFGUefex)

#else

#ifndef LMP_DUMP_CFG_UEFEX_H
#define LMP_DUMP_CFG_UEFEX_H

#include "dump_cfg_uef.h"

namespace LAMMPS_NS {

class DumpCFGUefex : public DumpCFGUef {
 public:
  DumpCFGUefex(LAMMPS*, int, char**);
  virtual void init_style();
  virtual void write_header(bigint);
 protected:
  int ifix_uefex;
};

}

#endif
#endif

/* ERROR/WARNING messages:
  
E: Can't use dump cfg/uefex without defining a fix nve/uefex

Self-explanatory.

*/
