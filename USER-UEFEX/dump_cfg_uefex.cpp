/* ----------------------------------------------------------------------
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


#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "dump_cfg.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "error.h"
#include "dump_cfg_uef.h"
#include "dump_cfg_uefex.h"
#include "fix_nve_uefex.h"

using namespace LAMMPS_NS;

enum{INT,DOUBLE,STRING,BIGINT};   // same as in DumpCustom

#define UNWRAPEXPAND 10.0
#define ONEFIELD 32
#define DELTA 1048576

/* ---------------------------------------------------------------------- 
 * base method is mostly fine, just need to find the FixNHUef 
 * ----------------------------------------------------------------------*/

DumpCFGUefex::DumpCFGUefex(LAMMPS *lmp, int narg, char **arg) :
  DumpCFGUef(lmp, narg, arg)
{}

void DumpCFGUefex::init_style()
{
  DumpCFG::init_style();

  int i=0;
  for (i=0; i<modify->nfix; i++)
  {
    if (strcmp(modify->fix[i]->style,"nve/uefex")==0)
      break;
  }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use dump cfg/uefex without defining a fix nve/uefex");
  ifix_uefex=i;
  ifix_uef=i;
}


void DumpCFGUefex::write_header(bigint n)
{
  double rot[3][3];
  ((FixNVEUefex*) modify->fix[ifix_uefex])->get_rot(rot);
  fprintf(fp,"# === Rotation Matrix === \n");
  fprintf(fp,"# rot(0,0) = %g \n",rot[0][0]);
  fprintf(fp,"# rot(0,1) = %g \n",rot[0][1]);
  fprintf(fp,"# rot(0,2) = %g \n",rot[0][2]);
  fprintf(fp,"# rot(1,0) = %g \n",rot[1][0]);
  fprintf(fp,"# rot(1,1) = %g \n",rot[1][1]);
  fprintf(fp,"# rot(1,2) = %g \n",rot[1][2]);
  fprintf(fp,"# rot(2,0) = %g \n",rot[2][0]);
  fprintf(fp,"# rot(2,1) = %g \n",rot[2][1]);
  fprintf(fp,"# rot(2,2) = %g \n",rot[2][2]);

  DumpCFGUef::write_header(n);

}
