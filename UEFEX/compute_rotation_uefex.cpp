/*
  ----------------------------------------------------------------------
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
  https://www.lammps.com 
  Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

  Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.  This software is distributed under
  the GNU General Public License.

  See the README file in the top-level LAMMPS directory.

  Contributing author: T. Murashima (Tohoku Univ, JPN)
  -------------------------------------------------------------------------
*/

#include "mpi.h"
#include <cstring>
#include <cstdlib>
#include "compute_rotation_uefex.h"
#include "fix_nve_uefex.h"
#include "fix_npt_uef.h"
#include "fix_nvt_uef.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/*
  ----------------------------------------------------------------------
  Default values for the ext flags
  ----------------------------------------------------------------------
*/
ComputeRotationUefex::ComputeRotationUefex(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  vector_flag = 1;
  extvector = 1;
  size_vector = 9;
  timeflag = 1;
  peflag = 1;
  vector=new double[9];
}

ComputeRotationUefex::~ComputeRotationUefex()
{
  delete [] vector;
}

/*
  ----------------------------------------------------------------------
  Check for the uefex fix
  ----------------------------------------------------------------------
*/
void ComputeRotationUefex::init()
{
  int i=0;
  for (i=0; i<modify->nfix; i++)
    {
      if (strcmp(modify->fix[i]->style,"nve/uefex")==0){
	uef_flag=0;
	break;
      }
      if (strcmp(modify->fix[i]->style,"nvt/uef")==0){
	uef_flag=1;
	break;
      }
      if (strcmp(modify->fix[i]->style,"npt/uef")==0){
	uef_flag=2;
	break;
      }
    }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use compute rotation/uefex without defining fix nve/uefex or fix nvt/uef or fix npt/uef");
  
  ifix_uef=i;
}

/* ----------------------------------------------------------------------
   Get the rotation tensor
   ------------------------------------------------------------------------- */

void ComputeRotationUefex::compute_vector()
{
  //  invoked_vector = update->ntimestep;
  //  if(update->vflag_global != invoked_vector)
  //    error->all(FLERR,"Rotation was not tallied on needed timestep");

  if(uef_flag==0){
    ((FixNVEUefex*) modify->fix[ifix_uef])->get_rot(rot);
  }else if(uef_flag==1){
    ((FixNVTUef*) modify->fix[ifix_uef])->get_rot(rot);
  }else if(uef_flag==2){
    ((FixNPTUef*) modify->fix[ifix_uef])->get_rot(rot);
  }
  
  vector[0]=rot[0][0];
  vector[1]=rot[0][1];
  vector[2]=rot[0][2];
  vector[3]=rot[1][0];
  vector[4]=rot[1][1];
  vector[5]=rot[1][2];
  vector[6]=rot[2][0];
  vector[7]=rot[2][1];
  vector[8]=rot[2][2];
  //for(int k=0;k<9;k++)printf("test[%d]== %lf ",k,vector[k]);
}
