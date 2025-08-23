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
*/
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include "compute_temp.h"
#include "compute_temp_uef.h"
#include "compute_temp_uefex.h"
#include "fix_nve_uefex.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "error.h"
#include "domain.h"


#include "atom.h"
#include "force.h"
#include "comm.h"
#include "group.h"

using namespace LAMMPS_NS;

/*
  ----------------------------------------------------------------------
  Base constructor is fine
  ----------------------------------------------------------------------
*/
ComputeTempUefex::ComputeTempUefex(LAMMPS *lmp, int narg, char **arg) :
  ComputeTempUef(lmp, narg, arg) 
{
  rot_flag=true;
  tempflag=1;
  tempbias=1;
}

/*
  ----------------------------------------------------------------------
  Check for the uef fix
  ----------------------------------------------------------------------
*/
void ComputeTempUefex::init()
{
  ComputeTemp::init();
  int i=0;
  for (i=0; i<modify->nfix; i++)
    {
      if (strcmp(modify->fix[i]->style,"nve/uefex")==0)
	break;
    }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use compute temp/uefex without defining a fix nve/uefex");
  ifix_uef=i;

}

// Murashima 2019/01/02
double ComputeTempUefex::compute_scalar()
{
  double lamda[3],vstream[3],vthermal[3];
  invoked_scalar = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h_rate=domain->h_rate;
  double *h_ratelo=domain->h_ratelo;


  double t = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];
      if (rmass)
        t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
              vthermal[2]*vthermal[2]) * rmass[i];
      else
        t += (vthermal[0]*vthermal[0] + vthermal[1]*vthermal[1] +
              vthermal[2]*vthermal[2]) * mass[type[i]];
    }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;

  return scalar;
}

// Murashima 2019/01/02
void ComputeTempUefex::compute_vector()
{
  double lamda[3],vstream[3],vthermal[3];

  invoked_vector = update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;



  //Note: The thermal stress is obtained here on the LAB-frame!


  double massone,t[6];
  for (int i = 0; i < 6; i++) t[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->x2lamda(x[i],lamda);
      vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + h_rate[4]*lamda[2] + h_ratelo[0];
      vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
      vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      vthermal[0] = v[i][0] - vstream[0];
      vthermal[1] = v[i][1] - vstream[1];
      vthermal[2] = v[i][2] - vstream[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * vthermal[0]*vthermal[0];
      t[1] += massone * vthermal[1]*vthermal[1];
      t[2] += massone * vthermal[2]*vthermal[2];
      t[3] += massone * vthermal[0]*vthermal[1];
      t[4] += massone * vthermal[0]*vthermal[2];
      t[5] += massone * vthermal[1]*vthermal[2];
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;


  if (rot_flag)
  {
    double rot[3][3];
    (dynamic_cast<FixNVEUefex*>( modify->fix[ifix_uef]))->get_rot(rot);
    virial_rot(vector,rot);
  }
  // rot_flag, virial_rot are found in compute_temp_uef
}


// Murashima 2019/01/02
// remove_bias & restore_bias are modified from those found in compute_temp_deform.cpp.
void ComputeTempUefex::remove_bias(int i, double *v)
{
  double lamda[3];
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;
  double **x = atom->x;


  domain->x2lamda(atom->x[i],lamda);
  // vbias is stored in compute.h

  vbias[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] + h_rate[4]*lamda[2] + h_ratelo[0];
  vbias[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] + h_ratelo[1];
  vbias[2] = h_rate[2]*lamda[2] + h_ratelo[2];

  v[0] -= vbias[0];
  v[1] -= vbias[1];
  v[2] -= vbias[2];
}
void ComputeTempUefex::restore_bias(int i, double *v)
{
  v[0] += vbias[0];
  v[1] += vbias[1];
  v[2] += vbias[2];
}
