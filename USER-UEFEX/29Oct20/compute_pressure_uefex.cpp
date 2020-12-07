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

#include "mpi.h"
#include <cstring>
#include <cstdlib>
#include "compute_pressure.h"
#include "compute_pressure_uef.h"
#include "compute_pressure_uefex.h"
#include "fix_nve_uefex.h"
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
ComputePressureUefex::ComputePressureUefex(LAMMPS *lmp, int narg, char **arg) :
  ComputePressureUef(lmp, narg, arg)
{
  ext_flags[0] = true;
  ext_flags[1] = true;
  ext_flags[2] = true;
  in_fix=false;
}

/*
  ----------------------------------------------------------------------
  Check for the uefex fix
  ----------------------------------------------------------------------
*/
void ComputePressureUefex::init()
{
  ComputePressure::init();
  int i=0;
  for (i=0; i<modify->nfix; i++)
    {
      if (strcmp(modify->fix[i]->style,"nve/uefex")==0)
	break;
    }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use compute pressure/uefex without defining a fix nve/uefex");
  
  ifix_uef=i;
  ((FixNVEUefex*) modify->fix[ifix_uef])->get_ext_flags(ext_flags);

  if (strcmp(temperature->style,"temp/uefex") != 0)
    error->warning(FLERR,"The temperature used in compute pressure/uefex is not of style temp/uefex");

}

/*
  ----------------------------------------------------------------------
  Compute pressure in the directions i corresponding to ext_flag[i]=true
  ----------------------------------------------------------------------
*/
double ComputePressureUefex::compute_scalar()
{

  temperature->compute_scalar(); // Streaming velocity is considered in compute_temp_uefex.cpp

  // if all pressures are external the scalar is found as normal
  if (ext_flags[0] && ext_flags[1] && ext_flags[2])
    return ComputePressure::compute_scalar();
  /*
    Note: Pressure (scalar) is rotationally invariant.
  */

  // otherwise compute the full tensor and average desired components
  compute_vector();// Pressure (tensor)
  addstep(update->ntimestep+1);

  int k =0;
  scalar = 0;
  if (ext_flags[0])
    {
      scalar += vector[0];
      k++;
    }
  if (ext_flags[1])
    {
      scalar += vector[1];
      k++;
    }
  if (ext_flags[2])
    {
      scalar += vector[2];
      k++;
    }

  scalar /= k;
  return scalar;
}

/* ----------------------------------------------------------------------
   Compute the pressure tensor in the rotated coordinate system
   ------------------------------------------------------------------------- */
void ComputePressureUefex::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  if (force->kspace && kspace_virial && force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' for "
	       "tensor components with kspace_style msm");

  // invoke temperature if it hasn't been already

  double *ke_tensor;

  //printf("keflag %d\n",keflag);
  
  
  if (keflag) {
    if (temperature->invoked_vector != update->ntimestep)temperature->compute_vector();
    ke_tensor = temperature->vector; // temp/uefex obtains thermal tensor of LAB frame value.
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6,3);
    /*
      The above virial_compute is defined in compute_pressure.cpp
      This function obtains LAMMPS(UT) frame value.
    */

    double r[3][3];
    ( (FixNVEUefex*) modify->fix[ifix_uef])->get_rot(r);
    virial_rot(virial,r); // LAMMPS(UT) to LAB, virial_rot is defined in compute_pressure_uef.cpp
    /*
      Now, virial is LAB frame value
     */

    if (keflag) {
      for (int i = 0; i < 6; i++)
	{
	  vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
	  //printf("compute_pressure_uefex.cpp %lf %lf\n",ke_tensor[0],virial[0]);
	  
	}
      
    } else
      for (int i = 0; i < 6; i++)
        vector[i] = virial[i] * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4,2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
  // Note: virial is defined as
  // [ 0 3 4 ]
  // [ 3 1 5 ]
  // [ 4 5 2 ]
  // This definition is different from the box shape tensor.
  // [ 0 5 4 ]
  // [ 5 1 3 ]
  // [ 4 3 2 ]
  
}
