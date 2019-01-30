/*
  ----------------------------------------------------------------------
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
  www.cs.sandia.gov/~sjplimp/lammps.html
  Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

  Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.  This software is distributed under
  the GNU General Public License.

  See the README file in the top-level LAMMPS directory.

  Contributing author: T. Murashima (Tohoku Univ, JPN)
  -------------------------------------------------------------------------
*/

#include <cstring>
#include <cstdlib>
#include <cmath>
#include "fix_nve_uefex.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "comm.h"
#include "citeme.h"
#include "irregular.h"
#include "modify.h"
#include "compute.h"
#include "kspace.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "output.h"
#include "timer.h"
#include "neighbor.h"
#include "compute_pressure_uefex.h"
#include "compute_temp_uefex.h"
#include "uef_utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{ISO,ANISO,TRICLINIC};
enum{NO_REMAP, X_REMAP, V_REMAP}; // same as fix_deform.cpp

// citation info

static const char cite_user_uef_package[] =
  "USER-UEF package:\n\n"
  "@Article{NicholsonRutledge16,\n"
  "author = {David A. Nicholson and Gregory C. Rutledge},\n"
  "title = {Molecular simulation of flow-enhanced nucleation in n-eicosane melt\
s under steady shear and uniaxial extension},\n"
  "journal = {The Journal of Chemical Physics},\n"
  "volume = {145},\n"
  "number = {24},\n"
  "pages = {244903},\n"
  "year = {2016}\n"
  "}\n\n";

static const char cite_user_uefex_package[] =
  "USER-UEFEX package:\n\n"
  "@Article{MurashimaHagitaKawakatsu18,\n"
  "author = {Takahiro Murashima, Katsumi Hagita, and Toshihiro Kawakatsu},\n"
  "title = {Elongational Viscosity of Weakly Entangled Polymer Melt via Coarse-Grained Molecular Dynamics Simulation},\n"
  "journal = {Nihon Reoroji Gakkaishi (The Journal of Society of Rheology, Japan)},\n"
  "volume = {46},\n"
  "number = {5},\n"
  "pages = {207-220},\n"
  "year = {2018}\n"
  "}\n\n";


FixNVEUefex::FixNVEUefex(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  rfix(NULL),id_dilate(NULL),irregular(NULL),id_temp(NULL),id_press(NULL)
{

  if(lmp->citeme){
    lmp->citeme->add(cite_user_uef_package);
    lmp->citeme->add(cite_user_uefex_package);
  }

  // fix_nh.cpp
  restart_global = 1;
  dynamic_group_allow = 1;
  time_integrate = 1;
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 0;

  allremap = 1;
  id_dilate = NULL;

  //initialization
  erate[0] = erate[1] = 0;
  //default values 
  strain[0]=strain[1]= 0;
  ext_flags[0]=ext_flags[1]=ext_flags[2] = true;

  // parse remaining input
  bool erate_flag = false;
  int iarg = 3;
  while (iarg <narg)
    {
      if (strcmp(arg[iarg],"erate")==0) {
	if (iarg+3 > narg) error->all(FLERR,"Illegal fix nve/uefex command");
	erate[0] = force->numeric(FLERR,arg[iarg+1]);
	erate[1] = force->numeric(FLERR,arg[iarg+2]);
	erate_flag = true;
	iarg += 3;
      }
      else if (strcmp(arg[iarg],"strain")==0) {
	if (iarg+3 > narg) error->all(FLERR,"Illegal fix nve/uefex command");
	strain[0] = force->numeric(FLERR,arg[iarg+1]);
	strain[1] = force->numeric(FLERR,arg[iarg+2]);
	iarg += 3;
      }
      else
	error->all(FLERR,"Illegal fix nve/uefex command");

    }
  if (!erate_flag)
    error->all(FLERR,"Keyword erate must be set for fix nve/uefex command");


  if (!domain->triclinic)
    error->all(FLERR,"Simulation box must be triclinic for fix nve/uefex");

  //check for conditions that impose a deviatoric stress

  double erate_tmp[3];
  erate_tmp[0]=erate[0];
  erate_tmp[1]=erate[1];
  erate_tmp[2]=-erate[0]-erate[1];

  // conditions that produce a deviatoric stress have already
  // been eliminated.
  //deviatoric_flag=0;

  // need pre_exchange and irregular migration
  pre_exchange_flag = 1;
  irregular = new Irregular(lmp);

  // flag that I change the box here (in case of nvt)
  box_change_shape = 1;

  // initialize the UEFBox class which computes the box at each step
  uefbox = new UEF_utils::UEFBox();
  uefbox->set_strain(strain[0],strain[1]);

  // reset fixedpoint to the stagnation point. I don't allow fixedpoint 
  // to be set by the user.
  fixedpoint[0] = domain->boxlo[0];
  fixedpoint[1] = domain->boxlo[1];
  fixedpoint[2] = domain->boxlo[2];

  // Create temp and pressure computes for uefex

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");
  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp/uefex";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tcomputeflag = 1;

  n = strlen(id) + 7;
  id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_press");
  newarg = new char*[4];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure/uefex";
  newarg[3] = id_temp;
  modify->add_compute(4,newarg);
  delete [] newarg;
  pcomputeflag = 1;
  nevery = 1;


}

/* ----------------------------------------------------------------------
 * Erase the UEFBox object and get rid of the pressure compute if the nvt 
 * version is being used. Everything else will be done in base destructor
 * ---------------------------------------------------------------------- */
FixNVEUefex::~FixNVEUefex()
{
  delete [] id_dilate;
  delete [] rfix;
  delete irregular;
  delete uefbox;
  modify->delete_compute(id_temp);
  delete [] id_temp;
  modify->delete_compute(id_press);
  delete [] id_press;


}

/* ----------------------------------------------------------------------
 * Make the end_of_step() routine callable
 * ---------------------------------------------------------------------- */
int FixNVEUefex::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  if (pre_exchange_flag) mask |= PRE_EXCHANGE;
  mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
 * Error checking. Set the temperature & pressure pointers 
 * ---------------------------------------------------------------------- */
void FixNVEUefex::init()
{
  //FixNH::init();

  // recheck that dilate group has not been deleted

  if (allremap == 0) {
    int idilate = group->find(id_dilate);
    if (idilate == -1)
      error->all(FLERR,"Fix nvt/npt/nph dilate group ID does not exist");
    dilate_group_bit = group->bitmask[idilate];
  }

  int icompute = modify->find_compute(id_temp);
  temperature = modify->compute[icompute];
  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;
  icompute = modify->find_compute(id_press);
  pressure = modify->compute[icompute];
  

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dthalf = 0.5 * update->dt;
  dto = dthalf;

  // Murashiam 2018/12/31
  //dtv = update->dt; set in fix_nh.cpp
  dt_inv=1.0/dtv;
  remapflag=V_REMAP;

  if (strstr(update->integrate_style,"respa")) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
    dto = 0.5*step_respa[0];
  }

  delete [] rfix;
  nrigid = 0;
  rfix = NULL;

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->rigid_flag) rfix[nrigid++] = i;
  }




  // find conflict with fix/deform or other box chaging fixes
  for (int i=0; i < modify->nfix; i++)
    {
      if (strcmp(modify->fix[i]->id,id) != 0)
	if (modify->fix[i]->box_change_shape != 0)
	  error->all(FLERR,"Can't use another fix which changes box shape with fix nve/uefex");
    }


  
  if (strcmp(pressure->style,"pressure/uefex") != 0)
    error->all(FLERR,"Using fix nve/uefex without a compute pressure/uefex");


  if (strcmp(temperature->style,"temp/uefex") != 0)
    error->all(FLERR,"Using fix nve/uefex without a compute temp/uefex");

  
}


void FixNVEUefex::setup(int j)
{
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  double tol = 1e-4;
  // ensure the box is ok for uef
  bool isok = true;
  isok &= nearly_equal(domain->h[0],box[0][0],tol);
  isok &= nearly_equal(domain->h[1],box[1][1],tol);
  isok &= nearly_equal(domain->h[2],box[2][2],tol);
  isok &= nearly_equal(domain->xy,box[0][1],tol);
  isok &= nearly_equal(domain->yz,box[1][2],tol);
  isok &= nearly_equal(domain->xz,box[0][2],tol);
  if (!isok)
    error->all(FLERR,"Initial box is not close enough to the expected uef box");

  uefbox->get_rot(rot);
  ((ComputeTempUefex*) temperature)->yes_rot();
  ((ComputePressureUefex*) pressure)->in_fix = true;
  ((ComputePressureUefex*) pressure)->update_rot();

}

/* ----------------------------------------------------------------------
 * rotate -> initial integration step -> rotate back
 * ---------------------------------------------------------------------- */
void FixNVEUefex::initial_integrate(int vflag)
{
  inv_rotate_x(rot);
  inv_rotate_v(rot);
  inv_rotate_f(rot);
  ((ComputeTempUefex*) temperature)->no_rot();// LAB
  nve_v();
  nve_x();
  rotate_x(rot);
  rotate_v(rot);
  rotate_f(rot);
  ((ComputeTempUefex*) temperature)->yes_rot();// LAMMPS(UT)
}

/* ----------------------------------------------------------------------
 * rotate -> final integration step -> rotate back
 * ---------------------------------------------------------------------- */
void FixNVEUefex::final_integrate()
{
  // update rot here since it must directly follow the virial calculation
  ((ComputePressureUefex*) pressure)->update_rot();
  inv_rotate_v(rot);
  inv_rotate_f(rot);
  ((ComputeTempUefex*) temperature)->no_rot();// LAB
  nve_v();
  rotate_v(rot);
  rotate_f(rot);
  ((ComputeTempUefex*) temperature)->yes_rot();// LAMMPS(UT)
}


/* ----------------------------------------------------------------------
 * rotate -> initial integration step -> rotate back (RESPA)
 * ---------------------------------------------------------------------- */
void FixNVEUefex::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  if(ilevel==0)initial_integrate(vflag);
  else final_integrate();
}

/* ----------------------------------------------------------------------
 * at outer level: call this->final_integrate()
 * at other levels: rotate -> 2nd verlet step -> rotate back
 * ---------------------------------------------------------------------- */
void FixNVEUefex::final_integrate_respa(int ilevel, int iloop)
{
  // set timesteps by level
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];
  // outermost level - update eta_dot and omega_dot, apply via final_integrate
  // all other levels - NVE update of v
  if (ilevel == nlevels_respa-1) final_integrate();
  else 
    {
      inv_rotate_v(rot);
      inv_rotate_f(rot);
      nve_v(); // LAB frame
      rotate_v(rot);
      rotate_f(rot);
    }
}

/* ----------------------------------------------------------------------
   SLLOD velocity update in time-reversible (by Nicholson) increments 
   v -> exp(-edot*dt/2)*v
   v -> v +f/m*dt
   v -> exp(-edot*dt/2)*v
   -----------------------------------------------------------------------*/
void FixNVEUefex::nve_v()
{
  double dtfm;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ex = erate[0]*dtf/2;
  double ey = erate[1]*dtf/2;
  double ez = -ex-ey;
  double e0 = exp(-ex);
  double e1 = exp(-ey);
  double e2 = exp(-ez);
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
        v[i][0] += dtfm*f[i][0];
        v[i][1] += dtfm*f[i][1];
        v[i][2] += dtfm*f[i][2];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
        v[i][0] += dtfm*f[i][0];
        v[i][1] += dtfm*f[i][1];
        v[i][2] += dtfm*f[i][2];
        v[i][0] *= e0;
        v[i][1] *= e1;
        v[i][2] *= e2;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Don't actually move atoms in remap(), just change the box
   -----------------------------------------------------------------------*/
void FixNVEUefex::remap()
{
  double vol = domain->xprd * domain->yprd * domain->zprd;
  double box[3][3];
  uefbox->get_box(box,vol);
  domain->boxhi[0] = domain->boxlo[0]+box[0][0];
  domain->boxhi[1] = domain->boxlo[1]+box[1][1];
  domain->boxhi[2] = domain->boxlo[2]+box[2][2];
  domain->xy = box[0][1];
  domain->xz = box[0][2];
  domain->yz = box[1][2];
  domain->set_global_box();
  domain->set_local_box();
  uefbox->get_rot(rot);
}

/* ----------------------------------------------------------------------
   SLLOD position update in time-reversible (i think) increments 
   x -> exp(edot*dt/2)*x
   x -> x + v*dt
   x -> exp(edot*dt/2)*x
   -----------------------------------------------------------------------*/
void FixNVEUefex::nve_x()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ex = erate[0]*dtv;
  strain[0] += ex;
  double e0 = exp(ex*0.5);
  double ey = erate[1]*dtv;
  strain[1] += ey;
  double e1 = exp(ey*0.5);
  double ez = -ex -ey;
  double e2 = exp(ez*0.5);
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x update by full step only for atoms in group
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] *= e0;
      x[i][1] *= e1;
      x[i][2] *= e2;
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
      x[i][0] *= e0;
      x[i][1] *= e1;
      x[i][2] *= e2;
    }
  }

  uefbox->step_deform(ex,ey);

  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  domain->boxhi[0] = domain->boxlo[0]+box[0][0];
  domain->boxhi[1] = domain->boxlo[1]+box[1][1];
  domain->boxhi[2] = domain->boxlo[2]+box[2][2];
  domain->xy = box[0][1];
  domain->xz = box[0][2];
  domain->yz = box[1][2];
  domain->set_global_box();
  domain->set_local_box();
  uefbox->get_rot(rot);

  set_h_rate(); // Murashima 2019/01/07
}

/* ----------------------------------------------------------------------
 * Do the lattice reduction if necessary.
 -----------------------------------------------------------------------*/
void FixNVEUefex::pre_exchange()
{
  // only need to reset things if the lattice needs to be reduced
  if (uefbox->reduce())
    {
      // go to lab frame
      inv_rotate_x(rot);
      inv_rotate_v(rot);
      inv_rotate_f(rot);
      // get & set the new box and rotation matrix
      double vol = domain->xprd * domain->yprd * domain->zprd;
      double box[3][3];
      uefbox->get_box(box,vol);
      domain->boxhi[0] = domain->boxlo[0]+box[0][0];
      domain->boxhi[1] = domain->boxlo[1]+box[1][1];
      domain->boxhi[2] = domain->boxlo[2]+box[2][2];
      domain->xy = box[0][1];
      domain->xz = box[0][2];
      domain->yz = box[1][2];
      domain->set_global_box();
      domain->set_local_box();
      uefbox->get_rot(rot);

      // rotate to the new upper triangular frame
      rotate_v(rot);
      rotate_x(rot);
      rotate_f(rot);

      // this is a generalization of what is done in domain->image_flip(...)
      int ri[3][3];
      uefbox->get_inverse_cob(ri);
      imageint *image = atom->image;
      int nlocal = atom->nlocal;
      for (int i=0; i<nlocal; i++) {
	int iold[3],inew[3];
	iold[0] = (image[i] & IMGMASK) - IMGMAX;
	iold[1] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
	iold[2] = (image[i] >> IMG2BITS) - IMGMAX;
	inew[0] = ri[0][0]*iold[0] + ri[0][1]*iold[1] + ri[0][2]*iold[2];
	inew[1] = ri[1][0]*iold[0] + ri[1][1]*iold[1] + ri[1][2]*iold[2];
	inew[2] = ri[2][0]*iold[0] + ri[2][1]*iold[1] + ri[2][2]*iold[2];
	image[i] = ((imageint) (inew[0] + IMGMAX) & IMGMASK) |
	  (((imageint) (inew[1] + IMGMAX) & IMGMASK) << IMGBITS) |
	  (((imageint) (inew[2] + IMGMAX) & IMGMASK) << IMG2BITS);
      }
    

      
      // put all atoms in the new box
      double **x = atom->x;
      for (int i=0; i<nlocal; i++) domain->remap(x[i],image[i]);

      // move atoms to the right processors
      domain->x2lamda(atom->nlocal);
      irregular->migrate_atoms();
      domain->lamda2x(atom->nlocal);
    }
}

/* ---------------------------------------------------------------------- 
 * The following are routines to rotate between the lab and upper triangular
 * (UT) frames. For most of the time the simulation is in the UT frame. 
 * To get to the lab frame, apply the inv_rotate_[..](rot) and to 
 * get back to the UT frame apply rotate_[..](rot). 
 *
 * Note: the rotate_x() functions also apply a shift to/from the fixedpoint
 * to make the integration a little simpler.
 * ---------------------------------------------------------------------- */
void FixNVEUefex::rotate_x(double r[3][3])
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double xn[3];
  for (int i=0;i<nlocal;i++)
    {
      if (mask[i] & groupbit)
	{
	  xn[0]=r[0][0]*x[i][0]+r[0][1]*x[i][1]+r[0][2]*x[i][2];
	  xn[1]=r[1][0]*x[i][0]+r[1][1]*x[i][1]+r[1][2]*x[i][2];
	  xn[2]=r[2][0]*x[i][0]+r[2][1]*x[i][1]+r[2][2]*x[i][2];
	  x[i][0]=xn[0]+domain->boxlo[0]; 
	  x[i][1]=xn[1]+domain->boxlo[1]; 
	  x[i][2]=xn[2]+domain->boxlo[2];
	}
    }
}

void FixNVEUefex::inv_rotate_x(double r[3][3])
{
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double xn[3];
  for (int i=0;i<nlocal;i++)
    {
      if (mask[i] & groupbit)
	{
	  x[i][0] -= domain->boxlo[0];
	  x[i][1] -= domain->boxlo[1];
	  x[i][2] -= domain->boxlo[2];
	  xn[0]=r[0][0]*x[i][0]+r[1][0]*x[i][1]+r[2][0]*x[i][2];
	  xn[1]=r[0][1]*x[i][0]+r[1][1]*x[i][1]+r[2][1]*x[i][2];
	  xn[2]=r[0][2]*x[i][0]+r[1][2]*x[i][1]+r[2][2]*x[i][2];
	  x[i][0]=xn[0];
	  x[i][1]=xn[1];
	  x[i][2]=xn[2];
	}
    }
}

void FixNVEUefex::inv_rotate_pos(double r[3][3],double p[3])
{
  double xn[3];
  p[0] -= domain->boxlo[0];
  p[1] -= domain->boxlo[1];
  p[2] -= domain->boxlo[2];
  xn[0]=r[0][0]*p[0]+r[1][0]*p[1]+r[2][0]*p[2];
  xn[1]=r[0][1]*p[0]+r[1][1]*p[1]+r[2][1]*p[2];
  xn[2]=r[0][2]*p[0]+r[1][2]*p[1]+r[2][2]*p[2];
  p[0]=xn[0];
  p[1]=xn[1];
  p[2]=xn[2];
}

void FixNVEUefex::rotate_v(double r[3][3])
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double vn[3];
  for (int i=0;i<nlocal;i++)
    {
      if (mask[i] & groupbit)
	{
	  vn[0]=r[0][0]*v[i][0]+r[0][1]*v[i][1]+r[0][2]*v[i][2];
	  vn[1]=r[1][0]*v[i][0]+r[1][1]*v[i][1]+r[1][2]*v[i][2];
	  vn[2]=r[2][0]*v[i][0]+r[2][1]*v[i][1]+r[2][2]*v[i][2];
	  v[i][0]=vn[0]; v[i][1]=vn[1]; v[i][2]=vn[2];
	}
    }
}

void FixNVEUefex::inv_rotate_v(double r[3][3])
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double vn[3];
  for (int i=0;i<nlocal;i++)
    {
      if (mask[i] & groupbit)
	{
	  vn[0]=r[0][0]*v[i][0]+r[1][0]*v[i][1]+r[2][0]*v[i][2];
	  vn[1]=r[0][1]*v[i][0]+r[1][1]*v[i][1]+r[2][1]*v[i][2];
	  vn[2]=r[0][2]*v[i][0]+r[1][2]*v[i][1]+r[2][2]*v[i][2];
	  v[i][0]=vn[0]; v[i][1]=vn[1]; v[i][2]=vn[2];
	}
    }
}

void FixNVEUefex::rotate_f(double r[3][3])
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double fn[3];
  for (int i=0;i<nlocal;i++)
    {
      if (mask[i] & groupbit)
	{
	  fn[0]=r[0][0]*f[i][0]+r[0][1]*f[i][1]+r[0][2]*f[i][2];
	  fn[1]=r[1][0]*f[i][0]+r[1][1]*f[i][1]+r[1][2]*f[i][2];
	  fn[2]=r[2][0]*f[i][0]+r[2][1]*f[i][1]+r[2][2]*f[i][2];
	  f[i][0]=fn[0]; f[i][1]=fn[1]; f[i][2]=fn[2];
	}
    }
}

void FixNVEUefex::inv_rotate_f(double r[3][3])
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  double fn[3];
  for (int i=0;i<nlocal;i++)
    {
      if (mask[i] & groupbit)
	{
	  fn[0]=r[0][0]*f[i][0]+r[1][0]*f[i][1]+r[2][0]*f[i][2];
	  fn[1]=r[0][1]*f[i][0]+r[1][1]*f[i][1]+r[2][1]*f[i][2];
	  fn[2]=r[0][2]*f[i][0]+r[1][2]*f[i][1]+r[2][2]*f[i][2];
	  f[i][0]=fn[0]; f[i][1]=fn[1]; f[i][2]=fn[2];
	}
    }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
   ------------------------------------------------------------------------- */

void FixNVEUefex::write_restart(FILE *fp)
{
  int nsize = size_restart_global();

  double *list;
  memory->create(list,nsize,"nve_uefex:list");

  pack_restart_data(list);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),nsize,fp);
  }

  memory->destroy(list);
}

/*
  ----------------------------------------------------------------------
  calculate the number of data to be packed
  ----------------------------------------------------------------------
*/
int FixNVEUefex::size_restart_global()
{
  return 2;
}

/*
  ----------------------------------------------------------------------
  pack restart data
  ----------------------------------------------------------------------
*/
int FixNVEUefex::pack_restart_data(double *list)
{
  int n = 0;
  list[n++] = strain[0];
  list[n++] = strain[1];
  return n;
}

/*
  ----------------------------------------------------------------------
  read and set the strains
  ----------------------------------------------------------------------
*/
void FixNVEUefex::restart(char *buf)
{
  int n = size_restart_global();
  double *list = (double *) buf;
  strain[0] = list[n-2];
  strain[1] = list[n-1];
  uefbox->set_strain(strain[0],strain[1]);
}

/*
  ----------------------------------------------------------------------
  If the step writes a restart, reduce the box beforehand. This makes sure
  the unique box shape can be found once the restart is read and that
  all of the atoms lie within the box.
  This may only be necessary for RESPA runs, but I'm leaving it in anyway.
  ----------------------------------------------------------------------
*/
void FixNVEUefex::end_of_step()
{
  if (update->ntimestep==output->next_restart)
    {
      pre_exchange();
      domain->x2lamda(atom->nlocal);
      set_h_rate(); // Murashima 2019/01/02
      domain->pbc();
      timer->stamp();
      comm->exchange();
      comm->borders();
      domain->lamda2x(atom->nlocal+atom->nghost);
      timer->stamp(Timer::COMM);
      neighbor->build(1);
      timer->stamp(Timer::NEIGH);
    }
}

/*
  ----------------------------------------------------------------------
  reduce the simulation box after a run is complete. otherwise it won't
  be possible to resume from a write_restart since the initialization of
  the simulation box requires reduced simulation box
  ----------------------------------------------------------------------
*/
void FixNVEUefex::post_run()
{
  pre_exchange();
  domain->x2lamda(atom->nlocal);
  set_h_rate();// Murashima 2019/01/02
  domain->pbc();
  timer->stamp();
  comm->exchange();
  comm->borders();
  domain->lamda2x(atom->nlocal+atom->nghost);
  timer->stamp(Timer::COMM);
  neighbor->build(1);
  timer->stamp(Timer::NEIGH);
}

/*
  ----------------------------------------------------------------------
  public read for rotation matrix
  ----------------------------------------------------------------------
*/
void FixNVEUefex::get_rot(double r[3][3])
{
  uefbox->get_rot(rot); // 2017/11/08 Murashima
  r[0][0] = rot[0][0];
  r[0][1] = rot[0][1];
  r[0][2] = rot[0][2];
  r[1][0] = rot[1][0];
  r[1][1] = rot[1][1];
  r[1][2] = rot[1][2];
  r[2][0] = rot[2][0];
  r[2][1] = rot[2][1];
  r[2][2] = rot[2][2];
}

/*
  ----------------------------------------------------------------------
  public read for ext flags
  ----------------------------------------------------------------------
*/
void FixNVEUefex::get_ext_flags(bool* e)
{
  e[0] = ext_flags[0];
  e[1] = ext_flags[1];
  e[2] = ext_flags[2];
}

/*
  ----------------------------------------------------------------------
  public read for simulation box
  ----------------------------------------------------------------------
*/
void FixNVEUefex::get_box(double b[3][3])
{
  double box[3][3];
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol);
  b[0][0] = box[0][0];
  b[0][1] = box[0][1];
  b[0][2] = box[0][2];
  b[1][0] = box[1][0];
  b[1][1] = box[1][1];
  b[1][2] = box[1][2];
  b[2][0] = box[2][0];
  b[2][1] = box[2][1];
  b[2][2] = box[2][2];
}

/*
  ----------------------------------------------------------------------
  comparing floats
  it's imperfect, but should work provided no infinities
  ----------------------------------------------------------------------
*/
bool FixNVEUefex::nearly_equal(double a, double b, double epsilon)
{
  double absa = fabs(a);
  double absb = fabs(b);
  double diff = fabs(a-b);
  if (a == b) return true;
  else if ( (absa+absb) < epsilon)
    return diff < epsilon*epsilon;
  else
    return diff/(absa+absb) < epsilon;
}


// Murashima 2019/01/02 - 2019/01/07
void FixNVEUefex::set_h_rate()
{
  double ex,ey,ez,er[6],rot[3][3],trot[3][3];
  double er9[3][3],h9[3][3],h_rate9[3][3],h_rate9_lab[3][3];
  double h9_lab[3][3];
  double e9_lab[3][3];
  double er9_shear[3][3];
  double box[3][3];
  uefbox->get_rot(rot); // Murashima 2019/01/02
  mat_trans(rot,trot); // Murashima 2019/01/05
  double vol = domain->xprd * domain->yprd * domain->zprd;
  uefbox->get_box(box,vol); // Murashima 2019/01/02
  ex=erate[0];
  ey=erate[1];
  ez=-ex-ey;
  er[0]=ex;
  er[1]=ey;
  er[2]=ez;
  er[3]=0.0;
  er[4]=0.0;
  er[5]=0.0;
  // er9 is symmetric tensor
  symmat_rot(er,rot,er9); // LAB -> UT
  //
  // for general flow, symmat_rot should be modified to mat_rot
  //

  // Murashima 2019/01/07
  er9_shear[0][0]=er9[0][0];
  er9_shear[1][1]=er9[1][1];
  er9_shear[2][2]=er9[2][2];
  er9_shear[1][2]=2.0*er9[1][2];// er9[1][2]+er9[2][1]
  er9_shear[0][2]=2.0*er9[0][2];// er9[0][2]+er9[2][0]
  er9_shear[0][1]=2.0*er9[0][1];// er9[0][1]+er9[1][0]
  er9_shear[2][1]=0.0;
  er9_shear[2][0]=0.0;
  er9_shear[1][0]=0.0;

  // h9 is upper trianglular matrix

  h9[0][0]=box[0][0];
  h9[1][1]=box[1][1];
  h9[2][2]=box[2][2];
  h9[1][2]=box[1][2];
  h9[0][2]=box[0][2];
  h9[0][1]=box[0][1];
  h9[2][1]=0.0;
  h9[2][0]=0.0;
  h9[1][0]=0.0;


  // Murashima 2019/01/07
  h_rate9[0][0]=er9_shear[0][0]*box[0][0];
  h_rate9[1][0]=0.0;
  h_rate9[2][0]=0.0;
  h_rate9[0][1]=er9_shear[0][1]*box[1][1];
  h_rate9[1][1]=er9_shear[1][1]*box[1][1];
  h_rate9[2][1]=0.0;
  h_rate9[0][2]=er9_shear[0][2]*box[2][2];
  h_rate9[1][2]=er9_shear[1][2]*box[2][2];
  h_rate9[2][2]=er9_shear[2][2]*box[2][2];

  domain->h_rate[0]=h_rate9[0][0];
  domain->h_rate[1]=h_rate9[1][1];
  domain->h_rate[2]=h_rate9[2][2];
  domain->h_rate[3]=h_rate9[1][2];
  domain->h_rate[4]=h_rate9[0][2];
  domain->h_rate[5]=h_rate9[0][1];


  
  domain->h_ratelo[0]=domain->h_ratelo[1]=domain->h_ratelo[2]=0.0;
  // boxlo is the stagnation point!
}

// Murashima 2019/01/02
void FixNVEUefex::mat_mul(const double a[3][3],const double b[3][3],double c[3][3])
{
  for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
	{
	  c[i][j]=0.0;
	  for(int k=0;k<3;k++)c[i][j]+=a[i][k]*b[k][j];
	}
    }
}

// Murashima 2019/01/02
void FixNVEUefex::mat_trans(const double a[3][3],double at[3][3])
{
  for(int i=0;i<3;i++)
    {
      at[0][i]=a[i][0];
      at[1][i]=a[i][1];
      at[2][i]=a[i][2];
    }
}

// Murashima 2019/01/02
// a:LAB -> ao:LAMMPS(UT)
// ao=rot.a.trot
void FixNVEUefex::symmat_rot(const double a[6],const double rot[3][3],double ao[3][3])
{
  double trot[3][3],a9[3][3],t[3][3];
  mat_trans(rot,trot);
  a9[0][0]=a[0];
  a9[1][1]=a[1];
  a9[2][2]=a[2];
  a9[1][2]=a9[2][1]=a[3];
  a9[0][2]=a9[2][0]=a[4];
  a9[0][1]=a9[1][0]=a[5];

  mat_mul(a9,trot,t);// t=a9.trot
  mat_mul(rot,t,ao);// ao=rot.t
}

// Murashima 2019/01/05
void FixNVEUefex::mat_rot(const double a[3][3],const double rot[3][3],double ao[3][3])
{
  double trot[3][3],t[3][3];
  mat_trans(rot,trot);
  mat_mul(a,trot,t);// t=a.trot
  mat_mul(rot,t,ao);// ao=rot.t
}

// Murashima 2019/01/05
void FixNVEUefex::mat_irot(const double a[3][3],const double rot[3][3],double ao[3][3])
{
  double trot[3][3],t[3][3];
  mat_trans(rot,trot);
  mat_mul(a,rot,t);// t=a.rot
  mat_mul(trot,t,ao);// ao=trot.t
}



// Murashima 2018/12/25
void FixNVEUefex::get_erate(double er[3])
{
  double ex,ey,ez;
  ex=erate[0];
  ey=erate[1];
  ez=-ex-ey;
  er[0]=ex;
  er[1]=ey;
  er[2]=ez;
}

// Murashima 2018/12/25
void FixNVEUefex::get_strain(double er[3])
{
  double ex,ey,ez;
  ex=strain[0];
  ey=strain[1];
  ez=-ex-ey;
  er[0]=ex;
  er[1]=ey;
  er[2]=ez;
}

// Murashima 2018/12/25
void FixNVEUefex::inv_rotate() // LAMMPS(UT) to LAB
{
  uefbox->get_rot(rot);
  inv_rotate_x(rot);
  inv_rotate_v(rot);
  inv_rotate_f(rot);
  ((ComputeTempUefex*) temperature)->no_rot();
}

// Murashima 2018/12/25
void FixNVEUefex::rotate() // LAB to LAMMPS(UT)
{
  uefbox->get_rot(rot);
  rotate_x(rot);
  rotate_v(rot);
  rotate_f(rot);
  ((ComputeTempUefex*) temperature)->yes_rot();
}


// Murashima 2018/12/28
/*
  ----------------------------------------------------------------------
  Transform the pressure tensor to the rotated coordinate system
  [P]rot = Q.[P].Q^t
  [P]    : LAMMPS(UT)
  [P]rot : LAB
  -------------------------------------------------------------------------
*/
/*
  void FixNVEUefex::virial_rot(double *x, const double r[3][3])
  {

  double t[3][3];
  // Original UEF
  // [00 10 20 ] [ 0 3 4 ] [00 01 02 ]
  // [01 11 21 ] [ 3 1 5 ] [10 11 12 ]
  // [02 12 22 ] [ 4 5 2 ] [20 21 22 ]
  //
  // Here
  // [00 10 20 ] [ 0 5 4 ] [00 01 02 ]
  // [01 11 21 ] [ 5 1 3 ] [10 11 12 ]
  // [02 12 22 ] [ 4 3 2 ] [20 21 22 ]
  // This is consistent with domain.cpp
  //
  for (int k = 0; k<3; ++k) 
  {
  t[0][k] = x[0]*r[0][k] + x[5]*r[1][k] + x[4]*r[2][k];
  t[1][k] = x[5]*r[0][k] + x[1]*r[1][k] + x[3]*r[2][k];
  t[2][k] = x[4]*r[0][k] + x[3]*r[1][k] + x[2]*r[2][k];
  }
  x[0] = r[0][0]*t[0][0] + r[1][0]*t[1][0] + r[2][0]*t[2][0];
  x[5] = r[0][0]*t[0][1] + r[1][0]*t[1][1] + r[2][0]*t[2][1];
  x[4] = r[0][0]*t[0][2] + r[1][0]*t[1][2] + r[2][0]*t[2][2];
  x[1] = r[0][1]*t[0][1] + r[1][1]*t[1][1] + r[2][1]*t[2][1];
  x[3] = r[0][1]*t[0][2] + r[1][1]*t[1][2] + r[2][1]*t[2][2];
  x[2] = r[0][2]*t[0][2] + r[1][2]*t[1][2] + r[2][2]*t[2][2];
  }
*/

// Murashima 2018/12/28
/*
  ----------------------------------------------------------------------
  Transform the symmetric tensor to the rotated coordinate system
  [P] = Q^t.[P]rot.Q
  [P]    : LAMMPS(UT)
  [P]rot : LAB
  -------------------------------------------------------------------------
*/

void FixNVEUefex::inv_virial_rot(double *x, const double r[3][3])
{

  double t[3][3];
  
  // Original UEF 
  // [00 01 02 ] [ 0 3 4 ] [00 10 20 ]
  // [10 11 12 ] [ 3 1 5 ] [01 11 21 ]
  // [20 21 22 ] [ 4 5 2 ] [02 12 22 ]
  //
  // Here
  // [00 01 02 ] [ 0 5 4 ] [00 10 20 ]
  // [10 11 12 ] [ 5 1 3 ] [01 11 21 ]
  // [20 21 22 ] [ 4 3 2 ] [02 12 22 ]
  // This is consistent with domain.cpp
  //
  for (int k = 0; k<3; ++k) 
    {
      t[0][k] = x[0]*r[k][0] + x[5]*r[k][1] + x[4]*r[k][2];
      t[1][k] = x[5]*r[k][0] + x[1]*r[k][1] + x[3]*r[k][2];
      t[2][k] = x[4]*r[k][0] + x[3]*r[k][1] + x[2]*r[k][2];
    }
  x[0] = r[0][0]*t[0][0] + r[0][1]*t[1][0] + r[0][2]*t[2][0];
  x[5] = r[0][0]*t[0][1] + r[0][1]*t[1][1] + r[0][2]*t[2][1];
  x[4] = r[0][0]*t[0][2] + r[0][1]*t[1][2] + r[0][2]*t[2][2];
  x[1] = r[1][0]*t[0][1] + r[1][1]*t[1][1] + r[1][2]*t[2][1];
  x[3] = r[1][0]*t[0][2] + r[1][1]*t[1][2] + r[1][2]*t[2][2];
  x[2] = r[2][0]*t[0][2] + r[2][1]*t[1][2] + r[2][2]*t[2][2];
}



// Murashima 2018/12/31
void FixNVEUefex::mat_mul6(const double a[6],const double b[6],double c[6])
{
  // Multiply upper triangle matrix
  // [0 5 4]   [xx xy xz]
  // [* 1 3] = [*  yy yz]
  // [* * 2]   [*  *  zz]
  //
  // A B
  // [xx xy xz][xx xy xz]
  // [ 0 yy yz][ 0 yy yz]
  // [ 0  0 zz][ 0  0 zz]
  // =
  // [xx(0).xx(0)  xx(0).xy(5)+xy(5).yy(1)  xx(0).xz(4)+xy(5).yz(3)+xz(4).zz(2)]
  // [0            yy(1).yy(1)              yy(1).yz(3)+yz(3).zz(2)            ]
  // [0            0                        zz(2).zz(2)                        ]
  // =
  // [0.0 0.5+5.1 0.4+5.3+4.2]
  // [0   1.1     1.3+3.2    ]
  // [0   0       2.2        ]
  // C
  c[0]=a[0]*b[0];
  c[1]=a[1]*b[1];
  c[2]=a[2]*b[2];
  c[3]=a[1]*b[3]+a[3]*b[2];
  c[4]=a[0]*b[4]+a[5]*b[3]+a[4]*b[2];
  c[5]=a[0]*b[5]+a[5]*b[1];
}

