/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/*
  This code is based on dump_xyz.cpp
  Modified by T. Murashima (Tohoku Univ, JPN)
 */


#include <string.h>
#include "dump_xyz_uefex.h" // 2017/12/20 Murashima
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include "modify.h"
#include "domain.h"       // 2017/12/20 Murashima
#include "fix_nve_uefex.h" // 2017/12/20 Murashima
#include "uef_utils.h"

using namespace LAMMPS_NS;

#define ONELINE 128
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpXYZUefex::DumpXYZUefex(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  typenames(NULL)
{
  if (narg != 5) error->all(FLERR,"Illegal dump xyz/uefex command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump xyz/uefex filename");

  size_one = 5;

  buffer_allow = 1;
  buffer_flag = 1;
  sort_flag = 1;
  sortcol = 0;

  if (format_default) delete [] format_default;

  char *str = (char *) "%s %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);

  ntypes = atom->ntypes;
  typenames = NULL;


}

/* ---------------------------------------------------------------------- */

DumpXYZUefex::~DumpXYZUefex()
{
  delete[] format_default;
  format_default = NULL;

  if (typenames) {
    for (int i = 1; i <= ntypes; i++)
      delete [] typenames[i];
    delete [] typenames;
    typenames = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZUefex::init_style()
{
  // format = copy of default or user-specified line format

  delete [] format;
  char *str;
  if (format_line_user) str = format_line_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // initialize typenames array to be backward compatible by default
  // a 32-bit int can be maximally 10 digits plus sign

  if (typenames == NULL) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[12];
      sprintf(typenames[itype],"%d",itype);
    }
  }

  // setup function ptr

  if (buffer_flag == 1) write_choice = &DumpXYZUefex::write_string;
  else write_choice = &DumpXYZUefex::write_lines;

  // open single file, one time only

  if (multifile == 0) openfile();

  
  // 2017/12/20 Murashima

  int i=0;
  for (i=0; i<modify->nfix; i++)
  {
    if (strcmp(modify->fix[i]->style,"nve/uefex")==0)
      break;
  }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use dump xyz/uefex without defining a fix nve/uefex");
  ifix_uefex=i;

}

/* ---------------------------------------------------------------------- */

int DumpXYZUefex::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR, "Dump modify element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++)
        delete [] typenames[i];

      delete [] typenames;
      typenames = NULL;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      int n = strlen(arg[itype]) + 1;
      typenames[itype] = new char[n];
      strcpy(typenames[itype],arg[itype]);
    }

    return ntypes+1;
  }

  return 0;
}



/* ---------------------------------------------------------------------- */

void DumpXYZUefex::write_header(bigint n)
{

  double box[3][3],rot[3][3];
  ((FixNVEUefex*) modify->fix[ifix_uefex])->get_box(box);
  ((FixNVEUefex*) modify->fix[ifix_uefex])->get_rot(rot);

  for(int i=0;i<3;i++)
    {
      for(int j=i+1;j<3;j++)
	{
	  double t=rot[i][j];
	  rot[i][j]=rot[j][i];
	  rot[j][i]=t;
	}
    }
  UEF_utils::mul_m2(rot,box);
  

  
  if (me == 0) {
    fprintf(fp,BIGINT_FORMAT "\n",n);
    fprintf(fp,"Atoms. Timestep: " BIGINT_FORMAT " ",update->ntimestep);
    fprintf(fp,"Lattice=\" %g %g %g %g %g %g %g %g %g \" \n",
	    box[0][0],box[1][0],box[2][0],
	    box[0][1],box[1][1],box[2][1],
	    box[0][2],box[1][2],box[2][2]
	    );
  }
}

/* ---------------------------------------------------------------------- */

void DumpXYZUefex::pack(tagint *ids)
{

  double rot[3][3];
  ((FixNVEUefex*) modify->fix[ifix_uefex])->get_rot(rot);
  ((FixNVEUefex*) modify->fix[ifix_uefex])->inv_rotate_x(rot);



  int m,n;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      if (ids) ids[n++] = tag[i];
    }

  ((FixNVEUefex*) modify->fix[ifix_uefex])->rotate_x(rot);

}


/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpXYZUefex::convert_string(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset],format,
                      typenames[static_cast<int> (mybuf[m+1])],
                      mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */


void DumpXYZUefex::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */
void DumpXYZUefex::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}


/* ---------------------------------------------------------------------- */

void DumpXYZUefex::write_lines(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            typenames[static_cast<int> (mybuf[m+1])],
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }
}

