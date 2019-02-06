/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: T. Murashima (Tohoku Univ, JPN)
   ------------------------------------------------------------------------- */


#ifdef FIX_CLASS
FixStyle(nve/uefex,FixNVEUefex)
#else

#ifndef LMP_FIX_NVE_UEFEX_H
#define LMP_FIX_NVE_UEFEX_H

//#include "fix_nh.h"
#include "fix.h"

namespace LAMMPS_NS {
  namespace UEF_utils {
    class UEFBox;
  }
  
  class FixNVEUefex : public Fix {
  public:
    FixNVEUefex(class LAMMPS *, int, char **);
    virtual ~FixNVEUefex();
    virtual int setmask(); 
    virtual void init();
    virtual void setup(int); 
    virtual void pre_exchange();
    void write_restart(FILE *); 
    virtual int pack_restart_data(double*);
    virtual void restart(char *);
    virtual void end_of_step();
    virtual void initial_integrate(int);
    virtual void final_integrate();
    virtual void initial_integrate_respa(int, int, int);
    virtual void final_integrate_respa(int, int);
    virtual void post_run();

    // fix_nh
    int modify_param(int, char **);
    //void reset_target(double);
    //void reset_dt();
    //virtual void *extract(const char*, int &);
    //double memory_usage();
  
    void get_rot(double[3][3]);
    void get_ext_flags(bool*);
    void get_box(double[3][3]);
    void rotate_x(double [3][3]);
    void inv_rotate_x(double[3][3]);
    void rotate_v(double[3][3]);
    void inv_rotate_v(double[3][3]);
    void rotate_f(double[3][3]);
    void inv_rotate_f(double[3][3]);
    void get_erate(double[3]); // Murashima 2018/12/25
    void get_strain(double[3]); // Murashima 2018/12/25
    void inv_rotate(); // Murashima 2018/12/25
    void rotate(); // Murashima 2018/12/25
    void inv_rotate_pos(double[3][3],double [3]);// Murashima 2019/01/25

    //void virial_rot(double *x, const double r[3][3]);// Murashima 2018/12/28 LAMMPS -> LAB
    void inv_virial_rot(double *x, const double r[3][3]);// Murashima 2019/01/02 LAB -> LAMMPS
    //void get_h_rate(double[6]); // Murashima 2018/12/28
    //void save_h_pre(); // Murashima 2018/12/31
    void set_h_rate(); // Murashima 2019/01/02
    void mat_mul6(const double a[6],const double b[6],double c[6]); // Murashima 2019/01/02
    void mat_mul(const double a[3][3],const double b[3][3],double c[3][3]); // Murashima 2019/01/02

    void symmat_rot(const double a[6],const double rot[3][3],double ao[3][3]); // Murashima 2019/01/02
    void mat_rot(const double a[3][3],const double rot[3][3],double ao[3][3]); // Murashima 2019/01/05
    void mat_irot(const double a[3][3],const double rot[3][3],double ao[3][3]); // Murashima 2019/01/05
    //void symmat_irot(const double a[6],const double rot[3][3],double ao[3][3]);
    void mat_trans(const double a[3][3],double at[3][3]);
    //void set_h_rate_from_box(); // Murashima 2019/01/02
    int remapflag; // same as fix_deform.h


  

  protected:
    virtual void remap();
    virtual void nve_x();
    virtual void nve_v();
    virtual int size_restart_global(); // fix_nh
    //char **arg_kludge(int&, char**, int&);
    //int  narg_kludge(int&, char**); // add by Murashima 2017/11/08

    double strain[2],erate[2]; // strain/strain rate : [e_x, e_y]
    // always assume traceless e_z = -e_x-e_y

    int rem;                   //this is for the narg kluge

    UEF_utils::UEFBox *uefbox;            // interface for the special simulation box

    double rot[3][3];          // rotation matrix
    bool ext_flags[3];         // flags for external "free surfaces"
    bool nearly_equal(double,double,double);
    //bool rotate_output;      // experimental feature. Too many issues for now

    double h_rate[6]; // 2018/12/30
    double h_pre[6]; // Murashima 2018/12/31
    double h_now[6]; // Murashima 2018/12/31
    double h_now_inv[6]; // Murashima 2018/12/31
    double dt_inv;

    int dimension,which;
    double dtv,dtf,dthalf,dto;

    int nrigid;                      // number of rigid fixes
    int *rfix;                       // indices of rigid fixes

    int allremap;
    int dilate_group_bit;            // mask for dilation group
    char *id_dilate;                 // group name to dilate
  
  
    class Irregular *irregular;      // for migrating atoms after box flips
    int nlevels_respa;
    double *step_respa;
  
    char *id_temp,*id_press;
    class Compute *temperature,*pressure;
    int tcomputeflag,pcomputeflag;   // 1 = compute was created by fix
    // 0 = created externally

    int pre_exchange_flag;           // set if pre_exchange needed for box flips

    double fixedpoint[3];            // location of dilation fixed-point
  
  
  };

}

#endif
#endif

/* ERROR/WARNING messages:

   This is a base class for FixNVEUefex so it will inherit most of its error/warning messages along with the following:

   E: Illegal fix nve/uefex command

   Self-explanatory

   E: Keyword erate must be set for fix nve/uefex command

   Self-explanatory.

   E: Simulation box must be triclinic for fix nve/uefex

   Self-explanatory.

   E: Can't use another fix which changes box shape with fix nve/uefex

   The fix nve/uefex command must have full control over the box shape. You cannot use a simultaneous fix deform command, for example.

   E: Pressure ID for fix nve/uefex doesn't exist

   The compute pressure introduced via fix_modify does not exist

   E: Using fix nve/uefex without a compute pressure/uefex

   Self-explanatory.

   E: Using fix nve/uefex without a compute temp/uefex

   Self-explanatory.

   E: Initial box is not close enough to the expected uef box

   The initial box does not correspond to the shape required by the value of the strain keyword. If the default strain value of zero was used, the initial box is not cubic.

*/
