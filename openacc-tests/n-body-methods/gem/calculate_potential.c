/****************************************************************************
 * GEM -- electrostatics calculations and visualization                     *
 * Copyright (C) 2006  John C. Gordon                                       *
 *                                                                          *
 * This program is free software; you can redistribute it and/or modify     *
 * it under the terms of the GNU General Public License as published by     *
 * the Free Software Foundation; either version 2 of the License, or        *
 * (at your option) any later version.                                      *
 *                                                                          *
 * This program is distributed in the hope that it will be useful,          *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 * GNU General Public License for more details.                             *
 *                                                                          *
 * You should have received a copy of the GNU General Public License along  *
 * with this program; if not, write to the Free Software Foundation, Inc.,  *
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.              *
 ****************************************************************************/
#include <sys/time.h>
#include "calculations.h"
#include "defines.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fcntl.h>
#include "ctsar.h"
typedef long long int64;
static struct timeval start_time;
void init_time()
{
  gettimeofday(&start_time, NULL);
}

int64 get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (int64) (t.tv_sec - start_time.tv_sec) * 1000000
         + (t.tv_usec - start_time.tv_usec);
}

/*****************************************************************************
 *FUNCTION -- dist2                                                          *
 *                                                                           *
 *INPUTS:  x, y, z    (position of first point in space)                     *
 *         x2, y2, z2 (position of second point in space)                    *
 *                                                                           *
 *RETURNS: length of the vector from (x,y,z) to (x2, y2, z2) squared         *
 *             (good for length squared operations to avoid square root and  *
 *               then multiply) --jcg                                        *
 *                                                                           *
 *****************************************************************************/
/* inline float dist2(float x, float y, float z, float x2, float y2, float z2)
{
    /\* local variables *\/
    double xd = x - x2,
           yd = y - y2,
           zd = z - z2;

    return (xd*xd + yd*yd + zd*zd);

} /\* end function dist2 *\/ */
/* #define dist2(x, y, z, x2, y2, z2) \
        ((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2))

#define dist(x,y,z,x2,y2,z2) \
        sqrt(dist2(x, y, z, x2, y2, z2)) */

/*****************************************************************************
 *FUNCTION -- dist                                                           *
 *                                                                           *
 *INPUTS:  x, y, z    (position of first point in space)                     *
 *         x2, y2, z2 (position of second point in space)                    *
 *                                                                           *
 *RETURNS: length of the vector from (x,y,z) to (x2, y2, z2)                 *
 *                                                                           *
 *****************************************************************************/

/* inline float dist (float x, float y, float z, float x2, float y2, float z2)
{
    return sqrt(dist2(x, y, z, x2, y2, z2));
} /\* end function dist *\/ */
/*
 * There is a figure in this directory called "surfacepot.fig"
 * demonstrating the computation of these potentials and the
 * relative distances, etc. --jcg
 */

/* If we intend to continue work on this we have to __HAVE TO HAVE TO__
 * unify the positional calculations used in both methods, by writing
 * a subfunction that takes position (x, y, z) and vert and res and
 * returns the potential at that point (will also be useful for
 * arbitrary sample point distributions) --jcg
 */

/***************************************************************************
 * FUNCTION: calculate_potential_contribution -- calculates the potential  *
 *                                             contribution of one atom    *
 *                                             to the total potential of a *
 *                                             given point in space. Note  *
 *                                             that while this function is *
 *                                             slower (slightly) it is     *
 *                                             really easy to use.         *
 *                                                                         *
 * INPUTS: r            -- distance from point to charge                   *
 *         rprime       -- distance from point to surface along r          *
 *         d            -- distance from point to center                   *
 *         dprime       -- distance from point to surface along d          *
 *         A            -- electrostatic radius of the molecule            *
 *         r0           -- rnaut in the formulations                       *
 *         sal          -- salt concentration                              *
 *         charge       -- charge under consideration                      *
 *         diel_ext     -- external dielectric environment                 *
 *         diel_int     -- internal dielectric environment                 *
 *         defs         -- definitions to speed things up                  *
 *                                                                         *
 * RETURNS: the potential contribution                                     *
 *                                                                         *
 ***************************************************************************/
inline float calculate_potential_contribution (float d, float dprime, float r,
    float rprime, float r0, float A,
    float sal, float charge,
    float diel_int, float diel_ext,
    analytical_definitions_struct *defs)
{
  /* local variables */
  /*******************/
  float sum1 = 0,
        sum2 = 0,
        sum3 = 0,
        sum4 = 0,
        salt = 0,
        coulomb = 0,
        to_return = 0,
        one_over_one_plus_kappa_rprime;

  /* <><> add salt to the inside <><> */

  if (defs->phiType != TOTAL_POTENTIAL)
  {
    coulomb = charge / (d * diel_int);
  }

  if (defs->phiType & REACTION_POTENTIAL)
  {
    if (defs->region == 3)  /* Dist to surf > ion_excl_rad */
    {

      sum1 = (charge
              * (exp (-defs->kappa * (d - dprime))) / d)
             * (defs->one_plus_alpha / (1. + defs->kappa * dprime));

      sum2 = charge
             * (defs->alpha_by_one_minus_beta * exp (-defs->kappa * (r - rprime)) / (r * (1. + defs->kappa * rprime)));

      to_return = defs->inverse_one_plus_ab_by_diel_ext * (sum1 - sum2);
    }
    else if (defs->region == 2) /* 0 < dist_to_surf <= ion_excl_rad */
    {

      one_over_one_plus_kappa_rprime = 1. / (1. + (defs->kappa * rprime));

      /* electrostatic terms for region 2 */
      sum1 = defs->one_plus_alpha * charge / d;
      sum2 = defs->alpha_by_one_minus_beta * charge / r;

      /* salt terms for region 2 */
      sum3 = charge * (defs->alpha_by_one_minus_beta / rprime) * (1. - one_over_one_plus_kappa_rprime);
      sum4 = defs->one_plus_alpha * charge
             * ( 1.0 / (dprime + defs->kappa * dprime * dprime) - (1. / dprime));

      to_return = defs->inverse_one_plus_ab_by_diel_ext * (sum1 - sum2 + sum3 + sum4);
    }
    else /* region = 1 */
    {

      one_over_one_plus_kappa_rprime = 1. / (1. + (defs->kappa * rprime));

      sum1 =  1. / (d * diel_int);

      sum2 =  (1. / diel_int - 1. / diel_ext) / (defs->one_plus_alpha_beta * A);

      sum3 = defs->Asq / sqrt(defs->Asq_minus_rnautsq * defs->Asq_minus_rsq + defs->Asq_by_dsq);

      salt = defs->inverse_one_plus_ab_by_diel_ext
             * (defs->one_plus_alpha / dprime * (1. / (1. + defs->kappa * dprime) - 1.)
                + defs->alpha_by_one_minus_beta / rprime * (1. - one_over_one_plus_kappa_rprime));

      to_return = (sum1 - sum2 * sum3 + salt) * charge;
      /* to_return = d;//(sum1 - sum2*sum3 + salt) * atom_charge[j+atom_index[k]]; */
    }

    if (defs->phiType == REACTION_POTENTIAL)
    {
      to_return -= coulomb;
    }
  }
  else
  {
    to_return = coulomb;
  }


  return to_return;
}

typedef struct residue_noat
{
  int  natoms;   /* this residue number         */

  float x, y, z;
  float total_charge; /* total charge of the residue */
} residue_noat;
/* Calculates potential for surface out of a molecule -- this is the free
 * space solution. --JCG */
void calc_potential_single_step(residue *residues, int nres, vertx *vert, int nvert, float A, float proj_len, float diel_int, float diel_ext, float sal, float ion_exc_rad, int phiType, int *i, int step_size)
{
  /* local variables */
  /*******************/
  int eye, k, bound;
  float r,
        r0,
        rprime,
        d,
        dprime;

  analytical_definitions_struct defs;

  eye = *i;

  r = A + proj_len;
  rprime = A + ion_exc_rad;

  defs.Asq   = A * A;
  defs.Asq_minus_rsq = defs.Asq - (r * r);
  float Asq_minus_rsq = defs.Asq - (r * r);

  defs.kappa = 0.316 * sqrt(sal);
  defs.beta = diel_int / diel_ext;
  defs.alpha_by_one_minus_beta = ALPHA_OF * (1. - defs.beta);
  defs.alpha_beta = ALPHA_OF * defs.beta;
  defs.one_plus_alpha_beta = 1.0 + defs.alpha_beta;
  defs.one_plus_alpha = 1.0 + ALPHA_OF;
  defs.inverse_one_plus_ab_by_diel_ext = 1.0 / (defs.one_plus_alpha_beta * diel_ext);
  defs.phiType = phiType;
  float kappa = 0.316 * sqrt(sal);
  float beta = diel_int / diel_ext;
  float alpha_by_one_minus_beta = ALPHA_OF * (1. - defs.beta);
  float alpha_beta = ALPHA_OF * defs.beta;
  float one_plus_alpha_beta = 1.0 + defs.alpha_beta;
  float one_plus_alpha = 1.0 + ALPHA_OF;
  float inverse_one_plus_ab_by_diel_ext = 1.0 / (defs.one_plus_alpha_beta * diel_ext);
  float one_over_one_plus_kappa_rprime = 1. / (1. + (kappa * rprime));

  int total_atoms = 0, *atit, it;
  residue_noat *resa = (residue_noat*)malloc(sizeof(residue_noat) * nres);
  for (it = 0; it < nres; it++)
  {
    total_atoms += residues[it].natoms;
    resa[it].natoms = residues[it].natoms;
    resa[it].x = residues[it].x;
    resa[it].y = residues[it].y;
    resa[it].z = residues[it].z;
    resa[it].total_charge = residues[it].total_charge;
  }
  atom *ata = (atom*)ctsar_calloc_test(sizeof(atom), total_atoms);
  atit = (int*)ctsar_calloc_test(sizeof(int), nres);
  int outit = 0;
  for (it = 0; it < nres; it++)
  {
    atit[it] = outit;
    int j;
    for (j = 0; j < residues[it].natoms; j++)
    {
      ata[outit] = residues[it].atoms[j];
      outit++;
    }
  }

  if (r > rprime)
  {
    defs.region = 3;
  }
  else if (r > A)
  {
    defs.region = 2;
  }
  else
  {
    defs.region = 1;
  }
  int region = defs.region;

  bound = eye + step_size;

  /* for(it=0;it<acc_get_num_devices(acc_device_nvidia); it++){
      acc_set_device_num(it, acc_device_nvidia);
      acc_init(acc_device_nvidia);
  } */
  /* acc_set_device_num(0, acc_device_nvidia);
  acc_init(acc_device_nvidia); */
  /* acc_init(acc_device_nvidia); */
  /* fprintf(stderr,"num acc devices %d\n", acc_get_num_devices(acc_device_nvidia)); */
  /* #pragma omp parallel num_threads(8)
      {
          int tid = omp_get_thread_num();

              if(tid <= acc_get_num_devices(acc_device_nvidia)){
                  acc_set_device_num(tid, acc_device_nvidia);
                  printf("selected device number %d of %d, got %d\n", tid, acc_get_num_devices(acc_device_nvidia),acc_get_device_num(acc_device_nvidia));
              }
      } */

  init_time();
  if (nvert < bound)
  {
    bound = nvert;
  }


  float *pots = (float*)ctsar_calloc_test(sizeof(float), bound);
  /* float *pots2 = (float*)malloc(sizeof(float)*bound); */

  /* need to estimate r0, dist from atom to surface */
  /* if we intend to support negative projection lengths <><> */
  r0 = A / 2.; /* might be tricky ... hrmmm */

  defs.Asq_minus_rnautsq = defs.Asq - r0 * r0;
  float Asq = A * A;
  float Asq_minus_rnautsq = defs.Asq - r0 * r0;

  char * env;

  int iters;
  if (env = getenv("TEST_ITERATIONS"))
  {
    iters = atoi(env);
  }
  else
  {
    iters = 1;
  }
  ctsar * s = NULL;
  int safe_it;
  for (safe_it = 0; safe_it < iters; safe_it++)
  {
    if ((env = getenv("OMP_CTSAR_FALLBACK")) != NULL && atoi(env) != 0)
    {
      printf("fallback active\n");
      #pragma omp parallel for schedule(runtime) private(eye, d, dprime, k) firstprivate(defs)
      for (eye = 0; eye < bound; eye++)
      {
        float vert_x, vert_y, vert_z,
              vert_xprime, vert_yprime, vert_zprime;
        /* offset point for calculation of d */
        vert_x = vert[eye].x + proj_len * vert[eye].xNorm;
        vert_y = vert[eye].y + proj_len * vert[eye].yNorm;
        vert_z = vert[eye].z + proj_len * vert[eye].zNorm;

        /* offset point for calculation of d' only */
        vert_xprime = vert[eye].x + ion_exc_rad * vert[eye].xNorm;
        vert_yprime = vert[eye].y + ion_exc_rad * vert[eye].yNorm;
        vert_zprime = vert[eye].z + ion_exc_rad * vert[eye].zNorm;

        pots[eye] = 0;

        for (k = 0; k < nres; k++)
        {
          int natoms = resa[k].natoms;
          int j;

          for (j = 0; j < natoms; j++)
          {
            /* distance from point to charge */
            d = dist
                (
                  vert_x, vert_y, vert_z,
                  ata[atit[k] + j].x,
                  ata[atit[k] + j].y,
                  ata[atit[k] + j].z
                );

            /* distance from point to charge */
            dprime = dist
                     (
                       vert_xprime, vert_yprime, vert_zprime,
                       ata[atit[k] + j].x,
                       ata[atit[k] + j].y,
                       ata[atit[k] + j].z
                     );


            defs.Asq_by_dsq = defs.Asq * d * d;

            pots[eye] +=
              calculate_potential_contribution
              (
                d, dprime, r, rprime,
                r0, A, sal,
                ata[atit[k] + j].charge,
                diel_int, diel_ext, &defs
              );
          }
        }
      }
    }
    else
    {
      ctsar_init(&s, bound, CTSAR_STATIC, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL, NULL, NULL);
      /* ctsar_init_memories(NULL,sizeof(float)*bound); */

      /* int out_it;
         for(out_it=0; out_it < s->d_end; out_it++)
         {
         s = ctsar_next(bound,out_it);
         } */

      #pragma omp parallel default(shared)\
      private(eye)\
      private(eye,k,d,dprime)\
      firstprivate(defs,r,rprime)
      {
        {
          int tid = omp_get_thread_num();

          do
          {
            ctsar_next(s, bound);
            atom * cata = ctsar_reg_mem(s, ata, sizeof(atom), total_atoms,
                                        CTSAR_MEM_PERSIST | CTSAR_MEM_INPUT);
            int * catit = ctsar_reg_mem(s, atit, sizeof(int), nres,
                                        CTSAR_MEM_PERSIST | CTSAR_MEM_INPUT);
            vertx * cvert = ctsar_reg_mem(s, vert, sizeof(vertx), bound,
                                          CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT);
            analytical_definitions_struct * mydefs = ctsar_reg_mem(s, &defs, sizeof(analytical_definitions_struct), 1,
                CTSAR_MEM_PERSIST | CTSAR_MEM_INPUT);
            float* data = ctsar_reg_mem(s, pots, sizeof(float), bound,
                                        CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT);
            ctsar_start(s);
            int gts = CSTART(s, tid),
                gte = CEND(s, tid);
            /* analytical_definitions_struct *mydefs = &defs; */
            /* printf("running job: %d gts: %d gte: %d\n", tid,gts,gte); */
            /*acc_set_device_num(tid, acc_device_nvidia);*/
            #pragma acc kernels loop independent if(ctsar_get_type(s) == CTSAR_DEV_GPU) deviceptr(data) deviceptr(mydefs) deviceptr(cvert,cata, catit)
            for (eye = gts; eye < gte; eye++)
            {
              float vert_x, vert_y, vert_z,
                    vert_xprime, vert_yprime, vert_zprime,
                    pot = 0;

              /* offset point for calculation of d */
              vert_x = cvert[eye].x + proj_len * cvert[eye].xNorm;
              vert_y = cvert[eye].y + proj_len * cvert[eye].yNorm;
              vert_z = cvert[eye].z + proj_len * cvert[eye].zNorm;

              /* offset point for calculation of d' only */
              vert_xprime = cvert[eye].x + ion_exc_rad * cvert[eye].xNorm;
              vert_yprime = cvert[eye].y + ion_exc_rad * cvert[eye].yNorm;
              vert_zprime = cvert[eye].z + ion_exc_rad * cvert[eye].zNorm;
              data[eye] = 0;

              int j;
#pragma acc loop independent reduction(+:pot)
              for (j = 0; j < total_atoms; j++)
              {
                /* distance from point to charge */
                d = dist
                    (
                      vert_x, vert_y, vert_z,
                      cata[j].x,
                      cata[j].y,
                      cata[j].z
                    );

                /* distance from point to charge */
                dprime = dist
                         (
                           vert_xprime, vert_yprime, vert_zprime,
                           cata[j].x,
                           cata[j].y,
                           cata[j].z
                         );


                float Asq_by_dsq = mydefs[0].Asq * d * d;
                /* local variables */
                /*******************/
                float charge = cata[j].charge;
                float sum1 = 0,
                      sum2 = 0,
                      sum3 = 0,
                      sum4 = 0,
                      salt = 0,
                      coulomb = 0,
                      to_return = 0,
                      one_over_one_plus_kappa_rprime;

                /* <><> add salt to the inside <><> */

                if (mydefs[0].phiType != TOTAL_POTENTIAL)
                {
                  coulomb = charge / (d * diel_int);
                }

                if (mydefs[0].phiType & REACTION_POTENTIAL)
                {
                  if (mydefs[0].region == 3)  /* Dist to surf > ion_excl_rad */
                  {

                    sum1 = (charge
                            * (exp (-mydefs[0].kappa * (d - dprime))) / d)
                           * (mydefs[0].one_plus_alpha / (1. + mydefs[0].kappa * dprime));

                    sum2 = charge
                           * (mydefs[0].alpha_by_one_minus_beta * exp (-mydefs[0].kappa * (r - rprime)) / (r * (1. + mydefs[0].kappa * rprime)));

                    to_return = mydefs[0].inverse_one_plus_ab_by_diel_ext * (sum1 - sum2);
                  }
                  else if (mydefs[0].region == 2) /* 0 < dist_to_surf <= ion_excl_rad */
                  {

                    one_over_one_plus_kappa_rprime = 1. / (1. + (mydefs[0].kappa * rprime));

                    /* electrostatic terms for region 2 */
                    sum1 = mydefs[0].one_plus_alpha * charge / d;
                    sum2 = mydefs[0].alpha_by_one_minus_beta * charge / r;

                    /* salt terms for region 2 */
                    sum3 = charge * (mydefs[0].alpha_by_one_minus_beta / rprime) * (1. - one_over_one_plus_kappa_rprime);
                    sum4 = mydefs[0].one_plus_alpha * charge
                           * ( 1.0 / (dprime + mydefs[0].kappa * dprime * dprime) - (1. / dprime));

                    to_return = mydefs[0].inverse_one_plus_ab_by_diel_ext * (sum1 - sum2 + sum3 + sum4);
                  }
                  else /* region = 1 */
                  {

                    one_over_one_plus_kappa_rprime = 1. / (1. + (mydefs[0].kappa * rprime));

                    sum1 =  1. / (d * diel_int);

                    sum2 =  (1. / diel_int - 1. / diel_ext) / (mydefs[0].one_plus_alpha_beta * A);

                    sum3 = mydefs[0].Asq / sqrt(mydefs[0].Asq_minus_rnautsq * mydefs[0].Asq_minus_rsq + Asq_by_dsq);

                    salt = mydefs[0].inverse_one_plus_ab_by_diel_ext
                           * (mydefs[0].one_plus_alpha / dprime * (1. / (1. + mydefs[0].kappa * dprime) - 1.)
                              + mydefs[0].alpha_by_one_minus_beta / rprime * (1. - one_over_one_plus_kappa_rprime));

                    to_return = (sum1 - sum2 * sum3 + salt) * charge;
                    /* to_return = d;//(sum1 - sum2*sum3 + salt) * atom_charge[j+atom_index[k]]; */
                  }

                  if (mydefs[0].phiType == REACTION_POTENTIAL)
                  {
                    to_return -= coulomb;
                  }
                }
                else
                {
                  to_return = coulomb;
                }


                pot += to_return;
                /* data[eye]=5; */

              }
              data[eye] = pot;
            }
            ctsar_end(s);
          }
          while (ctsar_loop(s));
        }
      }
    }
  }
  *i = eye;

  {
    int j, k;
    for (eye = 0; eye < bound; eye++)
    {
      vert[eye].potential = pots[eye];
    }
  }
  printf("time:%lu\n", get_time());
}

