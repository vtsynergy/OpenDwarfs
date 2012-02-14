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

#ifndef __CALCULATIONS_H__
#define __CALCULATIONS_H__

#include "structures.h"
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

/* definitions for operations that project_grd can do */
#define DIFFERENTIATE (uchar)1
#define SET_VALUES    (uchar)2

/* some useful typedefs for shortening function definitions */

typedef struct
{
    cl_int phiType,
        region;

    cl_float one_plus_alpha,
           beta,
           alpha_beta,
           one_plus_alpha_beta,
           alpha_by_one_minus_beta,
           inverse_one_plus_ab_by_diel_ext,
           kappa,
           Asq,
           Asq_minus_rnautsq,
           Asq_minus_rsq,
           Asq_by_dsq;
           
} analytical_definitions_struct;


#ifdef __cplusplus
extern "C" {
#endif

/* calculates potential for a point */
/* void calc_potential_single_step (residue*, int, vertx*, int, double, double, double, double, double, double, int, int *, int); */
void calc_potential_single_step(residue *residues, 
                                int nres, 
                                vertx *vert, 
                                int nvert, 
                                float A, 
                                float proj_len, 
                                float diel_int, 
                                float diel_ext, 
                                float sal, 
                                float ion_exc_rad, 
                                int phiType, 
                                int *i, 
                                int step_size);

/* calculates potential for a grid of points */
void calc_potential_grid(residue*,int,vertx*,int,int*,float*,int,int,int,float[6],double,double,double,double,double,int);


/* finds bonds in the atoms */
int extrapolate_bonds (residue*, int, bond **);

/* estimates the "floating center" within a range */
vertx floating_center (double, vertx, residue *, int);

/* interpolates precalculated potentials in phimap format */
void project_grd(vertx *, float *, int, int, int, int, float[6], float max_rng, unsigned char op);

/* resamples a grid of dimension n1 to a grid of dimension n2 */
void resample_grid(float *, int, float[6], float *, int, float[6]);

/* samples the potential at a given physical location according to the given grid */
inline float sample_grid (float, float, float, float *, int, int, int, float[6]);

/* vector math stuff */
/*********************/

/* dot product for 2 3-d vectors */
double dot_prod   (double v1[3], double v2[3]);

/* magnitude of a vector (dist of vector with origin */
double magnitude  (double v[3]);

/* cross product between two vectors, gets stored in cross */
void   cross_prod (double v1[3], double v2[3], double cross[3]);

/* distance between 2 3d points */
inline double dist (double x, double y, double z, double x2, double y2, double z2);

/* distance between 2 3d points squared */
inline double dist2(double x, double y, double z, double x2, double y2, double z2);



/* gets the dimensions of a set of vertices in space */
int get_dimensions(vertx *vert, int nvert, float[4], float[4], float[4]);

/* sets r, g, b correctly for val in range */
void getColor (float *r, float *g, float *b, float val, float range[4]);

/* sets r, g, b correctly for a specific atom type */
void getAtomTypeColor (float *r, float *g, float *b, char *name);

/* estimates the electrostatic radius of the given molecule */
double estimate_A(residue *, int);

/* uses radix sort to sort a given array of objects */
void radix_sort (ushort keylen, ushort keyloc, ushort size, ulong N, void *sort);

/* determines minimum distance from grid point to mol surface */
void gridDistances (float** grid,
                    int** surfaceIdx,
                    const int gridxDim, 
                    const int gridyDim, 
                    const int gridzDim,
                    float gridWorldDims[6],
                    const vertx* molVerts,
                    const int numMolVerts);

/* this function calculates the surface charge */
int calculate_surface_charge(vertx *, int, float *, int, float, float[6], float);

void calc_charge_potential_single_step(residue *, int, vertx *, int, double, double, double, double, double, double, int, int *, int);

void calc_charge_single_step(residue *, int, vertx *, int, double, double, double, double, double, double, int, int *, int);
   
/* this function tells us the endian-ness of our system */
inline uchar big_endian (void);

int extract_statistics (residue *, int, vertx *, int, float[3], float[3], float[3], float[3], float[3]);

/************/
/* MACROS   */
/************/

/* returns the minimum, a or b */
#define min(a, b) (a < b)?a:b

/* returns the maximum, a or b */
#define max(a, b) (a > b)?a:b

/* clamps a between b and c */
#define clamp(a, b, c) (max(b, (min(a, c))))

#ifdef __cplusplus
}
#endif

#endif
