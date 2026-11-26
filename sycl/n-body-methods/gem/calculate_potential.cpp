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

#include "dwarf_sycl.hpp"   /* first: the app headers define min/max macros */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>
#include <malloc.h>
#include "calculations.h"
#include "defines.h"
int gem_valid = 1;   /* set by calc_potential_single_step: all potentials finite */
#define AOCL_ALIGNMENT 64



/* ---- kernel, translated from calculate_potential.cl (body verbatim; native_exp -> sycl::native::exp).
 * OpenCL NDRange {256, BLOCKS} local {256, 1}: OpenCL dim 0 = SYCL dim 1. ---- */
#define K_ALPHA_OF .571412f
#define K_REACTION_POTENTIAL 2
#define K_TOTAL_POTENTIAL 3
namespace dev {
typedef struct
{
	int phiType,
	    region;

	float one_plus_alpha,
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


//extern int printf(constant char* format, ...);

inline float dist2D(float x, float y, float z, float x2, float y2, float z2)
{
	
	float xd = x - x2,
	      yd = y - y2,
	      zd = z - z2;
	return (xd*xd + yd*yd + zd*zd);
} 

inline float distD (float x, float y, float z, float x2, float y2, float z2)
{
	return sycl::sqrt(dist2D(x, y, z, x2, y2, z2));
} 


} // namespace dev

static void calc_potential_single_step_dev(sycl::queue &q, sycl::nd_range<2> ndr, int nres, int nvert, float A, float proj_len, float diel_int, float diel_ext, float sal, float ion_exc_rad, int phiType, int eye_arg, int step_size, unsigned int * atom_addrs, unsigned int * atom_lengths, float r, float r0_arg, float rprime, float *res_c_s, float *res_x_s, float *res_y_s, float *res_z_s, float *at_c_s, float *at_x_s, float *at_y_s, float *at_z_s, float *vert_c_s, float *vert_x_s, float *vert_y_s, float *vert_z_s, float *vert_x_p_s, float *vert_y_p_s, float *vert_z_p_s)
{
	q.parallel_for<class gem_potential>(ndr, [=](sycl::nd_item<2> item) {
		int eye = eye_arg; float r0 = r0_arg;   /* private copies of the arguments the body modifies */
	
	
		int   j,
		      k,
		      //bound,
		      natoms;
		float d,
		      dprime,
		      charge;
		float sum1 = 0,
		      sum2 = 0,
		      sum3 = 0,
		      sum4 = 0,
		      salt = 0,
		      coulomb = 0,
		      to_return = 0,
		      one_over_one_plus_kappa_rprime;
		float potential = 0;



		dev::analytical_definitions_struct defs;
		//calculate bound
		//bound = eye+step_size < nvert ? nvert : eye + step_size;
		//calculate the vertx on which to act
		eye = eye + (item.get_global_range(1)*item.get_global_range(0)*0 /* get_global_id(2) of a 2-D range */)+(item.get_global_range(1)*item.get_global_id(0))+item.get_global_id(1);
		if(eye >= nvert)//TODO:double check this
			return;

		//this shouldn't be here, but the kernel wigs out without it, think there might be a maximum struct size for arguments... must look into this
		defs.Asq   = A*A;
		defs.Asq_minus_rsq = defs.Asq - (r*r);

		defs.kappa = 0.316f * sycl::sqrt(sal);
		defs.beta = diel_int/diel_ext;
		defs.alpha_by_one_minus_beta = K_ALPHA_OF*(1.0f - defs.beta);
		defs.alpha_beta = K_ALPHA_OF * defs.beta;
		defs.one_plus_alpha_beta = 1.0f + defs.alpha_beta;
		defs.one_plus_alpha = 1.0f + K_ALPHA_OF;
		defs.inverse_one_plus_ab_by_diel_ext = 1.0f/(defs.one_plus_alpha_beta*diel_ext);
		defs.phiType = phiType;

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

	
	
		r0 = A/2.0f;  
		vert_c_s[eye] = 0;
		defs.Asq_minus_rnautsq = defs.Asq - r0*r0;
		//for (; eye < bound; eye++)
		{
			//vert[eye].potential = 0;

			for (k = 0; k < nres; k++)
			{
				natoms = atom_lengths[k];

				for (j = 0; j < natoms; j++)
				{
				
					d = dev::distD
						(
						 vert_x_s[eye], vert_y_s[eye], vert_z_s[eye],
						 at_x_s[atom_addrs[k]+j],
						 at_y_s[atom_addrs[k]+j],
						 at_z_s[atom_addrs[k]+j]
						);

				
					dprime = dev::distD
						(
						 vert_x_p_s[eye], vert_y_p_s[eye], vert_z_p_s[eye],
						 at_x_s[atom_addrs[k]+j],
						 at_y_s[atom_addrs[k]+j],
						 at_z_s[atom_addrs[k]+j]
						);


					defs.Asq_by_dsq = A*A * d*d;

					charge = at_c_s[atom_addrs[k]+j];
				
				

				
					if (defs.phiType != K_TOTAL_POTENTIAL)
					{
						coulomb = charge/(d * diel_int);
					}


					if (defs.phiType & K_REACTION_POTENTIAL)
					{
						if (defs.region == 3)  
						{

							sum1 = (charge
									* (sycl::native::exp (-defs.kappa * (d - dprime))) / d)
								* (defs.one_plus_alpha / (1.0f + defs.kappa * dprime));

							sum2 = charge
								* (defs.alpha_by_one_minus_beta * sycl::native::exp (-defs.kappa * (r-rprime)) / (r * (1.0f + defs.kappa*rprime)));

							to_return = defs.inverse_one_plus_ab_by_diel_ext * (sum1 - sum2);
						}
						else if (defs.region == 2) 
						{

							one_over_one_plus_kappa_rprime = 1.0f/(1.0f + (defs.kappa * rprime));

						
							sum1 = defs.one_plus_alpha * charge / d;
							sum2 = defs.alpha_by_one_minus_beta * charge/r;

						
							sum3 = charge * (defs.alpha_by_one_minus_beta / rprime) * (1.0f - one_over_one_plus_kappa_rprime);
							sum4 = defs.one_plus_alpha * charge
								* ( 1.0f / (dprime + defs.kappa * dprime * dprime) - (1.0f / dprime));

							to_return = defs.inverse_one_plus_ab_by_diel_ext * (sum1 - sum2 + sum3 + sum4);
						}
						else 
						{

							one_over_one_plus_kappa_rprime = 1.0f/(1.0f + (defs.kappa * rprime));

							sum1 =  1.0f / (d * diel_int);

							sum2 =  (1.0f/diel_int - 1.0f/diel_ext) / (defs.one_plus_alpha_beta * A);

							sum3 = defs.Asq / sycl::sqrt(defs.Asq_minus_rnautsq*defs.Asq_minus_rsq + defs.Asq_by_dsq);

							salt = defs.inverse_one_plus_ab_by_diel_ext
								* (defs.one_plus_alpha/dprime*(1.0f/(1.0f+defs.kappa*dprime)-1.0f)  
										+defs.alpha_by_one_minus_beta/rprime*(1.0f-one_over_one_plus_kappa_rprime));

							to_return = (sum1 - sum2*sum3 + salt) * charge;
						}

						if (defs.phiType == K_REACTION_POTENTIAL)
						{
							to_return -= coulomb;
						}
					}
					else
					{
						to_return = coulomb;
					}

					//potential += to_return;
					vert_c_s[eye] += to_return;
					//vert_c_s[eye] = 5;
				}
			}
		}
		//vert_c_s[eye] = potential;
}).wait_and_throw();
}

typedef long long int64;
static struct timeval start_time;

void init_time() {
	gettimeofday(&start_time, NULL);
}

int64 get_time() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return (int64) (t.tv_sec - start_time.tv_sec) * 1000000
		+ (t.tv_usec - start_time.tv_usec);
}
// #define BLOCK_DIM_X 8
// #define BLOCK_DIM_Y 8
int  BLOCK_DIM_X=16;
int  BLOCK_DIM_Y=16;
//nthreads/block 64
//#define BLOCKS 390
#define BLOCKS 120
//multiple of 30

/* Calculates potential for surface out of a molecule -- this is the free
 * space solution. --JCG */
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
		int step_size)
{
	int it, eye, bound;
	//residue *residues_s;
	float *res_c, *res_x, *res_y, *res_z,
	      *at_c, *at_x, *at_y, *at_z,
	      *vert_c, *vert_x, *vert_y, *vert_z,
	      *vert_x_p, *vert_y_p, *vert_z_p;
	float *res_c_s, *res_x_s, *res_y_s, *res_z_s,
	       *at_c_s, *at_x_s, *at_y_s, *at_z_s,
	       *vert_c_s, *vert_x_s, *vert_y_s, *vert_z_s,
	       *vert_x_p_s, *vert_y_p_s, *vert_z_p_s;
	unsigned int *atom_addrs_s;
	unsigned int *atom_lengths_s;
	int natoms;
	unsigned int *atom_addrs;
	unsigned int *atom_lengths;

	static sycl::queue q = dsycl::make_queue();
	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */

	init_time();

	eye = *i;

	bound = eye + step_size;

	if (nvert < bound) bound = nvert;

	natoms = 0;
	for (it = 0; it < nres; it++)
	{
		natoms += residues[it].natoms;
	}

	    res_c = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nres);
	    res_x = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nres);
	    res_y = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nres);
	    res_z = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nres);
	    at_c = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * natoms);
	    at_x = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * natoms);
	    at_y = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * natoms);
	    at_z = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * natoms);
	    vert_c = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nvert);
	    vert_x = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nvert);
	    vert_y = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nvert);
	    vert_z = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nvert);
	    vert_x_p = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nvert);
	    vert_y_p = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nvert);
	    vert_z_p = (float *)  memalign ( AOCL_ALIGNMENT,sizeof(float) * nvert);
	    /* allocate temporary arrays for atoms and start addresses */
	    atom_addrs = (unsigned int *) memalign ( AOCL_ALIGNMENT,nres * sizeof(unsigned int));
	    atom_lengths = (unsigned int *) memalign ( AOCL_ALIGNMENT,nres * sizeof(unsigned int));

        //res_c = (float *) malloc(sizeof(float) * nres);
	    //res_x = (float *) malloc(sizeof(float) * nres);
	    //res_y = (float *) malloc(sizeof(float) * nres);
	    //res_z = (float *) malloc(sizeof(float) * nres);
	    //at_c = (float *) malloc(sizeof(float) * natoms);
	    //at_x = (float *) malloc(sizeof(float) * natoms);
	    //at_y = (float *) malloc(sizeof(float) * natoms);
	    //at_z = (float *) malloc(sizeof(float) * natoms);
	    //vert_c = (float *) malloc(sizeof(float) * nvert);
	    //vert_x = (float *) malloc(sizeof(float) * nvert);
	    //vert_y = (float *) malloc(sizeof(float) * nvert);
	    //vert_z = (float *) malloc(sizeof(float) * nvert);
	    //vert_x_p = (float *) malloc(sizeof(float) * nvert);
	    //vert_y_p = (float *) malloc(sizeof(float) * nvert);
	    //vert_z_p = (float *) malloc(sizeof(float) * nvert);
	    //
	    //atom_addrs = (unsigned int *)malloc(nres * sizeof(unsigned int));
	    //atom_lengths = (unsigned int *)malloc(nres * sizeof(unsigned int));


	/* copy the atom structure to the device */
	natoms = 0;
	for (it = 0; it < nres; it++)
	{
		atom_addrs[it] = natoms;

		res_c[it] = residues[it].total_charge;
		res_x[it] = residues[it].x;
		res_y[it] = residues[it].y;
		res_z[it] = residues[it].z;

		for(int j=0; j<residues[it].natoms; j++)
		{
			at_c[natoms+j] = residues[it].atoms[j].charge;
			at_x[natoms+j] = residues[it].atoms[j].x;
			at_y[natoms+j] = residues[it].atoms[j].y;
			at_z[natoms+j] = residues[it].atoms[j].z;
		}
		natoms += residues[it].natoms;
		atom_lengths[it] = residues[it].natoms;
	}

	// free(res_c);
	// free(res_x);
	// free(res_y);
	// free(res_z);
	// free(at_c );
	// free(at_x );
	// free(at_y );
	// free(at_z );


	for (it = 0; it < nvert; it++)
	{
		vert_c[it] = 0;
		/* offset point for calculation of d */
		vert_x[it] = vert[it].x + proj_len * vert[it].xNorm;
		vert_y[it] = vert[it].y + proj_len * vert[it].yNorm;
		vert_z[it] = vert[it].z + proj_len * vert[it].zNorm;

		/* offset point for calculation of d' only */
		vert_x_p[it] = vert[it].x + ion_exc_rad * vert[it].xNorm;
		vert_y_p[it] = vert[it].y + ion_exc_rad * vert[it].yNorm;
		vert_z_p[it] = vert[it].z + ion_exc_rad * vert[it].zNorm;
	}


	//TODO: fix this
	//free(vert_x);
	//free(vert_y);
	//free(vert_z);
	//free(vert_x_p);
	//free(vert_y_p);
	//free(vert_z_p);

	/* free temporary arrays for atoms and start addresses */
	// free(atom_addrs);
	// free(atom_lengths);


	float r, r0, rprime;
	r = A + proj_len;
	rprime = A + ion_exc_rad;


	/* need to estimate r0, dist from atom to surface */
	/* if we intend to support negative projection lengths <><> */
	r0 = A/2.;  /* might be tricky ... hrmmm */

	res_c_s = dsycl::alloc<float>(q, nres);
	res_x_s = dsycl::alloc<float>(q, nres);
	res_y_s = dsycl::alloc<float>(q, nres);
	res_z_s = dsycl::alloc<float>(q, nres);
	at_c_s = dsycl::alloc<float>(q, natoms);
	at_x_s = dsycl::alloc<float>(q, natoms);
	at_y_s = dsycl::alloc<float>(q, natoms);
	at_z_s = dsycl::alloc<float>(q, natoms);
	vert_c_s = dsycl::alloc<float>(q, nvert);
	vert_x_s = dsycl::alloc<float>(q, nvert);
	vert_y_s = dsycl::alloc<float>(q, nvert);
	vert_z_s = dsycl::alloc<float>(q, nvert);
	vert_x_p_s = dsycl::alloc<float>(q, nvert);
	vert_y_p_s = dsycl::alloc<float>(q, nvert);
	vert_z_p_s = dsycl::alloc<float>(q, nvert);
	atom_addrs_s = dsycl::alloc<unsigned int>(q, nres);
	atom_lengths_s = dsycl::alloc<unsigned int>(q, nres);


	dsycl::write(q, res_c_s, res_c, sizeof(cl_float)*nres);
	dsycl::write(q, res_x_s, res_x, sizeof(cl_float)*nres);
	dsycl::write(q, res_y_s, res_y, sizeof(cl_float)*nres);
	dsycl::write(q, res_z_s, res_z, sizeof(cl_float)*nres);
	dsycl::write(q, at_c_s, at_c, sizeof(cl_float)*natoms);
	dsycl::write(q, at_x_s, at_x, sizeof(cl_float)*natoms);
	dsycl::write(q, at_y_s, at_y, sizeof(cl_float)*natoms);
	dsycl::write(q, at_z_s, at_z, sizeof(cl_float)*natoms);
	dsycl::write(q, vert_c_s, vert_c, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_x_s, vert_x, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_y_s, vert_y, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_z_s, vert_z, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_x_p_s, vert_x_p, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_y_p_s, vert_y_p, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_z_p_s, vert_z_p, sizeof(cl_float)*nvert);
	dsycl::write(q, atom_addrs_s, atom_addrs, sizeof(cl_int)*nres);
	dsycl::write(q, atom_lengths_s, atom_lengths, sizeof(cl_int)*nres);

	/* [FIX-OCL-34] pinned launch: work-group = DWARF_LS = 256, BLOCKS=120 work-groups per launch; the kernel guards eye >= nvert */
	BLOCK_DIM_X = 16; BLOCK_DIM_Y = 16;
	const int arg_step_size = step_size;   /* kernel argument 10 is set once, before the loop */
	sycl::nd_range<2> ndr({BLOCKS, DWARF_WG_ALG}, {1, DWARF_WG_ALG});

	/* common timer: repetitions around [all buffer uploads, the vertex-block kernel loop, the download of the potentials];
	 * the uploads in front of the kernel-argument setup above are the first (untimed) copy of the same data. */
	dwarf_init("gem", "n-body-methods", "sycl");
	dwarf_set_size(nvert, natoms);
	int eye0 = eye, step0 = step_size, rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	eye = eye0; step_size = step0;
	dwarf_region_begin();
	dsycl::write(q, res_c_s, res_c, sizeof(cl_float)*nres);
	dsycl::write(q, res_x_s, res_x, sizeof(cl_float)*nres);
	dsycl::write(q, res_y_s, res_y, sizeof(cl_float)*nres);
	dsycl::write(q, res_z_s, res_z, sizeof(cl_float)*nres);
	dsycl::write(q, at_c_s, at_c, sizeof(cl_float)*natoms);
	dsycl::write(q, at_x_s, at_x, sizeof(cl_float)*natoms);
	dsycl::write(q, at_y_s, at_y, sizeof(cl_float)*natoms);
	dsycl::write(q, at_z_s, at_z, sizeof(cl_float)*natoms);
	dsycl::write(q, vert_c_s, vert_c, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_x_s, vert_x, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_y_s, vert_y, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_z_s, vert_z, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_x_p_s, vert_x_p, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_y_p_s, vert_y_p, sizeof(cl_float)*nvert);
	dsycl::write(q, vert_z_p_s, vert_z_p, sizeof(cl_float)*nvert);
	dsycl::write(q, atom_addrs_s, atom_addrs, sizeof(cl_int)*nres);
	dsycl::write(q, atom_lengths_s, atom_lengths, sizeof(cl_int)*nres);
	for (;  eye/*+((BLOCK_DIM_X*BLOCK_DIM_Y)*BLOCKS)*/ < bound; eye+=((BLOCK_DIM_X*BLOCK_DIM_Y)*BLOCKS))
		//for(it=0; it<3; it++)
	{
		fprintf(stdout,"finished first %d of %d\n", eye, nvert);
		calc_potential_single_step_dev(q, ndr, nres, nvert, A, proj_len, diel_int, diel_ext, sal, ion_exc_rad, phiType, eye, arg_step_size,
				atom_addrs_s, atom_lengths_s, r, r0, rprime,
				res_c_s, res_x_s, res_y_s, res_z_s, at_c_s, at_x_s, at_y_s, at_z_s,
				vert_c_s, vert_x_s, vert_y_s, vert_z_s, vert_x_p_s, vert_y_p_s, vert_z_p_s);

		if(eye + step_size > bound)
			step_size = bound - eye;


		/* copy back the calculation result */
	}

	if(eye > bound)
		eye = bound;
	dsycl::read(q, vert_c, vert_c_s, nvert * sizeof(cl_float));

	q.wait_and_throw();
	dwarf_region_end(rep);
	}
	printf("runtime:%llu\n", get_time());


	gem_valid = 1;
	for(it = 0; it < nvert; it++)
	{
		vert[it].potential = vert_c[it];
		if(!(vert_c[it] == vert_c[it]) || vert_c[it] > 1e30f || vert_c[it] < -1e30f) gem_valid = 0;   /* [FIX-OCL-23] no validity check before: NaN/Inf potentials */
	}

	*i = eye;
}

