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

#include "calculations.h"
#include "defines.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>
#include <malloc.h>

#include "dwarf_cuda.h"
#include "dwarf_timer.h"
int gem_valid = 1;   /* set by calc_potential_single_step: all potentials finite */
#define AOCL_ALIGNMENT 64

#include "calculate_potential_kernel.cuh"

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


/*
 * Converts the contents of a file into a string
 */
	std::string
convertToString(const char *filename)
{
	size_t size;
	char*  str;
	std::string s;

	std::fstream f(filename, (std::fstream::in | std::fstream::binary));

	if(f.is_open())
	{
		size_t fileSize;
		f.seekg(0, std::fstream::end);
		size = fileSize = f.tellg();
		f.seekg(0, std::fstream::beg);

		str = new char[size+1];
		if(!str)
		{
			f.close();
			return NULL;
		}

		f.read(str, fileSize);
		f.close();
		str[size] = '\0';

		s = str;

		return s;
	}
	return NULL;
}


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
#define GEM_H2D(dst, src, bytes) do { CUDA_CHECK(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice)); CUDA_CHECK(cudaDeviceSynchronize()); } while (0)
#define GEM_UPLOAD_ALL() do { \
	GEM_H2D(res_c_s, res_c, sizeof(float)*nres); GEM_H2D(res_x_s, res_x, sizeof(float)*nres); \
	GEM_H2D(res_y_s, res_y, sizeof(float)*nres); GEM_H2D(res_z_s, res_z, sizeof(float)*nres); \
	GEM_H2D(at_c_s, at_c, sizeof(float)*natoms); GEM_H2D(at_x_s, at_x, sizeof(float)*natoms); \
	GEM_H2D(at_y_s, at_y, sizeof(float)*natoms); GEM_H2D(at_z_s, at_z, sizeof(float)*natoms); \
	GEM_H2D(vert_c_s, vert_c, sizeof(float)*nvert); GEM_H2D(vert_x_s, vert_x, sizeof(float)*nvert); \
	GEM_H2D(vert_y_s, vert_y, sizeof(float)*nvert); GEM_H2D(vert_z_s, vert_z, sizeof(float)*nvert); \
	GEM_H2D(vert_x_p_s, vert_x_p, sizeof(float)*nvert); GEM_H2D(vert_y_p_s, vert_y_p, sizeof(float)*nvert); \
	GEM_H2D(vert_z_p_s, vert_z_p, sizeof(float)*nvert); \
	GEM_H2D(atom_addrs_s, atom_addrs, sizeof(int)*nres); GEM_H2D(atom_lengths_s, atom_lengths, sizeof(int)*nres); } while (0)
	/* clCreateBuffer */
	CUDA_CHECK(cudaMalloc((void**)&res_c_s, sizeof(float)*nres));
	CUDA_CHECK(cudaMalloc((void**)&res_x_s, sizeof(float)*nres));
	CUDA_CHECK(cudaMalloc((void**)&res_y_s, sizeof(float)*nres));
	CUDA_CHECK(cudaMalloc((void**)&res_z_s, sizeof(float)*nres));
	CUDA_CHECK(cudaMalloc((void**)&at_c_s, sizeof(float)*natoms));
	CUDA_CHECK(cudaMalloc((void**)&at_x_s, sizeof(float)*natoms));
	CUDA_CHECK(cudaMalloc((void**)&at_y_s, sizeof(float)*natoms));
	CUDA_CHECK(cudaMalloc((void**)&at_z_s, sizeof(float)*natoms));
	CUDA_CHECK(cudaMalloc((void**)&vert_c_s, sizeof(float)*nvert));
	CUDA_CHECK(cudaMalloc((void**)&vert_x_s, sizeof(float)*nvert));
	CUDA_CHECK(cudaMalloc((void**)&vert_y_s, sizeof(float)*nvert));
	CUDA_CHECK(cudaMalloc((void**)&vert_z_s, sizeof(float)*nvert));
	CUDA_CHECK(cudaMalloc((void**)&vert_x_p_s, sizeof(float)*nvert));
	CUDA_CHECK(cudaMalloc((void**)&vert_y_p_s, sizeof(float)*nvert));
	CUDA_CHECK(cudaMalloc((void**)&vert_z_p_s, sizeof(float)*nvert));
	CUDA_CHECK(cudaMalloc((void**)&atom_addrs_s, sizeof(float)*nres));
	CUDA_CHECK(cudaMalloc((void**)&atom_lengths_s, sizeof(float)*nres));

	GEM_UPLOAD_ALL();

	/* [FIX-OCL-34] pinned launch: work-group = DWARF_LS = 256 (16x16 work-items per group, 1-D in effect), BLOCKS=120 work-groups per
	 * launch (= OpenMP num_teams(120) thread_limit(256)); was halved according to a device query. The kernel guards eye >= nvert. */
	BLOCK_DIM_X = 16; BLOCK_DIM_Y = 16;


	/* OpenCL NDRange {256, BLOCKS} with local {256, 1} -> grid (1, BLOCKS) of 256-thread blocks */
	dim3 gridDims(1, BLOCKS);
	dim3 blockDims(DWARF_WG_ALG, 1);


	/* common timer: repetitions around [all buffer uploads, the vertex-block kernel loop, the download of the potentials];
	 * the uploads in front of the kernel-argument setup above are the first (untimed) copy of the same data. */
	dwarf_init("gem", "n-body-methods", "cuda");
	dwarf_set_size(nvert, natoms);
	int eye0 = eye, step0 = step_size, rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	eye = eye0; step_size = step0;
	dwarf_region_begin();
	GEM_UPLOAD_ALL();
	for (;  eye/*+((BLOCK_DIM_X*BLOCK_DIM_Y)*BLOCKS)*/ < bound; eye+=((BLOCK_DIM_X*BLOCK_DIM_Y)*BLOCKS))
		//for(it=0; it<3; it++)
	{
		fprintf(stdout,"finished first %d of %d\n", eye, nvert);
		/* copy the vert data of the current block to device */
		calc_potential_single_step_dev<<<gridDims, blockDims>>>(nres, nvert, (float)A, (float)proj_len, (float)diel_int, (float)diel_ext, (float)sal,
				(float)ion_exc_rad, phiType, eye, step_size, atom_addrs_s, atom_lengths_s, r, r0, rprime,
				res_c_s, res_x_s, res_y_s, res_z_s, at_c_s, at_x_s, at_y_s, at_z_s,
				vert_c_s, vert_x_s, vert_y_s, vert_z_s, vert_x_p_s, vert_y_p_s, vert_z_p_s);
		/* wait for the kernel call to finish execution */
		CUDA_SYNC_LAUNCH();

		if(eye + step_size > bound)
			step_size = bound - eye;


		/* copy back the calculation result */
	}

	if(eye > bound)
		eye = bound;
	CUDA_CHECK(cudaMemcpy(vert_c, vert_c_s, nvert * sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());
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

