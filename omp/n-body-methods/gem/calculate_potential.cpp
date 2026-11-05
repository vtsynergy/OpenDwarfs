/* GEM electrostatic potential (n-body dwarf) -- OpenMP target-offload translation of
 * opencl/n-body-methods/gem/calculate_potential.cpp + calculate_potential.cl.
 *
 * calc_potential_single_step() keeps the OpenCL host logic: flatten the residues/atoms/vertices into
 * structure-of-arrays buffers, upload them, process the vertices in blocks of BLOCK_DIM_X*BLOCK_DIM_Y*BLOCKS
 * (= 256 x 120 = 30720) work-items per kernel launch, read the potentials back. The kernel body
 * (calc_potential_single_step_dev) is the OpenCL kernel text: one work-item per vertex, a double loop over all
 * residues/atoms. All the other GEM sources (surface generation, I/O, statistics) are the unchanged host code.
 * Timed region: device allocation + uploads + all kernel launches + the download of the vertex potentials.
 * Validation (added): all potentials must be finite; fidelity to the OpenCL app is checked through the printed
 * RMSD / potential minimum / maximum and the vertex dump (-dumpVertices).
 */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>
#include <malloc.h>
#include <omp.h>
#include "dwarf_timer.h"
#include "calculations.h"
#include "defines.h"

int gem_valid = 1;

#define BLOCK_DIM_X_C 16
#define BLOCK_DIM_Y_C 16
#define BLOCKS 120   /* multiple of 30 */

#pragma omp declare target
#define ALPHA_OF .571412f /* "Incorporating variable dielectric environments into the generalized Born model" Sigalov, Scheffe, Onufriev J. Chem. Phys. */
#define REACTION_POTENTIAL 2
#define TOTAL_POTENTIAL 3



inline float dist2D(float x, float y, float z, float x2, float y2, float z2)
{
	/* local variables */
	float xd = x - x2,
	      yd = y - y2,
	      zd = z - z2;
	return (xd*xd + yd*yd + zd*zd);
} /* end function dist2 */

inline float distD (float x, float y, float z, float x2, float y2, float z2)
{
	return sqrtf(dist2D(x, y, z, x2, y2, z2));
} /* end function dist */

/* Calculates potential for surface out of a molecule -- this is the free
 * space solution. --JCG */

#pragma omp end declare target

static void run_kernel_block(int nres, int nvert, float A, float proj_len, float diel_int, float diel_ext, float sal,
                             float ion_exc_rad, int phiType, int eye_base, int step_size, unsigned int *atom_addrs,
                             unsigned int *atom_lengths, float r, float r0, float rprime, float *res_c_s, float *res_x_s,
                             float *res_y_s, float *res_z_s, float *at_c_s, float *at_x_s, float *at_y_s, float *at_z_s,
                             float *vert_c_s, float *vert_x_s, float *vert_y_s, float *vert_z_s, float *vert_x_p_s,
                             float *vert_y_p_s, float *vert_z_p_s, int natoms) {
  const int gsize0 = BLOCK_DIM_X_C * BLOCK_DIM_Y_C, gsize1 = BLOCKS;
#pragma omp target teams distribute parallel for collapse(2) map(present, alloc : atom_addrs[0 : nres], atom_lengths[0 : nres], at_c_s[0 : natoms], at_x_s[0 : natoms], at_y_s[0 : natoms], at_z_s[0 : natoms], vert_c_s[0 : nvert], vert_x_s[0 : nvert], vert_y_s[0 : nvert], vert_z_s[0 : nvert], vert_x_p_s[0 : nvert], vert_y_p_s[0 : nvert], vert_z_p_s[0 : nvert])
  for (int gid1 = 0; gid1 < gsize1; gid1++) {
    for (int gid0 = 0; gid0 < gsize0; gid0++) {
      int eye = eye_base + gsize0 * gid1 + gid0;
      if (eye < nvert) {
        /* ---- calc_potential_single_step_dev, body of the OpenCL kernel ---- */
        int j, k, natoms_k;
        float d, dprime, charge;
        float sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, salt = 0, coulomb = 0, to_return = 0, one_over_one_plus_kappa_rprime;
        analytical_definitions_struct defs;
        float r0_l = r0;
        (void)step_size;
        (void)A;
        {

	/* local variables */
	/*******************/
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



	analytical_definitions_struct defs;
	//calculate bound
	//bound = eye+step_size < nvert ? nvert : eye + step_size;
	//calculate the vertx on which to act
	/* eye: global vertex index, computed by the caller as in the OpenCL kernel
	 * (eye_base + global_size0*global_id1 + global_id0) and already checked against nvert */

	//this shouldn't be here, but the kernel wigs out without it, think there might be a maximum struct size for arguments... must look into this
	defs.Asq   = A*A;
	defs.Asq_minus_rsq = defs.Asq - (r*r);

	defs.kappa = 0.316f * sqrtf(sal);
	defs.beta = diel_int/diel_ext;
	defs.alpha_by_one_minus_beta = ALPHA_OF*(1.0f - defs.beta);
	defs.alpha_beta = ALPHA_OF * defs.beta;
	defs.one_plus_alpha_beta = 1.0f + defs.alpha_beta;
	defs.one_plus_alpha = 1.0f + ALPHA_OF;
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

	/* need to estimate r0, dist from atom to surface */
	/* if we intend to support negative projection lengths <><> */
	r0_l = A/2.0f;  /* might be tricky ... hrmmm */
	vert_c_s[eye] = 0;
	defs.Asq_minus_rnautsq = defs.Asq - r0_l*r0_l;
	//for (; eye < bound; eye++)
	{
		//vert[eye].potential = 0;

		for (k = 0; k < nres; k++)
		{
			natoms_k = atom_lengths[k];

			for (j = 0; j < natoms_k; j++)
			{
				/* distance from point to charge */
				d = distD
					(
					 vert_x_s[eye], vert_y_s[eye], vert_z_s[eye],
					 at_x_s[atom_addrs[k]+j],
					 at_y_s[atom_addrs[k]+j],
					 at_z_s[atom_addrs[k]+j]
					);

				/* distance from point to charge */
				dprime = distD
					(
					 vert_x_p_s[eye], vert_y_p_s[eye], vert_z_p_s[eye],
					 at_x_s[atom_addrs[k]+j],
					 at_y_s[atom_addrs[k]+j],
					 at_z_s[atom_addrs[k]+j]
					);


				defs.Asq_by_dsq = A*A * d*d;

				charge = at_c_s[atom_addrs[k]+j];
				/* local variables */
				/*******************/

				/* <><> add salt to the inside <><> */
				if (defs.phiType != TOTAL_POTENTIAL)
				{
					coulomb = charge/(d * diel_int);
				}


				if (defs.phiType & REACTION_POTENTIAL)
				{
					if (defs.region == 3)  /* Dist to surf > ion_excl_rad */
					{

						sum1 = (charge
								* (expf (-defs.kappa * (d - dprime))) / d)
							* (defs.one_plus_alpha / (1.0f + defs.kappa * dprime));

						sum2 = charge
							* (defs.alpha_by_one_minus_beta * expf (-defs.kappa * (r-rprime)) / (r * (1.0f + defs.kappa*rprime)));

						to_return = defs.inverse_one_plus_ab_by_diel_ext * (sum1 - sum2);
					}
					else if (defs.region == 2) /* 0 < dist_to_surf <= ion_excl_rad */
					{

						one_over_one_plus_kappa_rprime = 1.0f/(1.0f + (defs.kappa * rprime));

						/* electrostatic terms for region 2 */
						sum1 = defs.one_plus_alpha * charge / d;
						sum2 = defs.alpha_by_one_minus_beta * charge/r;

						/* salt terms for region 2 */
						sum3 = charge * (defs.alpha_by_one_minus_beta / rprime) * (1.0f - one_over_one_plus_kappa_rprime);
						sum4 = defs.one_plus_alpha * charge
							* ( 1.0f / (dprime + defs.kappa * dprime * dprime) - (1.0f / dprime));

						to_return = defs.inverse_one_plus_ab_by_diel_ext * (sum1 - sum2 + sum3 + sum4);
					}
					else /* region = 1 */
					{

						one_over_one_plus_kappa_rprime = 1.0f/(1.0f + (defs.kappa * rprime));

						sum1 =  1.0f / (d * diel_int);

						sum2 =  (1.0f/diel_int - 1.0f/diel_ext) / (defs.one_plus_alpha_beta * A);

						sum3 = defs.Asq / sqrtf(defs.Asq_minus_rnautsq*defs.Asq_minus_rsq + defs.Asq_by_dsq);

						salt = defs.inverse_one_plus_ab_by_diel_ext
							* (defs.one_plus_alpha/dprime*(1.0f/(1.0f+defs.kappa*dprime)-1.0f)  
									+defs.alpha_by_one_minus_beta/rprime*(1.0f-one_over_one_plus_kappa_rprime));

						to_return = (sum1 - sum2*sum3 + salt) * charge;
					}

					if (defs.phiType == REACTION_POTENTIAL)
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

        }
      }
    }
  }
}

void calc_potential_single_step(residue *residues, int nres, vertx *vert, int nvert, float A, float proj_len,
                                float diel_int, float diel_ext, float sal, float ion_exc_rad, int phiType, int *i,
                                int step_size) {
  int it, eye, bound;
  int natoms;
  eye = *i;
  bound = eye + step_size;
  if (nvert < bound) bound = nvert;
  natoms = 0;
  for (it = 0; it < nres; it++) natoms += residues[it].natoms;

  float *res_c = (float *)malloc(sizeof(float) * nres), *res_x = (float *)malloc(sizeof(float) * nres),
        *res_y = (float *)malloc(sizeof(float) * nres), *res_z = (float *)malloc(sizeof(float) * nres),
        *at_c = (float *)malloc(sizeof(float) * natoms), *at_x = (float *)malloc(sizeof(float) * natoms),
        *at_y = (float *)malloc(sizeof(float) * natoms), *at_z = (float *)malloc(sizeof(float) * natoms),
        *vert_c = (float *)malloc(sizeof(float) * nvert), *vert_x = (float *)malloc(sizeof(float) * nvert),
        *vert_y = (float *)malloc(sizeof(float) * nvert), *vert_z = (float *)malloc(sizeof(float) * nvert),
        *vert_x_p = (float *)malloc(sizeof(float) * nvert), *vert_y_p = (float *)malloc(sizeof(float) * nvert),
        *vert_z_p = (float *)malloc(sizeof(float) * nvert);
  unsigned int *atom_addrs = (unsigned int *)malloc(nres * sizeof(unsigned int)),
               *atom_lengths = (unsigned int *)malloc(nres * sizeof(unsigned int));

  /* flatten the atom structure (copy to the "device" arrays) */
  natoms = 0;
  for (it = 0; it < nres; it++) {
    atom_addrs[it] = natoms;
    res_c[it] = residues[it].total_charge;
    res_x[it] = residues[it].x;
    res_y[it] = residues[it].y;
    res_z[it] = residues[it].z;
    for (int j = 0; j < residues[it].natoms; j++) {
      at_c[natoms + j] = residues[it].atoms[j].charge;
      at_x[natoms + j] = residues[it].atoms[j].x;
      at_y[natoms + j] = residues[it].atoms[j].y;
      at_z[natoms + j] = residues[it].atoms[j].z;
    }
    natoms += residues[it].natoms;
    atom_lengths[it] = residues[it].natoms;
  }
  for (it = 0; it < nvert; it++) {
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
  float r, r0, rprime;
  r = A + proj_len;
  rprime = A + ion_exc_rad;
  r0 = A / 2.; /* the kernel recomputes r0 = A/2 itself */

  /* the residue arrays are uploaded in the OpenCL code although the kernel never reads them: keep that traffic */
  dwarf_init("gem", "n-body-methods", "omp");
  dwarf_set_size(nvert, natoms);
  int eye0 = eye, step0 = step_size;
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    eye = eye0; step_size = step0;
    memset(vert_c, 0, sizeof(float) * nvert);
    dwarf_region_begin();
#pragma omp target enter data map(to : res_c[0 : nres], res_x[0 : nres], res_y[0 : nres], res_z[0 : nres], at_c[0 : natoms], at_x[0 : natoms], at_y[0 : natoms], at_z[0 : natoms], vert_c[0 : nvert], vert_x[0 : nvert], vert_y[0 : nvert], vert_z[0 : nvert], vert_x_p[0 : nvert], vert_y_p[0 : nvert], vert_z_p[0 : nvert], atom_addrs[0 : nres], atom_lengths[0 : nres])
    for (; eye < bound; eye += BLOCK_DIM_X_C * BLOCK_DIM_Y_C * BLOCKS) {
      if (rep == dwarf_first_rep()) fprintf(stdout, "finished first %d of %d\n", eye, nvert);
      run_kernel_block(nres, nvert, A, proj_len, diel_int, diel_ext, sal, ion_exc_rad, phiType, eye, step_size, atom_addrs,
                       atom_lengths, r, r0, rprime, res_c, res_x, res_y, res_z, at_c, at_x, at_y, at_z, vert_c, vert_x,
                       vert_y, vert_z, vert_x_p, vert_y_p, vert_z_p, natoms);
      if (eye + step_size > bound) step_size = bound - eye;
    }
    if (eye > bound) eye = bound;
#pragma omp target exit data map(from : vert_c[0 : nvert]) map(delete : res_c[0 : nres], res_x[0 : nres], res_y[0 : nres], res_z[0 : nres], at_c[0 : natoms], at_x[0 : natoms], at_y[0 : natoms], at_z[0 : natoms], vert_x[0 : nvert], vert_y[0 : nvert], vert_z[0 : nvert], vert_x_p[0 : nvert], vert_y_p[0 : nvert], vert_z_p[0 : nvert], atom_addrs[0 : nres], atom_lengths[0 : nres])
    dwarf_region_end(rep);
  }
  gem_valid = 1;
  for (it = 0; it < nvert; it++) {
    vert[it].potential = vert_c[it];
    if (!std::isfinite(vert_c[it])) gem_valid = 0;
  }
  *i = eye;
  free(res_c); free(res_x); free(res_y); free(res_z); free(at_c); free(at_x); free(at_y); free(at_z);
  free(vert_c); free(vert_x); free(vert_y); free(vert_z); free(vert_x_p); free(vert_y_p); free(vert_z_p);
  free(atom_addrs); free(atom_lengths);
}
