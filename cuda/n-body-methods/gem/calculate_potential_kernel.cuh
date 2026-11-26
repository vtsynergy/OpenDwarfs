/* CUDA translation of calculate_potential.cl. OpenCL native_exp -> GEM_NATIVE_EXP (see below). */
#ifndef GEM_NATIVE_EXP
#define GEM_NATIVE_EXP __expf
#endif

#define GEM_ALPHA_OF .571412f /* "Incorporating variable dielectric environments into the generalized Born model" Sigalov, Scheffe, Onufriev J. Chem. Phys. */
#define GEM_REACTION_POTENTIAL 2
#define GEM_TOTAL_POTENTIAL 3
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

} gem_dev_defs_struct;



__device__ static inline float dist2D(float x, float y, float z, float x2, float y2, float z2)
{
	/* local variables */
	float xd = x - x2,
	      yd = y - y2,
	      zd = z - z2;
	return (xd*xd + yd*yd + zd*zd);
} /* end function dist2 */

__device__ static inline float distD (float x, float y, float z, float x2, float y2, float z2)
{
	return sqrtf(dist2D(x, y, z, x2, y2, z2));
} /* end function dist */

/* Calculates potential for surface out of a molecule -- this is the free
 * space solution. --JCG */
__global__ void calc_potential_single_step_dev(
		int nres,
		int nvert,
		float A,
		float proj_len,
		float diel_int,
		float diel_ext,
		float sal,
		float ion_exc_rad,
		int phiType,
		int eye,
		int step_size,
		unsigned int * atom_addrs,
		unsigned int * atom_lengths,
		float r,
		float r0,
		float rprime,
		float *res_c_s, float *res_x_s, float *res_y_s, float *res_z_s,
		float *at_c_s, float *at_x_s, float *at_y_s, float *at_z_s,
		float *vert_c_s, float *vert_x_s, float *vert_y_s, float *vert_z_s,
		float *vert_x_p_s, float *vert_y_p_s, float *vert_z_p_s)
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



	gem_dev_defs_struct defs;
	//calculate bound
	//bound = eye+step_size < nvert ? nvert : eye + step_size;
	//calculate the vertx on which to act
	/* OpenCL: global {256,BLOCKS}, local {256,1}; eye + global_size(0)*global_id(1) + global_id(0) (global_id(2) = 0).
	 * CUDA: grid (1, BLOCKS), block 256: global_id(0) = threadIdx.x, global_id(1) = blockIdx.y, global_size(0) = 256. */
	eye = eye + (int)((gridDim.x*blockDim.x)*(blockIdx.y*blockDim.y + threadIdx.y)) + (int)(blockIdx.x*blockDim.x + threadIdx.x);
	if(eye >= nvert)//TODO:double check this
		return;

	//this shouldn't be here, but the kernel wigs out without it, think there might be a maximum struct size for arguments... must look into this
	defs.Asq   = A*A;
	defs.Asq_minus_rsq = defs.Asq - (r*r);

	defs.kappa = 0.316f * sqrtf(sal);
	defs.beta = diel_int/diel_ext;
	defs.alpha_by_one_minus_beta = GEM_ALPHA_OF*(1.0f - defs.beta);
	defs.alpha_beta = GEM_ALPHA_OF * defs.beta;
	defs.one_plus_alpha_beta = 1.0f + defs.alpha_beta;
	defs.one_plus_alpha = 1.0f + GEM_ALPHA_OF;
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
	r0 = A/2.0f;  /* might be tricky ... hrmmm */
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
				if (defs.phiType != GEM_TOTAL_POTENTIAL)
				{
					coulomb = charge/(d * diel_int);
				}


				if (defs.phiType & GEM_REACTION_POTENTIAL)
				{
					if (defs.region == 3)  /* Dist to surf > ion_excl_rad */
					{

						sum1 = (charge
								* (GEM_NATIVE_EXP(-defs.kappa * (d - dprime))) / d)
							* (defs.one_plus_alpha / (1.0f + defs.kappa * dprime));

						sum2 = charge
							* (defs.alpha_by_one_minus_beta * GEM_NATIVE_EXP(-defs.kappa * (r-rprime)) / (r * (1.0f + defs.kappa*rprime)));

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

					if (defs.phiType == GEM_REACTION_POTENTIAL)
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
