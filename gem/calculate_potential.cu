__device__ inline float dist2D(float x, float y, float z, float x2, float y2, float z2)
{
     /* local variables */
      float xd = x - x2,
             yd = y - y2,
             zd = z - z2;
      //return (__fadd_rn(__fmul_rn(xd,xd) , __fadd_rn( __fmul_rn(yd,yd) ,  __fmul_rn(zd,zd))));
      return (xd*xd + yd*yd + zd*zd);

} /* end function dist2 */

__device__ inline float distD (float x, float y, float z, float x2, float y2, float z2)
{
      return sqrt(dist2D(x, y, z, x2, y2, z2));
} /* end function dist */

__device__ float calculate_potential_contribution (float d, float dprime, float charge, analytical_definitions_struct defs)
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
          
   float temp1;
          
         
   /* <><> add salt to the inside <><> */

   if (defs.phiType != TOTAL_POTENTIAL) 
   {
        coulomb = charge/(d * diel_int_dev);
	//coulomb = __fdividef(charge, __fmul_rn(d, diel_int_dev));
   }

   if (defs.phiType & REACTION_POTENTIAL)
   {
      if (defs.region == 3)  /* Dist to surf > ion_excl_rad */
      {
   	/*
         temp1 = __fdividef(__fmul_rn(charge, (exp (-defs.kappa * (d - dprime)))) , d);
         sum1 = __fmul_rn(temp1, __fdividef(defs.one_plus_alpha , (1. + __fmul_rn(defs.kappa, dprime))));
         sum2 = __fmul_rn(charge, (defs.alpha_by_one_minus_beta * exp (__fmul_rn(-defs.kappa , (r_dev-rprime_dev))) / __fmul_rn(r_dev , (1. + defs.kappa*rprime_dev))));
         */
         sum1 = (charge 
                    * (exp (-defs.kappa * (d - dprime))) / d)
                       * (defs.one_plus_alpha / (1. + defs.kappa * dprime));
   
         sum2 = charge 
                  * (defs.alpha_by_one_minus_beta * exp (-defs.kappa * (r_dev-rprime_dev)) / (r_dev * (1. + defs.kappa*rprime_dev)));
	 
	 
         to_return = defs.inverse_one_plus_ab_by_diel_ext * (sum1 - sum2);
      }
      else if (defs.region == 2) /* 0 < dist_to_surf <= ion_excl_rad */
      {
   
         one_over_one_plus_kappa_rprime = 1./(1. + (defs.kappa * rprime_dev));
	 //one_over_one_plus_kappa_rprime = 1./__fadd_rn(1. , __fmul_rn(defs.kappa , rprime_dev));
   
         /* electrostatic terms for region 2 */
         /*
         sum1 = __fdividef(__fmul_rn(charge, defs.one_plus_alpha), d);
         sum2 = __fdividef(__fmul_rn(defs.alpha_by_one_minus_beta, charge), r_dev);
         */
         sum1 = defs.one_plus_alpha * charge / d;
         sum2 = defs.alpha_by_one_minus_beta * charge/r_dev;
   	
   	 
         /* salt terms for region 2 */
         /*
         sum3 = __fmul_rn(charge, __fmul_rn(__fdividef(defs.alpha_by_one_minus_beta , rprime_dev) , (1. - one_over_one_plus_kappa_rprime)));
         temp1 = __fmul_rn(defs.one_plus_alpha, charge);
         sum4 = __fmul_rn(temp1, ( 1.0 / (dprime + __fmul_rn(__fmul_rn(defs.kappa, dprime), dprime) - (1. / dprime))));
         */
         sum3 = charge * (defs.alpha_by_one_minus_beta / rprime_dev) * (1. - one_over_one_plus_kappa_rprime);
         sum4 = defs.one_plus_alpha * charge
                  * ( 1.0 / (dprime + defs.kappa * dprime * dprime) - (1. / dprime));
       	 
   	
         to_return = defs.inverse_one_plus_ab_by_diel_ext * (sum1 - sum2 + sum3 + sum4);
      }
      else /* region = 1 */
      {

         one_over_one_plus_kappa_rprime = 1./(1. + (defs.kappa * rprime_dev));
   
         sum1 =  1. / (d * diel_int_dev);
   
         sum2 =  (1./diel_int_dev - 1./diel_ext_dev) / (defs.one_plus_alpha_beta * A_dev);
   
         sum3 = defs.Asq / sqrt(defs.Asq_minus_rnautsq*defs.Asq_minus_rsq + defs.Asq_by_dsq);
	/*
	 salt = __fmul_rn(defs.inverse_one_plus_ab_by_diel_ext,
                   (__fmul_rn(__fdividef(defs.one_plus_alpha, dprime), (1./(1.+defs.kappa*dprime)-1.))  
                     + __fmul_rn(__fdividef(defs.alpha_by_one_minus_beta, rprime_dev), (1.-one_over_one_plus_kappa_rprime))));
        */
                     
         salt = defs.inverse_one_plus_ab_by_diel_ext
                  * (defs.one_plus_alpha/dprime*(1./(1.+defs.kappa*dprime)-1.)  
                     +defs.alpha_by_one_minus_beta/rprime_dev*(1.-one_over_one_plus_kappa_rprime));
   	 
   		
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
   
    //if (to_return != to_return)
   //{
   //  printf("Found a NaN contribution %f\n", to_return);
   //  printf("region = %i\n", defs.region);
   //  printf("sum1, sum2, sum3, sum4 %f %f %f %f\n", sum1, sum2, sum3, sum4);

   //  if (defs.region == 1)
   //  {
   //     printf("sum under root = %f\n", defs.Asq_minus_rnautsq*defs.Asq_minus_rsq + defs.Asq_by_dsq);
   //     printf("A-r by A-r0 = %f\n", defs.Asq_minus_rnautsq*defs.Asq_minus_rsq);
   //     printf("Asq_by_dsq = %f\n", defs.Asq_by_dsq);
   //     printf("A = %f\n", A);
   //     printf("r = %f\n", r);
   //     printf("r0 = %f\n", r0);
   //     printf("d = %f\n", d);
   //  }

   //  fflush(stdout);

   //}

   return to_return;
}

/* Calculates potential for surface out of a molecule -- this is the free
 * space solution. --JCG */
__global__ void calc_potential_single_step_dev(int eye, int step_size, partitioned_open_struct_dev *devStructPtr)
{
	
	int i, j=0, k, m, natoms, i2=0, i3=0, ni2=0, nk=0, use_charge_approx;
	float d=0.0, d2=0.0, dprime=0.0;

	//partitioned_open_struct_dev *devStructPtr;
	
	analytical_definitions_struct defs;
	
	extern __shared__ float4 atom_vert[];
	
	//devStructPtr=(partitioned_open_struct_dev *)devStruct;
	
	//i=blockIdx.x*blockDim.x;
	eye = eye + __mul24(blockIdx.x, __mul24(blockDim.x,  blockDim.y)) + __mul24(blockDim.x, threadIdx.y) + threadIdx.x;
	if(eye >= nvert_dev)//TODO:double check this
    		return;

	atom_vert[threadIdx.x].x=devStructPtr->xyz_vert[eye/*blockIdx.x*blockDim.x+threadIdx.x*/].x;
	atom_vert[threadIdx.x].y=devStructPtr->xyz_vert[eye/*blockIdx.x*blockDim.x+threadIdx.x*/].y;
	atom_vert[threadIdx.x].z=devStructPtr->xyz_vert[eye/*blockIdx.x*blockDim.x+threadIdx.x*/].z;
	atom_vert[threadIdx.x].w=0.0;
		//getchar();
	
	__syncthreads();
	
	/*for(j=0;j<64;j++)
	{
		printf("atom %f xyz %f", atom_vert[j].x, devStructPtr->xyz_vert[j].x);
		getchar();
	}
	__syncthreads();
	*/
	
	
	defs.Asq   = A_dev*A_dev;
        defs.Asq_minus_rsq = defs.Asq - (r_dev*r_dev);

   	defs.kappa = 0.316 * sqrt(sal_dev);
   	defs.beta = diel_int_dev/diel_ext_dev;
   	defs.alpha_by_one_minus_beta = ALPHA_OF*(1. - defs.beta);
   	defs.alpha_beta = ALPHA_OF * defs.beta;
   	defs.one_plus_alpha_beta = 1.0 + defs.alpha_beta;
   	defs.one_plus_alpha = 1.0 + ALPHA_OF;
   	defs.inverse_one_plus_ab_by_diel_ext = 1.0/(defs.one_plus_alpha_beta*diel_ext_dev);
   	defs.phiType = phiType_dev;
   
   	if (r_dev > rprime_dev)
   	{
      		defs.region = 3;
   	}
   	else if (r_dev > A_dev)
   	{
      		defs.region = 2;
   	}
   	else
  	{
      		defs.region = 1;
   	}
   	
   	
   	defs.Asq_minus_rnautsq = defs.Asq - r0_dev*r0_dev;
   	
             for (i2 = 0; i2<devStructPtr->n_lvl2; i2++)
             {
                /* check if charge approx should be used */
                use_charge_approx = NO;
                if (CHARGE_APPROX)
                {
                   d2 = dist2D
                        (
                              atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                              devStructPtr->xyz_lvl2[i2].x,
                              devStructPtr->xyz_lvl2[i2].y,
                              devStructPtr->xyz_lvl2[i2].z
                        );
                        
                   use_charge_approx = (d2 > CHARGE_APPROX_L2_DIST2);
                }

                /* charge approximation */
                if (use_charge_approx)
                {
 	           //devStructPtr->com_cnt[eye]++;
                   for (m = 0; m < CHARGE_APPROX; m++)
                   {
                      if (fabs(devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].w) > 0.01) 
                      {                                        
                         d = distD
                            (
                              atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                              devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].x,
                              devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].y,
                              devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].z
                            );
                         /* distance from point to charge */
                         dprime = distD
                            (
                              atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                              devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].x,
                              devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].y,
                              devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].z
                            );   
                         /* need to estimate r0, dist from atom to surface */
                         /* if we intend to support negative projection lengths <><> */
                         defs.Asq_by_dsq = A_dev*A_dev * d*d;
                         devStructPtr->vert_c[eye] +=
                            calculate_potential_contribution(d, dprime, devStructPtr->xyzq_lvl2[CHARGE_APPROX*i2+m].w, defs);
                      }  /* end if charged */
                   }  /* end for m CHARGE_APPROX loop*/
                }  /* end if level 2 charge approx */
                else  /* inside level 2 distance, check level 1 */
                {
                   if (i2 < devStructPtr->n_lvl2 - 1) { nk = devStructPtr->i_lvl2[i2+1]; }
                   else { nk = devStructPtr->nresidues; } 
			#pragma unroll 
                   for (k = devStructPtr->i_lvl2[i2]; k < nk; k++)
    // #endif
    //		for(k=0;k<devStructPtr->nresidues;k++)
                   {
                      /* check if charge approx should be used */
                      use_charge_approx = NO;
                      if (CHARGE_APPROX)
                      {
                         d2 = dist2D
                              (
                                    atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                    devStructPtr->residues[k].x,
                                    devStructPtr->residues[k].y,
                                    devStructPtr->residues[k].z
                              );
                         use_charge_approx = (d2 > CHARGE_APPROX_DIST2);
                      }
                      /* charge approximation */
                      if (use_charge_approx)
                      {
 	              // devStructPtr->com_cnt[eye]++;
 	              	if(CHARGE_APPROX == 1)
 	              	{
 	              		if ((devStructPtr->residues[k].qc[0] > .01) 
                        	       || (devStructPtr->residues[k].qc[0] < -.01))   /* charged */
                        	{
                        		d = distD
                                	  (
                                	    atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                	    devStructPtr->residues[k].xc[0],
                                	    devStructPtr->residues[k].yc[0],
                                	    devStructPtr->residues[k].zc[0]
                                	  );
                               		/* distance from point to charge */
                               		dprime = distD
                                  	(
	                                     atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
        	                             devStructPtr->residues[k].xc[0],
        	                             devStructPtr->residues[k].yc[0],
        	                             devStructPtr->residues[k].zc[0]
        	                          );   
                                  
                                 
        	                       defs.Asq_by_dsq = A_dev*A_dev * d*d;
        	                    devStructPtr->vert_c[eye] += calculate_potential_contribution(d, dprime, devStructPtr->residues[k].qc[0], defs);
                                       // printf("pot is %f eye is %d", devStructPtr->vert_c[eye], eye);
                                       // getchar();
                            }  /* end if charged */
 	              }
 	              if(CHARGE_APPROX == 2)
 	              {
 	              	if ((devStructPtr->residues[k].qc[0] > .01) 
                               || (devStructPtr->residues[k].qc[0] < -.01))   /* charged */
                        {
                        	d = distD
                                  (
                                    atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                    devStructPtr->residues[k].xc[0],
                                    devStructPtr->residues[k].yc[0],
                                    devStructPtr->residues[k].zc[0]
                                  );
                               /* distance from point to charge */
                               dprime = distD
                                  (
                                     atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                     devStructPtr->residues[k].xc[0],
                                     devStructPtr->residues[k].yc[0],
                                     devStructPtr->residues[k].zc[0]
                                  );   
                                  
                                 
                               defs.Asq_by_dsq = A_dev*A_dev * d*d;
                          	devStructPtr->vert_c[eye] += calculate_potential_contribution(d, dprime, devStructPtr->residues[k].qc[0], defs);
                                       // printf("pot is %f eye is %d", devStructPtr->vert_c[eye], eye);
                                       // getchar();
                        }  /* end if charged */
                        if ((devStructPtr->residues[k].qc[1] > .01) 
                               || (devStructPtr->residues[k].qc[1] < -.01))   /* charged */
                        {
                        	d = distD
                                  (
                                    atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                    devStructPtr->residues[k].xc[1],
                                    devStructPtr->residues[k].yc[1],
                                    devStructPtr->residues[k].zc[1]
                                  );
                               /* distance from point to charge */
                               dprime = distD
                                  (
                                     atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                     devStructPtr->residues[k].xc[1],
                                     devStructPtr->residues[k].yc[1],
                                     devStructPtr->residues[k].zc[1]
                                  );   
                                  
                                 
                               defs.Asq_by_dsq = A_dev*A_dev * d*d;
                            devStructPtr->vert_c[eye] += calculate_potential_contribution(d, dprime, devStructPtr->residues[k].qc[1], defs);
                                       // printf("pot is %f eye is %d", devStructPtr->vert_c[eye], eye);
                                       // getchar();
                            }  /* end if charged */
 	              }
                      	/*for (m = 0; m < CHARGE_APPROX; m++)
                         {
                            if ((devStructPtr->residues[k].qc[m] > .01) 
                               || (devStructPtr->residues[k].qc[m] < -.01))   /* charged */
                        /*    {
                               d = distD
                                  (
                                    atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                    devStructPtr->residues[k].xc[m],
                                    devStructPtr->residues[k].yc[m],
                                    devStructPtr->residues[k].zc[m]
                                  );
                               /* distance from point to charge */
                        /*       dprime = distD
                                  (
                                     atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                     devStructPtr->residues[k].xc[m],
                                     devStructPtr->residues[k].yc[m],
                                     devStructPtr->residues[k].zc[m]
                                  );   
                                  
                                 
                               defs.Asq_by_dsq = A_dev*A_dev * d*d;
                            devStructPtr->vert_c[eye] +=
                                  calculate_potential_contribution(d, dprime, devStructPtr->residues[k].qc[m], defs);
                                       // printf("pot is %f eye is %d", devStructPtr->vert_c[eye], eye);
                                       // getchar();
                            }  /* end if charged */
                       //  }  /* end for m CHARGE_APPROX loop*/
                      }  /* end if use charge approx */
                      else  /* do not use charge approx */
                      {
                         natoms = devStructPtr->residues[k].natoms;
			#pragma unroll 
                         for (j = 0; j < natoms; j++)
                         {
 	                  //   devStructPtr->com_cnt[eye]++;
                             /* distance from point to charge */
                             d = distD
                                   (
                                       atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                       devStructPtr->atoms_d[devStructPtr->atom_addr_d[k]+j].x,
                                       devStructPtr->atoms_d[devStructPtr->atom_addr_d[k]+j].y,
                                       devStructPtr->atoms_d[devStructPtr->atom_addr_d[k]+j].z
                                   );  
                             /* distance from point to charge */
                             dprime = distD
                                   (
                                       atom_vert[threadIdx.x].x, atom_vert[threadIdx.x].y, atom_vert[threadIdx.x].z,
                                       devStructPtr->atoms_d[devStructPtr->atom_addr_d[k]+j].x,
                                       devStructPtr->atoms_d[devStructPtr->atom_addr_d[k]+j].y,
                                       devStructPtr->atoms_d[devStructPtr->atom_addr_d[k]+j].z
                                   );  
                             
                             defs.Asq_by_dsq = A_dev*A_dev * d*d;
                              devStructPtr->vert_c[eye] +=
                                   calculate_potential_contribution
                                   	(d, dprime, devStructPtr->atoms_d[devStructPtr->atom_addr_d[k]+j].charge, defs);
   
                         }  /* end for j atoms */
                      }  /* end else do not use residue charge approx */
                   }  /* end for k residues */
                }  /* end else do not use level 2 charge approx */
             }  /* end for i2 level 2 components */
    //      }  /* end else do not use level 3 charge approx */
   //    }  /* end for i3 level 3 components */
}
}

