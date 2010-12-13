/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */



//#pragma once

//#include <stdio.h>

//#include "sparse_formats.h"
//#include "utils.h"
//#include "texture.h"

////////////////////////////////////////////////////////////////////////
// CSR SpMV kernels based on a scalar model (one thread per row)
///////////////////////////////////////////////////////////////////////
//
// spmv_csr_scalar_device
//   Straightforward translation of standard CSR SpMV to CUDA
//   where each thread computes y[i] += A[i,:] * x 
//   (the dot product of the i-th row of A with the x vector)
//
// spmv_csr_scalar_tex_device
//   Same as spmv_csr_scalar_device, except x is accessed via texture cache.
//

//template <typename unsigned int, typename float, bool UseCache>
__kernel void
spmv_csr_scalar_kernel(const unsigned int num_rows,
                       __global unsigned int * Ap, 
                       __global unsigned int * Aj, 
                       __global float * Ax, 
                       __global float * x, 
                       __global float * y)
{
    //Lazy, fix later
    //UseCache = false

    // row index

//((blockDim.x * (blockIdx.x + blockIdx.y*gridDim.x) + threadIdx.x))

//(get_local_size(0) * (get_group_id(0)+get_group_id(1)*get_num_groups(0))+get_local_id(0))

    __const unsigned int row = (get_local_size(0)*(get_group_id(0)+get_group_id(1)*get_num_groups(0))+get_local_id(0));
    
    if(row < num_rows){     
        float sum = y[row];

        __const unsigned int row_start = Ap[row];
        __const unsigned int row_end   = Ap[row+1];
    
        for (unsigned int jj = row_start; jj < row_end; jj++){             
            //sum += Ax[jj] * fetch_x<false>(Aj[jj], x);
	    //needs fixing, but then what dosent?
            sum += Ax[jj] * Aj[jj];       
        }

        y[row] = sum;
    }
}

//template <typename unsigned int, typename float>
//__kernel void 
//spmv_csr_scalar_device(const csr_matrix<unsigned int,float>& d_csr, 
//                            const float * d_x, 
//                                  float * d_y)
//{
//    const unsigned int BLOCK_SIZE = 256; 
    
//    const dim3 grid = make_large_grid(d_csr.num_rows, BLOCK_SIZE);

//    spmv_csr_scalar_kernel<unsigned int, float, false> <<<grid, BLOCK_SIZE>>> 
//        (d_csr.num_rows, d_csr.Ap, d_csr.Aj, d_csr.Ax, d_x, d_y);   
//}

//template <typename unsigned int, typename float>
//void spmv_csr_scalar_tex_device(const csr_matrix<unsigned int,float>& d_csr, 
//                                const float * d_x, 
//                                      float * d_y)
//{
//    const unsigned int BLOCK_SIZE = 256;
    
//    const dim3 grid = make_large_grid(d_csr.num_rows, BLOCK_SIZE);
    
//    bind_x(d_x);

//    spmv_csr_scalar_kernel<unsigned int, float, true> <<<grid, BLOCK_SIZE>>> 
//        (d_csr.num_rows, d_csr.Ap, d_csr.Aj, d_csr.Ax, d_x, d_y);   

//    unbind_x(d_x);
//}

