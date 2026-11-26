#include <stdio.h>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "dwarf_cuda.h"
#include "dwarf_timer.h"
#include "dwarf_launch.h"
#include <cmath>
/* host copies of the mesh, kept so every timed repetition can re-upload it */
static float *g_areas = 0, *g_normals = 0; static int *g_ese = 0;
#include <malloc.h>

#define AOCL_ALIGNMENT 64

/* float3 is CUDA's built-in 12-byte struct (same layout as the OpenCL host's own typedef) */

/*
 * Options
 *
 */
#define GAMMA 1.4
#define iterations 1

#define NDIM 3
#define NNB 4

#define RK 3	// 3rd order RK
#define ff_mach 1.2
#define deg_angle_of_attack 0.0f

#define VAR_DENSITY 0
#define VAR_MOMENTUM  1
#define VAR_DENSITY_ENERGY (VAR_MOMENTUM+NDIM)
#define NVAR (VAR_DENSITY_ENERGY+1)

#ifndef block_length
#define block_length 128
#endif

#include "cfd_kernel.cuh"


/*
 * Generic functions
 */
	template <typename T>
T* alloc(int N)
{
	T* mem = NULL;
	CUDA_CHECK(cudaMalloc((void**)&mem, sizeof(T)*N));
	return mem;
}

	template <typename T>
void dealloc(T* array)
{
	CUDA_CHECK(cudaFree(array));
}


	template <typename T>
void copy(T* dst, T* src, int N)
{
	CUDA_CHECK(cudaMemcpy(dst, src, N*sizeof(T), cudaMemcpyDeviceToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
}

	template <typename T>
void upload(T* dst, T* src, int N)
{
	CUDA_CHECK(cudaMemcpy(dst, src, sizeof(T) * N, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
}

	template <typename T>
void download(T* dst, T* src, int N)
{
	CUDA_CHECK(cudaMemcpy(dst, src, sizeof(T)*N, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());
}

void dump(float* variables, int nel, int nelr)
{
	float* h_variables;
    h_variables  = (float*)  memalign(AOCL_ALIGNMENT,nelr*NVAR*sizeof(float));
	download(h_variables, variables, nelr*NVAR);

	{
		std::ofstream file("density.txt");
		file << nel << " " << nelr << std::endl;
		for(int i = 0; i < nel; i++) file << h_variables[i + VAR_DENSITY*nelr] << std::endl;
	}


	{
		std::ofstream file("momentum.txt");
		file << nel << " " << nelr << std::endl;
		for(int i = 0; i < nel; i++)
		{
			for(int j = 0; j != NDIM; j++)
				file << h_variables[i + (VAR_MOMENTUM+j)*nelr] << " ";
			file << std::endl;
		}
	}

	{
		std::ofstream file("density_energy.txt");
		file << nel << " " << nelr << std::endl;
		for(int i = 0; i < nel; i++) file << h_variables[i + VAR_DENSITY_ENERGY*nelr] << std::endl;
	}
	free(h_variables);
}

inline void compute_flux_contribution(float* density, float3* momentum, float* density_energy,
		float pressure, float3* velocity,
		float3* fc_momentum_x, float3* fc_momentum_y, float3* fc_momentum_z,
		float3* fc_density_energy)
{
	fc_momentum_x->x = velocity->x*momentum->x + pressure;
	fc_momentum_x->y = velocity->x*momentum->y;
	fc_momentum_x->z = velocity->x*momentum->z;


	fc_momentum_y->x = fc_momentum_x->y;
	fc_momentum_y->y = velocity->y*momentum->y + pressure;
	fc_momentum_y->z = velocity->y*momentum->z;

	fc_momentum_z->x = fc_momentum_x->z;
	fc_momentum_z->y = fc_momentum_y->z;
	fc_momentum_z->z = velocity->z*momentum->z + pressure;

	float de_p = *density_energy+pressure;
	fc_density_energy->x = velocity->x*de_p;
	fc_density_energy->y = velocity->y*de_p;
	fc_density_energy->z = velocity->z*de_p;
}

int main(int argc, char** argv)
{
	dwarf_cuda_init(&argc, &argv);

	size_t global_size;

	float *ff_variable;
	float3 *ff_fc_momentum_x, *ff_fc_momentum_y, *ff_fc_momentum_z, *ff_fc_density_energy;

	if (argc < 2)
	{
		printf("Usage ./cfd <data input file>\n");
		return 0;
	}


	const char* data_file_name = argv[1];


	// set far field conditions and load them into constant memory on the gpu
	{
		float h_ff_variable[NVAR];
		const float angle_of_attack = (float)(3.1415926535897931 / 180.0) * (float)(deg_angle_of_attack);

		h_ff_variable[VAR_DENSITY] = (float)(1.4);

		float ff_pressure = (float)(1.0);
		float ff_speed_of_sound = sqrt(GAMMA*ff_pressure / h_ff_variable[VAR_DENSITY]);
		float ff_speed = (float)(ff_mach)*ff_speed_of_sound;

		float3 ff_velocity;
		ff_velocity.x = ff_speed*(float)(cos((float)angle_of_attack));
		ff_velocity.y = ff_speed*(float)(sin((float)angle_of_attack));
		ff_velocity.z = 0.0;

		h_ff_variable[VAR_MOMENTUM+0] = h_ff_variable[VAR_DENSITY] * ff_velocity.x;
		h_ff_variable[VAR_MOMENTUM+1] = h_ff_variable[VAR_DENSITY] * ff_velocity.y;
		h_ff_variable[VAR_MOMENTUM+2] = h_ff_variable[VAR_DENSITY] * ff_velocity.z;

		h_ff_variable[VAR_DENSITY_ENERGY] = h_ff_variable[VAR_DENSITY]*((float)(0.5)*(ff_speed*ff_speed)) + (ff_pressure / (float)(GAMMA-1.0));

		float3 h_ff_momentum;
		h_ff_momentum.x = *(h_ff_variable+VAR_MOMENTUM+0);
		h_ff_momentum.y = *(h_ff_variable+VAR_MOMENTUM+1);
		h_ff_momentum.z = *(h_ff_variable+VAR_MOMENTUM+2);
		float3 h_ff_fc_momentum_x;
		float3 h_ff_fc_momentum_y;
		float3 h_ff_fc_momentum_z;
		float3 h_ff_fc_density_energy;
		compute_flux_contribution(&h_ff_variable[VAR_DENSITY], &h_ff_momentum,
				&h_ff_variable[VAR_DENSITY_ENERGY], ff_pressure, &ff_velocity,
				&h_ff_fc_momentum_x, &h_ff_fc_momentum_y, &h_ff_fc_momentum_z,
				&h_ff_fc_density_energy);

		// copy far field conditions to the gpu (CL_MEM_COPY_HOST_PTR)
		ff_variable = alloc<float>(NVAR); upload<float>(ff_variable, h_ff_variable, NVAR);
		ff_fc_momentum_x = alloc<float3>(1); upload<float3>(ff_fc_momentum_x, &h_ff_fc_momentum_x, 1);
		ff_fc_momentum_y = alloc<float3>(1); upload<float3>(ff_fc_momentum_y, &h_ff_fc_momentum_y, 1);
		ff_fc_momentum_z = alloc<float3>(1); upload<float3>(ff_fc_momentum_z, &h_ff_fc_momentum_z, 1);
		ff_fc_density_energy = alloc<float3>(1); upload<float3>(ff_fc_density_energy, &h_ff_fc_density_energy, 1);
	}
	int nel;
	int nelr;

	// read in domain geometry
	float *areas;
	int *elements_surrounding_elements;
	float *normals;
	{
		std::ifstream file(data_file_name);

		file >> nel;

		nelr = block_length*((nel / block_length )+ std::min(1, nel % block_length));

		//float* h_areas = new float[nelr];
		//int* h_elements_surrounding_elements = new int[nelr*NNB];
		//float* h_normals = new float[nelr*NDIM*NNB];

		float* h_areas ;
		int* h_elements_surrounding_elements ;
		float* h_normals ;

        h_areas                         = (float*)  memalign(AOCL_ALIGNMENT,nelr*sizeof(float));
        h_elements_surrounding_elements = (int*)    memalign(AOCL_ALIGNMENT,nelr*NNB*sizeof(int));
        h_normals                       = (float *) memalign(AOCL_ALIGNMENT,nelr*NDIM*NNB*sizeof(float));


        //posix_memalign(&h_areas                         , AOCL_ALIGNMENT, nelr);
        //posix_memalign(&h_elements_surrounding_elements , AOCL_ALIGNMENT, nelr*NNB);
        //posix_memalign(&h_normals                       , AOCL_ALIGNMENT, nelr*NDIM*NNB);


		// read in data
		for(int i = 0; i < nel; i++)
		{
			file >> h_areas[i];
			for(int j = 0; j < NNB; j++)
			{
				file >> h_elements_surrounding_elements[i + j*nelr];
				if(h_elements_surrounding_elements[i+j*nelr] < 0) h_elements_surrounding_elements[i+j*nelr] = -1;
				h_elements_surrounding_elements[i + j*nelr]--; //it's coming in with Fortran numbering

				for(int k = 0; k < NDIM; k++)
				{
					file >> h_normals[i + (j + k*NNB)*nelr];
					h_normals[i + (j + k*NNB)*nelr] = -h_normals[i + (j + k*NNB)*nelr];
				}
			}
		}

		// fill in remaining data
		int last = nel-1;
		for(int i = nel; i < nelr; i++)
		{
			h_areas[i] = h_areas[last];
			for(int j = 0; j < NNB; j++)
			{
				// duplicate the last element
				h_elements_surrounding_elements[i + j*nelr] = h_elements_surrounding_elements[last + j*nelr];
				/* [FIX-OCL-10] was h_normals[last + ...] = h_normals[last + ...] (a self-assignment): the padded elements nel..nelr-1 kept
				 * uninitialised normals although the comment says "duplicate the last element" */
				for(int k = 0; k < NDIM; k++) h_normals[i + (j + k*NNB)*nelr] = h_normals[last + (j + k*NNB)*nelr];
			}
		}

		areas = alloc<float>(nelr);
		upload<float>(areas, h_areas, nelr);

		elements_surrounding_elements = alloc<int>(nelr*NNB);
		upload<int>(elements_surrounding_elements, h_elements_surrounding_elements, nelr*NNB);

		normals = alloc<float>(nelr*NDIM*NNB);
		upload<float>(normals, h_normals, nelr*NDIM*NNB);

		g_areas = h_areas; g_ese = h_elements_surrounding_elements; g_normals = h_normals;  // kept for the timed re-uploads
	}

	/* one element per thread: DWARF_LS threads per block, ceil(nelr/DWARF_LS) blocks (kernels guard i < nelr) */
	global_size = DWARF_ROUNDUP(nelr);
	unsigned int nblocks = DWARF_NBLOCKS(nelr);
	(void)global_size;

	// Create arrays and set initial conditions
	float *variables = alloc<float>(nelr*NVAR);
	initialize_variables<<<nblocks, DWARF_LS>>>(nelr, variables, ff_variable);
	CUDA_SYNC_LAUNCH();

	float *old_variables = alloc<float>(nelr*NVAR);
	float *fluxes = alloc<float>(nelr*NVAR);
	float *step_factors = alloc<float>(nelr);
	float *fc_momentum_x = alloc<float>(nelr*NDIM);
	float *fc_momentum_y = alloc<float>(nelr*NDIM);
	float *fc_momentum_z = alloc<float>(nelr*NDIM);
	float *fc_density_energy = alloc<float>(nelr*NDIM);

	// make sure all memory is floatly allocated before we start timing
	initialize_variables<<<nblocks, DWARF_LS>>>(nelr, old_variables, ff_variable);
	CUDA_SYNC_LAUNCH();
	initialize_variables<<<nblocks, DWARF_LS>>>(nelr, fluxes, ff_variable);
	CUDA_SYNC_LAUNCH();
	std::cout << "About to memcopy" << std::endl;
	dealloc<float>(step_factors);
	float temp[nelr];
	for(int i = 0; i < nelr; i++)
		temp[i] = 0;
	step_factors = alloc<float>(nelr);
	upload<float>(step_factors, temp, nelr);
	// make sure CUDA isn't still doing something before we start timing

	CUDA_CHECK(cudaDeviceSynchronize());

	// these need to be computed the first time in order to compute time step
	std::cout << "Starting..." << std::endl;

	/* common timer: repetitions around [mesh upload, variable initialisation, the solver iterations, download of the
	 * solution]; JIT/build is above, writing the output files is below. */
	dwarf_init("cfd", "unstructured-grids", "cuda");
	dwarf_set_size(nel, nelr);
	float *h_result = (float*) malloc(sizeof(float) * nelr * NVAR);
	for(int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	dwarf_region_begin();
	upload<float>(areas, g_areas, nelr);
	upload<int>(elements_surrounding_elements, g_ese, nelr*NNB);
	upload<float>(normals, g_normals, nelr*NDIM*NNB);
	{
		float *targets[3] = { variables, old_variables, fluxes };
		for(int t = 0; t < 3; t++)
		{
			initialize_variables<<<nblocks, DWARF_LS>>>(nelr, targets[t], ff_variable);
			CUDA_SYNC_LAUNCH();
		}
		upload<float>(step_factors, temp, nelr);   // step_factors = 0
	}


	// Begin iterations
	for(int i = 0; i < iterations; i++)
	{
		copy<float>(old_variables, variables, nelr*NVAR);

		// for the first iteration we compute the time step
		compute_step_factor<<<nblocks, DWARF_LS>>>(nelr, variables, areas, step_factors);
		CUDA_SYNC_LAUNCH();
		for(int j = 0; j < RK; j++)
		{
			compute_flux_contributions<<<nblocks, DWARF_LS>>>(nelr, variables, fc_momentum_x, fc_momentum_y, fc_momentum_z, fc_density_energy);
			CUDA_SYNC_LAUNCH();
			compute_flux<<<nblocks, DWARF_LS>>>(nelr, elements_surrounding_elements, normals, variables, fc_momentum_x, fc_momentum_y,
					fc_momentum_z, fc_density_energy, fluxes, ff_variable, ff_fc_momentum_x, ff_fc_momentum_y, ff_fc_momentum_z, ff_fc_density_energy);
			CUDA_SYNC_LAUNCH();
			time_step<<<nblocks, DWARF_LS>>>(j, nelr, old_variables, variables, step_factors, fluxes);
			CUDA_SYNC_LAUNCH();
		}
	}

	CUDA_CHECK(cudaDeviceSynchronize());
	download<float>(h_result, variables, nelr*NVAR);
	dwarf_region_end(rep);
	}
	std::cout << "Finished" << std::endl;
	std::cout << "Saving solution..." << std::endl;
	dump(variables, nel, nelr);
	std::cout << "Saved solution..." << std::endl;
	std::cout << "Cleaning up..." << std::endl;

	dealloc<float>(areas);
	dealloc<int>(elements_surrounding_elements);
	dealloc<float>(normals);

	dealloc<float>(variables);
	dealloc<float>(old_variables);
	dealloc<float>(fluxes);
	dealloc<float>(step_factors);
	dealloc<float>(fc_momentum_x);
	dealloc<float>(fc_momentum_y);
	dealloc<float>(fc_momentum_z);
	dealloc<float>(fc_density_energy);

	std::cout << "Done..." << std::endl;
	/* [FIX-OCL-11] validity check (none before): the solution must be finite with positive density */
	{
		long bad = 0;
		for(int i = 0; i < nel; i++) { float d = h_result[i + VAR_DENSITY*nelr]; if(!(std::isfinite(d) && d > 0.0f)) bad++; }
		printf("%s: %ld cells with non-finite or non-positive density\n", bad ? "FAILURE" : "SUCCESS", bad);
		return dwarf_report(bad == 0);
	}
}
