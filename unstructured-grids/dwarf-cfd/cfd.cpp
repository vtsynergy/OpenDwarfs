#include <stdio.h>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iostream>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#define CHKERR(err, str) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "CL Error %d: %s\n", err, str); \
        exit(1); \
    }

#define USEGPU 0

typedef struct{
    float x;
    float y;
    float z;
} float3;

/*
 * Options
 *
 */
#define GAMMA 1.4
#define iterations 2000

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

const char *KernelSourceFile = "cfd_kernel.cl";

/*
 * Generic functions
 */
template <typename T>
cl_mem alloc(cl_context context, int N)
{
    int err;
	cl_mem mem = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(T)*N, NULL, &err);
	CHKERR(err, "Unable to allocate memory");
	return mem;
}

template <typename T>
void dealloc(cl_mem array)
{
    int err = clReleaseMemObject(array);
	CHKERR(err, "Unable to release memory");
}


template <typename T>
void copy(cl_command_queue commands, cl_mem dst, cl_mem src, int N)
{
    int err = clEnqueueCopyBuffer(commands, src, dst, 0, 0, N*sizeof(T), 0, NULL, NULL);
	CHKERR(err, "Unable to copy memory");
}

template <typename T>
void upload(cl_command_queue commands, cl_mem dst, T* src, int N)
{
    int err = clEnqueueWriteBuffer(commands, dst, CL_TRUE, 0, sizeof(T) * N, src, 0, NULL, NULL);
	CHKERR(err, "Unable to write memory to device");
}

template <typename T>
void download(cl_command_queue commands, T* dst, cl_mem src, int N)
{
    int err = clEnqueueReadBuffer(commands, src, CL_TRUE, 0, sizeof(T)*N, dst, 0, NULL, NULL);
	CHKERR(err, "Unable to read memory from device");
}

long getKernelSize()
{
    FILE* fp;
    long size;

    /*Determine the size of the file*/
    fp = fopen(KernelSourceFile, "r");
    size = 0;
    while(fgetc(fp) != EOF)
        size++;
    /*if(!fp)
        return 0;
    fseek(fp, 0, SEEK_END);

    size = ftell(fp) - begin;*/
    fclose(fp);

    return size + 1;
}

void getKernelSource(char* output, long size)
{
    FILE* fp;

    fp = fopen(KernelSourceFile, "r");
    if(!fp)
        return;

    fread(output, sizeof(char), size - 1, fp);
    output[size - 1] = 0;
    fclose(fp);
}

void dump(cl_command_queue commands, cl_mem variables, int nel, int nelr)
{
	float* h_variables = new float[nelr*NVAR];
	download(commands, h_variables, variables, nelr*NVAR);

	{
		std::ofstream file("density");
		file << nel << " " << nelr << std::endl;
		for(int i = 0; i < nel; i++) file << h_variables[i + VAR_DENSITY*nelr] << std::endl;
	}


	{
		std::ofstream file("momentum");
		file << nel << " " << nelr << std::endl;
		for(int i = 0; i < nel; i++)
		{
			for(int j = 0; j != NDIM; j++)
				file << h_variables[i + (VAR_MOMENTUM+j)*nelr] << " ";
			file << std::endl;
		}
	}

	{
		std::ofstream file("density_energy");
		file << nel << " " << nelr << std::endl;
		for(int i = 0; i < nel; i++) file << h_variables[i + VAR_DENSITY_ENERGY*nelr] << std::endl;
	}
	delete[] h_variables;
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
    cl_int err;

    size_t global_size;
    size_t local_size;

    cl_platform_id platform_id;
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel_compute_flux;
    cl_kernel kernel_compute_flux_contributions;
    cl_kernel kernel_compute_step_factor;
    cl_kernel kernel_time_step;
    cl_kernel kernel_initialize_variables;

    cl_mem ff_variable;
    cl_mem ff_fc_momentum_x;
    cl_mem ff_fc_momentum_y;
    cl_mem ff_fc_momentum_z;
    cl_mem ff_fc_density_energy;

    if (argc < 2)
	{
		printf("specify data file name\n");
		return 0;
	}
	const char* data_file_name = argv[1];

    // Retrieve an OpenCL platform
    err = clGetPlatformIDs(1, &platform_id, NULL);
    CHKERR(err, "Failed to get a platform!");

    // Connect to a compute device
    err = clGetDeviceIDs(platform_id, USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    CHKERR(err, "Failed to create a device group!");

    // Create a compute context
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");

    // Create a command queue
    commands = clCreateCommandQueue(context, device_id, 0, &err);
    CHKERR(err, "Failed to create a command queue!");



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

		// copy far field conditions to the gpu
		ff_variable = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float) * NVAR, h_ff_variable, &err);
		CHKERR(err, "Unable to allocate ff data");
		ff_fc_momentum_x = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float3), &h_ff_fc_momentum_x, &err);
		CHKERR(err, "Unable to allocate ff data");
		ff_fc_momentum_y = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float3), &h_ff_fc_momentum_y, &err);
		CHKERR(err, "Unable to allocate ff data");
		ff_fc_momentum_z = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float3), &h_ff_fc_momentum_z, &err);
		CHKERR(err, "Unable to allocate ff data");
		ff_fc_density_energy = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float3), &h_ff_fc_density_energy, &err);
		CHKERR(err, "Unable to allocate ff data");
	}
	int nel;
	int nelr;

	// read in domain geometry
	cl_mem areas;
	cl_mem elements_surrounding_elements;
	cl_mem normals;
    {
        std::ifstream file(data_file_name);

		file >> nel;

		nelr = block_length*((nel / block_length )+ std::min(1, nel % block_length));

		float* h_areas = new float[nelr];
		int* h_elements_surrounding_elements = new int[nelr*NNB];
		float* h_normals = new float[nelr*NDIM*NNB];


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
				for(int k = 0; k < NDIM; k++) h_normals[last + (j + k*NNB)*nelr] = h_normals[last + (j + k*NNB)*nelr];
			}
		}

		areas = alloc<float>(context, nelr);
		upload<float>(commands, areas, h_areas, nelr);

		elements_surrounding_elements = alloc<int>(context, nelr*NNB);
		upload<int>(commands, elements_surrounding_elements, h_elements_surrounding_elements, nelr*NNB);

		normals = alloc<float>(context, nelr*NDIM*NNB);
		upload<float>(commands, normals, h_normals, nelr*NDIM*NNB);

		delete[] h_areas;
		delete[] h_elements_surrounding_elements;
		delete[] h_normals;
    }

    // Get program source.
    long kernelSize = getKernelSize();
    char* kernelSource = new char[kernelSize];
    getKernelSource(kernelSource, kernelSize);

    // Create the compute program from the source buffer
    program = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, NULL, &err);
    CHKERR(err, "Failed to create a compute program!");

    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err == CL_BUILD_PROGRAM_FAILURE)
    {
        char *log;
        size_t logLen;
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLen);
        log = (char *) malloc(sizeof(char)*logLen);
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logLen, (void *) log, NULL);
        fprintf(stderr, "CL Error %d: Failed to build program! Log:\n%s", err, log);
        free(log);
        exit(1);
    }
    CHKERR(err, "Failed to build program!");
    delete[] kernelSource;

    // Create the compute kernel in the program we wish to run
    kernel_compute_flux = clCreateKernel(program, "compute_flux", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    // Create the reduce kernel in the program we wish to run
    kernel_compute_flux_contributions = clCreateKernel(program, "compute_flux_contributions", &err);
    CHKERR(err, "Failed to create a compute_flux_contributions kernel!");
    // Create the reduce kernel in the program we wish to run
    kernel_compute_step_factor = clCreateKernel(program, "compute_step_factor", &err);
    CHKERR(err, "Failed to create a compute_step_factor kernel!");
    // Create the reduce kernel in the program we wish to run
    kernel_time_step = clCreateKernel(program, "time_step", &err);
    CHKERR(err, "Failed to create a time_step kernel!");
    // Create the reduce kernel in the program we wish to run
    kernel_initialize_variables = clCreateKernel(program, "initialize_variables", &err);
    CHKERR(err, "Failed to create a initialize_variables kernel!");

	// Create arrays and set initial conditions
	cl_mem variables = alloc<cl_float>(context, nelr*NVAR);

    err = 0;
    err = clSetKernelArg(kernel_initialize_variables, 0, sizeof(int), &nelr);
    err |= clSetKernelArg(kernel_initialize_variables, 1, sizeof(cl_mem),&variables);
    err |= clSetKernelArg(kernel_initialize_variables, 2, sizeof(cl_mem),&ff_variable);
    CHKERR(err, "Failed to set kernel arguments!");
    // Get the maximum work group size for executing the kernel on the device
    //err = clGetKernelWorkGroupInfo(kernel_initialize_variables, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel_initialize_variables work group info!");
    local_size = 1;//std::min(local_size, (size_t)nelr);
    global_size = nelr;
    err = clEnqueueNDRangeKernel(commands, kernel_initialize_variables, 1, NULL, &global_size, NULL, 0, NULL, NULL);
    CHKERR(err, "Failed to execute kernel [kernel_initialize_variables]! 0");
    err = clFinish(commands);
    CHKERR(err, "Failed to execute kernel [kernel_test]!");


	cl_mem old_variables = alloc<float>(context, nelr*NVAR);
	cl_mem fluxes = alloc<float>(context, nelr*NVAR);
	cl_mem step_factors = alloc<float>(context, nelr);
    clFinish(commands);
	cl_mem fc_momentum_x = alloc<float>(context, nelr*NDIM);
	cl_mem fc_momentum_y = alloc<float>(context, nelr*NDIM);
	cl_mem fc_momentum_z = alloc<float>(context, nelr*NDIM);
	cl_mem fc_density_energy = alloc<float>(context, nelr*NDIM);
    clFinish(commands);

	// make sure all memory is floatly allocated before we start timing
    err = 0;
    err = clSetKernelArg(kernel_initialize_variables, 0, sizeof(int), &nelr);
    err |= clSetKernelArg(kernel_initialize_variables, 1, sizeof(cl_mem),&old_variables);
    err |= clSetKernelArg(kernel_initialize_variables, 2, sizeof(cl_mem),&ff_variable);
    CHKERR(err, "Failed to set kernel arguments!");
    // Get the maximum work group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(kernel_initialize_variables, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel_initialize_variables work group info!");
    err = clEnqueueNDRangeKernel(commands, kernel_initialize_variables, 1, NULL, &global_size, NULL, 0, NULL, NULL);
    CHKERR(err, "Failed to execute kernel [kernel_initialize_variables]! 1");
    clFinish(commands);

    err = 0;
    err = clSetKernelArg(kernel_initialize_variables, 0, sizeof(int), &nelr);
    err |= clSetKernelArg(kernel_initialize_variables, 1, sizeof(cl_mem),&fluxes);
    err |= clSetKernelArg(kernel_initialize_variables, 2, sizeof(cl_mem),&ff_variable);
    CHKERR(err, "Failed to set kernel arguments!");
    // Get the maximum work group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(kernel_compute_step_factor, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel_compute_step_factor work group info!");
    err = clEnqueueNDRangeKernel(commands, kernel_initialize_variables, 1, NULL, &global_size, NULL, 0, NULL, NULL);
    CHKERR(err, "Failed to execute kernel [kernel_initialize_variables]! 2");

    clFinish(commands);
    std::cout << "About to memcopy" << std::endl;
    err = clReleaseMemObject(step_factors);
    float temp[nelr];
    for(int i = 0; i < nelr; i++)
        temp[i] = 0;
    step_factors = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(float) * nelr, temp, &err);
    CHKERR(err, "Unable to memset step_factors");
	// make sure CUDA isn't still doing something before we start timing

    clFinish(commands);

	// these need to be computed the first time in order to compute time step
	std::cout << "Starting..." << std::endl;


//	unsigned int timer = 0;
//	CUT_SAFE_CALL( cutCreateTimer( &timer));
//	CUT_SAFE_CALL( cutStartTimer( timer));

	// Begin iterations
	for(int i = 0; i < iterations; i++)
	{
		copy<float>(commands, old_variables, variables, nelr*NVAR);

		// for the first iteration we compute the time step
        err = 0;
        err = clSetKernelArg(kernel_compute_step_factor, 0, sizeof(int), &nelr);
        err |= clSetKernelArg(kernel_compute_step_factor, 1, sizeof(cl_mem),&variables);
        err |= clSetKernelArg(kernel_compute_step_factor, 2, sizeof(cl_mem), &areas);
        err |= clSetKernelArg(kernel_compute_step_factor, 3, sizeof(cl_mem), &step_factors);
        CHKERR(err, "Failed to set kernel arguments!");
        // Get the maximum work group size for executing the kernel on the device
        err = clGetKernelWorkGroupInfo(kernel_compute_step_factor, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
        CHKERR(err, "Failed to retrieve kernel_compute_step_factor work group info!");
        err = clEnqueueNDRangeKernel(commands, kernel_compute_step_factor, 1, NULL, &global_size, NULL, 0, NULL, NULL);
        CHKERR(err, "Failed to execute kernel[kernel_compute_step_factor]!");
		clFinish(commands);

		for(int j = 0; j < RK; j++)
        {
            err = 0;
            err = clSetKernelArg(kernel_compute_flux_contributions, 0, sizeof(int), &nelr);
            err |= clSetKernelArg(kernel_compute_flux_contributions, 1, sizeof(cl_mem),&variables);
            err |= clSetKernelArg(kernel_compute_flux_contributions, 2, sizeof(cl_mem), &fc_momentum_x);
            err |= clSetKernelArg(kernel_compute_flux_contributions, 3, sizeof(cl_mem), &fc_momentum_y);
            err |= clSetKernelArg(kernel_compute_flux_contributions, 4, sizeof(cl_mem), &fc_momentum_z);
            err |= clSetKernelArg(kernel_compute_flux_contributions, 5, sizeof(cl_mem), &fc_density_energy);
            CHKERR(err, "Failed to set kernel arguments!");
            // Get the maximum work group size for executing the kernel on the device
            err = clGetKernelWorkGroupInfo(kernel_compute_flux_contributions, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
            CHKERR(err, "Failed to retrieve kernel_compute_flux_contributions work group info!");
            err = clEnqueueNDRangeKernel(commands, kernel_compute_flux_contributions, 1, NULL, &global_size, NULL, 0, NULL, NULL);
            CHKERR(err, "Failed to execute kernel [kernel_compute_flux_contributions]!");
			//compute_flux_contributions(nelr, variables, fc_momentum_x, fc_momentum_y, fc_momentum_z, fc_density_energy);
			clFinish(commands);

            err = 0;
            err = clSetKernelArg(kernel_compute_flux, 0, sizeof(int), &nelr);
            err |= clSetKernelArg(kernel_compute_flux, 1, sizeof(cl_mem), &elements_surrounding_elements);
            err |= clSetKernelArg(kernel_compute_flux, 2, sizeof(cl_mem), &normals);
            err |= clSetKernelArg(kernel_compute_flux, 3, sizeof(cl_mem), &variables);
            err |= clSetKernelArg(kernel_compute_flux, 4, sizeof(cl_mem), &fc_momentum_x);
            err |= clSetKernelArg(kernel_compute_flux, 5, sizeof(cl_mem), &fc_momentum_y);
            err |= clSetKernelArg(kernel_compute_flux, 6, sizeof(cl_mem), &fc_momentum_z);
            err |= clSetKernelArg(kernel_compute_flux, 7, sizeof(cl_mem), &fc_density_energy);
            err |= clSetKernelArg(kernel_compute_flux, 8, sizeof(cl_mem), &fluxes);
            err |= clSetKernelArg(kernel_compute_flux, 9, sizeof(cl_mem), &ff_variable);
            err |= clSetKernelArg(kernel_compute_flux, 10, sizeof(cl_mem), &ff_fc_momentum_x);
            err |= clSetKernelArg(kernel_compute_flux, 11, sizeof(cl_mem), &ff_fc_momentum_y);
            err |= clSetKernelArg(kernel_compute_flux, 12, sizeof(cl_mem), &ff_fc_momentum_z);
            err |= clSetKernelArg(kernel_compute_flux, 13, sizeof(cl_mem), &ff_fc_density_energy);
            CHKERR(err, "Failed to set kernel arguments!");
            // Get the maximum work group size for executing the kernel on the device
            err = clGetKernelWorkGroupInfo(kernel_compute_flux, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
            CHKERR(err, "Failed to retrieve kernel_compute_flux work group info!");
            err = clEnqueueNDRangeKernel(commands, kernel_compute_flux, 1, NULL, &global_size, NULL, 0, NULL, NULL);
            CHKERR(err, "Failed to execute kernel [kernel_compute_flux]!");
			clFinish(commands);

            err = 0;
            err = clSetKernelArg(kernel_time_step, 0, sizeof(int), &j);
            err |= clSetKernelArg(kernel_time_step, 1, sizeof(int), &nelr);
            err |= clSetKernelArg(kernel_time_step, 2, sizeof(cl_mem), &old_variables);
            err |= clSetKernelArg(kernel_time_step, 3, sizeof(cl_mem), &variables);
            err |= clSetKernelArg(kernel_time_step, 4, sizeof(cl_mem), &step_factors);
            err |= clSetKernelArg(kernel_time_step, 5, sizeof(cl_mem), &fluxes);
            CHKERR(err, "Failed to set kernel arguments!");
            // Get the maximum work group size for executing the kernel on the device
            err = clGetKernelWorkGroupInfo(kernel_time_step, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
            CHKERR(err, "Failed to retrieve kernel_time_step work group info!");
            err = clEnqueueNDRangeKernel(commands, kernel_time_step, 1, NULL, &global_size, NULL, 0, NULL, NULL);
            CHKERR(err, "Failed to execute kernel [kernel_time_step]!");
			clFinish(commands);
		}
	}

//	cudaThreadSynchronize();
    clFinish(commands);
    std::cout << "Finished" << std::endl;

//	CUT_SAFE_CALL( cutStopTimer(timer) );

//	std::cout  << (cutGetAverageTimerValue(timer)/1000.0)  / iterations << " seconds per iteration" << std::endl;

	std::cout << "Saving solution..." << std::endl;
	dump(commands, variables, nel, nelr);
	std::cout << "Saved solution..." << std::endl;


	std::cout << "Cleaning up..." << std::endl;

    clReleaseProgram(program);
    clReleaseKernel(kernel_compute_flux);
    clReleaseKernel(kernel_compute_flux_contributions);
    clReleaseKernel(kernel_compute_step_factor);
    clReleaseKernel(kernel_time_step);
    clReleaseKernel(kernel_initialize_variables);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

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

    return 0;
}
