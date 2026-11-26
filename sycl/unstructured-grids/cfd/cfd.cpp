#include <stdio.h>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string.h>
#include "dwarf_sycl.hpp"
#include <cmath>
/* host copies of the mesh, kept so every timed repetition can re-upload it */
static float *g_areas = 0, *g_normals = 0; static int *g_ese = 0;
#include <malloc.h>

#define AOCL_ALIGNMENT 64

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

const char *KernelSourceFile = "cfd_kernel.cl";


/*
 * Generic functions
 */
	template <typename T>
T* alloc(sycl::queue &q, int N) { return dsycl::alloc<T>(q, N); }

	template <typename T>
void dealloc(sycl::queue &q, T* array) { dsycl::release(q, array); }

	template <typename T>
void copy(sycl::queue &q, T* dst, T* src, int N) { q.memcpy(dst, src, N*sizeof(T)).wait_and_throw(); }

	template <typename T>
void upload(sycl::queue &q, T* dst, T* src, int N) { dsycl::write(q, dst, src, sizeof(T) * N); }

	template <typename T>
void download(sycl::queue &q, T* dst, T* src, int N) { dsycl::read(q, dst, src, sizeof(T)*N); }

void dump(sycl::queue &q, float* variables, int nel, int nelr)
{
	float* h_variables;
    h_variables  = (float*)  memalign(AOCL_ALIGNMENT,nelr*NVAR*sizeof(float));
	download(q, h_variables, variables, nelr*NVAR);

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


/* ---- kernels, translated from cfd_kernel.cl (bodies verbatim). The device float3 values are only accessed through
 * .x/.y/.z, so the 12-byte host struct is used on both sides (as the host-side ff buffers already are). ---- */
namespace dev {
static constexpr float KGAMMA = 1.4f;   /* the kernels use GAMMA = 1.4f (the host code uses the double 1.4) */
using sycl::sqrt;
inline void compute_velocity(float density, float3 momentum, float3* velocity) {
	(*velocity).x = momentum.x / density;
	(*velocity).y = momentum.y / density;
	(*velocity).z = momentum.z / density;
}
inline float compute_speed_of_sound(float density, float pressure) {
	return sqrt(KGAMMA*pressure/density);
}
inline float compute_speed_sqd(float3 velocity) {
	return velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
}
inline float compute_pressure(float density, float density_energy, float speed_sqd) {
	return (KGAMMA-1.0f)*(density_energy - 0.5f*density*speed_sqd);
}
} // namespace dev
static void k_compute_flux(sycl::queue &q, size_t global_size, int nelr, int* elements_surrounding_elements, float* normals, float* variables, float* fc_momentum_x, float* fc_momentum_y, float* fc_momentum_z, float* fc_density_energy, float* fluxes, float* ff_variable, float3* ff_fc_momentum_x, float3* ff_fc_momentum_y, float3* ff_fc_momentum_z, float3* ff_fc_density_energy)
{
	q.parallel_for<class cfd_compute_flux>(sycl::nd_range<1>(global_size, DWARF_LS), [=](sycl::nd_item<1> it) {
		const float smoothing_coefficient = 0.2f;
		const int i = it.get_global_id(0);
		if (i >= nelr) return;  /* [FIX-OCL-34] guard padded work-items (global = roundup(nelr,256)) */

		int j, nb;
		float3 normal; float normal_len;
		float factor;

		float density_i = variables[i + VAR_DENSITY*nelr];
		float3 momentum_i;
		momentum_i.x = variables[i + (VAR_MOMENTUM+0)*nelr];
		momentum_i.y = variables[i + (VAR_MOMENTUM+1)*nelr];
		momentum_i.z = variables[i + (VAR_MOMENTUM+2)*nelr];

		float density_energy_i = variables[i + VAR_DENSITY_ENERGY*nelr];

		float3 velocity_i;            		   dev::compute_velocity(density_i, momentum_i, &velocity_i);
		float speed_sqd_i                          = dev::compute_speed_sqd(velocity_i);
		float speed_i                              = sycl::sqrt(speed_sqd_i);
		float pressure_i                           = dev::compute_pressure(density_i, density_energy_i, speed_sqd_i);
		float speed_of_sound_i                     = dev::compute_speed_of_sound(density_i, pressure_i);
		float3 fc_i_momentum_x, fc_i_momentum_y, fc_i_momentum_z;
		float3 fc_i_density_energy;

		fc_i_momentum_x.x = fc_momentum_x[i + 0*nelr];
		fc_i_momentum_x.y = fc_momentum_x[i + 1*nelr];
		fc_i_momentum_x.z = fc_momentum_x[i + 2*nelr];

		fc_i_momentum_y.x = fc_momentum_y[i + 0*nelr];
		fc_i_momentum_y.y = fc_momentum_y[i + 1*nelr];
		fc_i_momentum_y.z = fc_momentum_y[i + 2*nelr];

		fc_i_momentum_z.x = fc_momentum_z[i + 0*nelr];
		fc_i_momentum_z.y = fc_momentum_z[i + 1*nelr];
		fc_i_momentum_z.z = fc_momentum_z[i + 2*nelr];

		fc_i_density_energy.x = fc_density_energy[i + 0*nelr];
		fc_i_density_energy.y = fc_density_energy[i + 1*nelr];
		fc_i_density_energy.z = fc_density_energy[i + 2*nelr];

		float flux_i_density = 0.0f;
		float3 flux_i_momentum;
		flux_i_momentum.x = 0.0f;
		flux_i_momentum.y = 0.0f;
		flux_i_momentum.z = 0.0f;
		float flux_i_density_energy = 0.0f;

		float3 velocity_nb;
		float density_nb, density_energy_nb;
		float3 momentum_nb;
		float3 fc_nb_momentum_x, fc_nb_momentum_y, fc_nb_momentum_z;
		float3 fc_nb_density_energy;
		float speed_sqd_nb, speed_of_sound_nb, pressure_nb;

		for(j = 0; j < NNB; j++)
		{
			nb = elements_surrounding_elements[i + j*nelr];
			normal.x = normals[i + (j + 0*NNB)*nelr];
			normal.y = normals[i + (j + 1*NNB)*nelr];
			normal.z = normals[i + (j + 2*NNB)*nelr];
			normal_len = sycl::sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);

			if(nb >= 0) 	// a legitimate neighbor
			{
				density_nb = variables[nb + VAR_DENSITY*nelr];
				momentum_nb.x = variables[nb + (VAR_MOMENTUM+0)*nelr];
				momentum_nb.y = variables[nb + (VAR_MOMENTUM+1)*nelr];
				momentum_nb.z = variables[nb + (VAR_MOMENTUM+2)*nelr];
				density_energy_nb = variables[nb + VAR_DENSITY_ENERGY*nelr];
				dev::compute_velocity(density_nb, momentum_nb, &velocity_nb);
				speed_sqd_nb                      = dev::compute_speed_sqd(velocity_nb);
				pressure_nb                       = dev::compute_pressure(density_nb, density_energy_nb, speed_sqd_nb);
				speed_of_sound_nb                 = dev::compute_speed_of_sound(density_nb, pressure_nb);

				fc_nb_momentum_x.x = fc_momentum_x[nb + 0*nelr];
				fc_nb_momentum_x.y = fc_momentum_x[nb + 1*nelr];
				fc_nb_momentum_x.z = fc_momentum_x[nb + 2*nelr];

				fc_nb_momentum_y.x = fc_momentum_y[nb + 0*nelr];
				fc_nb_momentum_y.y = fc_momentum_y[nb + 1*nelr];
				fc_nb_momentum_y.z = fc_momentum_y[nb + 2*nelr];

				fc_nb_momentum_z.x = fc_momentum_z[nb + 0*nelr];
				fc_nb_momentum_z.y = fc_momentum_z[nb + 1*nelr];
				fc_nb_momentum_z.z = fc_momentum_z[nb + 2*nelr];

				fc_nb_density_energy.x = fc_density_energy[nb + 0*nelr];
				fc_nb_density_energy.y = fc_density_energy[nb + 1*nelr];
				fc_nb_density_energy.z = fc_density_energy[nb + 2*nelr];

				// artificial viscosity
				factor = -normal_len*smoothing_coefficient*0.5f*(speed_i + sycl::sqrt(speed_sqd_nb) + speed_of_sound_i + speed_of_sound_nb);
				flux_i_density += factor*(density_i-density_nb);
				flux_i_density_energy += factor*(density_energy_i-density_energy_nb);
				flux_i_momentum.x += factor*(momentum_i.x-momentum_nb.x);
				flux_i_momentum.y += factor*(momentum_i.y-momentum_nb.y);
				flux_i_momentum.z += factor*(momentum_i.z-momentum_nb.z);

				// accumulate cell-centered fluxes
				factor = 0.5f*normal.x;
				flux_i_density += factor*(momentum_nb.x+momentum_i.x);
				flux_i_density_energy += factor*(fc_nb_density_energy.x+fc_i_density_energy.x);
				flux_i_momentum.x += factor*(fc_nb_momentum_x.x+fc_i_momentum_x.x);
				flux_i_momentum.y += factor*(fc_nb_momentum_y.x+fc_i_momentum_y.x);
				flux_i_momentum.z += factor*(fc_nb_momentum_z.x+fc_i_momentum_z.x);

				factor = 0.5f*normal.y;
				flux_i_density += factor*(momentum_nb.y+momentum_i.y);
				flux_i_density_energy += factor*(fc_nb_density_energy.y+fc_i_density_energy.y);
				flux_i_momentum.x += factor*(fc_nb_momentum_x.y+fc_i_momentum_x.y);
				flux_i_momentum.y += factor*(fc_nb_momentum_y.y+fc_i_momentum_y.y);
				flux_i_momentum.z += factor*(fc_nb_momentum_z.y+fc_i_momentum_z.y);

				factor = 0.5f*normal.z;
				flux_i_density += factor*(momentum_nb.z+momentum_i.z);
				flux_i_density_energy += factor*(fc_nb_density_energy.z+fc_i_density_energy.z);
				flux_i_momentum.x += factor*(fc_nb_momentum_x.z+fc_i_momentum_x.z);
				flux_i_momentum.y += factor*(fc_nb_momentum_y.z+fc_i_momentum_y.z);
				flux_i_momentum.z += factor*(fc_nb_momentum_z.z+fc_i_momentum_z.z);
			}
			else if(nb == -1)	// a wing boundary
			{
				flux_i_momentum.x += normal.x*pressure_i;
				flux_i_momentum.y += normal.y*pressure_i;
				flux_i_momentum.z += normal.z*pressure_i;
			}
			else if(nb == -2) // a far field boundary
			{
				factor = 0.5f*normal.x;
				flux_i_density += factor*(ff_variable[VAR_MOMENTUM+0]+momentum_i.x);
				flux_i_density_energy += factor*(ff_fc_density_energy[0].x+fc_i_density_energy.x);
				flux_i_momentum.x += factor*(ff_fc_momentum_x[0].x + fc_i_momentum_x.x);
				flux_i_momentum.y += factor*(ff_fc_momentum_y[0].x + fc_i_momentum_y.x);
				flux_i_momentum.z += factor*(ff_fc_momentum_z[0].x + fc_i_momentum_z.x);

				factor = 0.5f*normal.y;
				flux_i_density += factor*(ff_variable[VAR_MOMENTUM+1]+momentum_i.y);
				flux_i_density_energy += factor*(ff_fc_density_energy[0].y+fc_i_density_energy.y);
				flux_i_momentum.x += factor*(ff_fc_momentum_x[0].y + fc_i_momentum_x.y);
				flux_i_momentum.y += factor*(ff_fc_momentum_y[0].y + fc_i_momentum_y.y);
				flux_i_momentum.z += factor*(ff_fc_momentum_z[0].y + fc_i_momentum_z.y);

				factor = 0.5f*normal.z;
				flux_i_density += factor*(ff_variable[VAR_MOMENTUM+2]+momentum_i.z);
				flux_i_density_energy += factor*(ff_fc_density_energy[0].z+fc_i_density_energy.z);
				flux_i_momentum.x += factor*(ff_fc_momentum_x[0].z + fc_i_momentum_x.z);
				flux_i_momentum.y += factor*(ff_fc_momentum_y[0].z + fc_i_momentum_y.z);
				flux_i_momentum.z += factor*(ff_fc_momentum_z[0].z + fc_i_momentum_z.z);

			}
		}

		fluxes[i + VAR_DENSITY*nelr] = flux_i_density;
		fluxes[i + (VAR_MOMENTUM+0)*nelr] = flux_i_momentum.x;
		fluxes[i + (VAR_MOMENTUM+1)*nelr] = flux_i_momentum.y;
		fluxes[i + (VAR_MOMENTUM+2)*nelr] = flux_i_momentum.z;
		fluxes[i + VAR_DENSITY_ENERGY*nelr] = flux_i_density_energy;
}).wait_and_throw();
}

static void k_compute_flux_contributions(sycl::queue &q, size_t global_size, int nelr, float* variables, float* fc_momentum_x, float* fc_momentum_y, float* fc_momentum_z, float* fc_density_energy)
{
	q.parallel_for<class cfd_compute_flux_contributions>(sycl::nd_range<1>(global_size, DWARF_LS), [=](sycl::nd_item<1> it) {
		const int i = it.get_global_id(0);
		if (i >= nelr) return;  /* [FIX-OCL-34] guard padded work-items (global = roundup(nelr,256)) */

		float density_i = variables[i + VAR_DENSITY*nelr];
		float3 momentum_i;
		momentum_i.x = variables[i + (VAR_MOMENTUM+0)*nelr];
		momentum_i.y = variables[i + (VAR_MOMENTUM+1)*nelr];
		momentum_i.z = variables[i + (VAR_MOMENTUM+2)*nelr];
		float density_energy_i = variables[i + VAR_DENSITY_ENERGY*nelr];

		float3 velocity_i;             				dev::compute_velocity(density_i, momentum_i, &velocity_i);
		float speed_sqd_i                          = dev::compute_speed_sqd(velocity_i);
		float speed_i                              = sycl::sqrt(speed_sqd_i);
		float pressure_i                           = dev::compute_pressure(density_i, density_energy_i, speed_sqd_i);
		float speed_of_sound_i                     = dev::compute_speed_of_sound(density_i, pressure_i);
		float3 fc_i_momentum_x, fc_i_momentum_y, fc_i_momentum_z;
		float3 fc_i_density_energy;
		::compute_flux_contribution(&density_i, &momentum_i, &density_energy_i, pressure_i, &velocity_i, &fc_i_momentum_x, &fc_i_momentum_y, &fc_i_momentum_z, &fc_i_density_energy);

		fc_momentum_x[i + 0*nelr] = fc_i_momentum_x.x;
		fc_momentum_x[i + 1*nelr] = fc_i_momentum_x.y;
		fc_momentum_x[i + 2*nelr] = fc_i_momentum_x.z;

		fc_momentum_y[i + 0*nelr] = fc_i_momentum_y.x;
		fc_momentum_y[i + 1*nelr] = fc_i_momentum_y.y;
		fc_momentum_y[i + 2*nelr] = fc_i_momentum_y.z;


		fc_momentum_z[i + 0*nelr] = fc_i_momentum_z.x;
		fc_momentum_z[i + 1*nelr] = fc_i_momentum_z.y;
		fc_momentum_z[i + 2*nelr] = fc_i_momentum_z.z;

		fc_density_energy[i + 0*nelr] = fc_i_density_energy.x;
		fc_density_energy[i + 1*nelr] = fc_i_density_energy.y;
		fc_density_energy[i + 2*nelr] = fc_i_density_energy.z;
}).wait_and_throw();
}

static void k_compute_step_factor(sycl::queue &q, size_t global_size, int nelr, float* variables, float* areas, float* step_factors)
{
	q.parallel_for<class cfd_compute_step_factor>(sycl::nd_range<1>(global_size, DWARF_LS), [=](sycl::nd_item<1> it) {
		const int i = it.get_global_id(0);
		if (i >= nelr) return;  /* [FIX-OCL-34] guard padded work-items (global = roundup(nelr,256)) */

		float density = variables[i + VAR_DENSITY*nelr];
		float3 momentum;
		momentum.x = variables[i + (VAR_MOMENTUM+0)*nelr];
		momentum.y = variables[i + (VAR_MOMENTUM+1)*nelr];
		momentum.z = variables[i + (VAR_MOMENTUM+2)*nelr];

		float density_energy = variables[i + VAR_DENSITY_ENERGY*nelr];

		float3 velocity;       dev::compute_velocity(density, momentum, &velocity);
		float speed_sqd      = dev::compute_speed_sqd(velocity);
		float pressure       = dev::compute_pressure(density, density_energy, speed_sqd);
		float speed_of_sound = dev::compute_speed_of_sound(density, pressure);

		// dt = float(0.5) * sycl::sqrt(areas[i]) /  (||v|| + c).... but when we do time stepping, this later would need to be divided by the area, so we just do it all at once
		step_factors[i] = 0.5f / (sycl::sqrt(areas[i]) * (sycl::sqrt(speed_sqd) + speed_of_sound));
}).wait_and_throw();
}

static void k_time_step(sycl::queue &q, size_t global_size, int j, int nelr, float* old_variables, float* variables, float* step_factors, float* fluxes)
{
	q.parallel_for<class cfd_time_step>(sycl::nd_range<1>(global_size, DWARF_LS), [=](sycl::nd_item<1> it) {
		const int i = it.get_global_id(0);
		if (i >= nelr) return;  /* [FIX-OCL-34] guard padded work-items (global = roundup(nelr,256)) */

		float factor = step_factors[i]/(RK+1-j);

		variables[i + VAR_DENSITY*nelr] = old_variables[i + VAR_DENSITY*nelr] + factor*fluxes[i + VAR_DENSITY*nelr];
		variables[i + VAR_DENSITY_ENERGY*nelr] = old_variables[i + VAR_DENSITY_ENERGY*nelr] + factor*fluxes[i + VAR_DENSITY_ENERGY*nelr];
		variables[i + (VAR_MOMENTUM+0)*nelr] = old_variables[i + (VAR_MOMENTUM+0)*nelr] + factor*fluxes[i + (VAR_MOMENTUM+0)*nelr];
		variables[i + (VAR_MOMENTUM+1)*nelr] = old_variables[i + (VAR_MOMENTUM+1)*nelr] + factor*fluxes[i + (VAR_MOMENTUM+1)*nelr];
		variables[i + (VAR_MOMENTUM+2)*nelr] = old_variables[i + (VAR_MOMENTUM+2)*nelr] + factor*fluxes[i + (VAR_MOMENTUM+2)*nelr];
}).wait_and_throw();
}

static void k_initialize_variables(sycl::queue &q, size_t global_size, int nelr, float* variables, float* ff_variable)
{
	q.parallel_for<class cfd_initialize_variables>(sycl::nd_range<1>(global_size, DWARF_LS), [=](sycl::nd_item<1> it) {
		unsigned int i = it.get_global_id(0);
		if (i >= (unsigned int)nelr) return;  /* [FIX-OCL-34] guard padded work-items */
		for(int j = 0; j < NVAR; j++)
			variables[i + j*nelr] = ff_variable[j];
}).wait_and_throw();
}

int main(int argc, char** argv)
{
	dsycl::parse_args(&argc, &argv);

	int err;

	size_t global_size;

	float *ff_variable;
	float3 *ff_fc_momentum_x, *ff_fc_momentum_y, *ff_fc_momentum_z, *ff_fc_density_energy;

	if (argc < 2)
	{
		printf("Usage ./cfd <data input file>\n");
		return 0;
	}


	const char* data_file_name = argv[1];
	sycl::queue q = dsycl::make_queue();


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
		ff_variable = alloc<float>(q, NVAR);                 upload<float>(q, ff_variable, h_ff_variable, NVAR);
		ff_fc_momentum_x = alloc<float3>(q, 1);              upload<float3>(q, ff_fc_momentum_x, &h_ff_fc_momentum_x, 1);
		ff_fc_momentum_y = alloc<float3>(q, 1);              upload<float3>(q, ff_fc_momentum_y, &h_ff_fc_momentum_y, 1);
		ff_fc_momentum_z = alloc<float3>(q, 1);              upload<float3>(q, ff_fc_momentum_z, &h_ff_fc_momentum_z, 1);
		ff_fc_density_energy = alloc<float3>(q, 1);          upload<float3>(q, ff_fc_density_energy, &h_ff_fc_density_energy, 1);
	}
	int nel;
	int nelr;

	// read in domain geometry
	float* areas;
	int* elements_surrounding_elements;
	float* normals;
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

		areas = alloc<float>(q, nelr);
		upload<float>(q, areas, h_areas, nelr);

		elements_surrounding_elements = alloc<int>(q, nelr*NNB);
		upload<int>(q, elements_surrounding_elements, h_elements_surrounding_elements, nelr*NNB);

		normals = alloc<float>(q, nelr*NDIM*NNB);
		upload<float>(q, normals, h_normals, nelr*NDIM*NNB);

		g_areas = h_areas; g_ese = h_elements_surrounding_elements; g_normals = h_normals;  // kept for the timed re-uploads
	}

	// Create arrays and set initial conditions
	float *variables = alloc<float>(q, nelr*NVAR);
	float *old_variables = alloc<float>(q, nelr*NVAR);
	float *fluxes = alloc<float>(q, nelr*NVAR);
	float *step_factors = alloc<float>(q, nelr);
	float *fc_momentum_x = alloc<float>(q, nelr*NDIM);
	float *fc_momentum_y = alloc<float>(q, nelr*NDIM);
	float *fc_momentum_z = alloc<float>(q, nelr*NDIM);
	float *fc_density_energy = alloc<float>(q, nelr*NDIM);

	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */

	global_size = DWARF_ROUNDUP(nelr);  /* [FIX-OCL-34] pinned launch: local DWARF_LS, global roundup(nelr) */
	k_initialize_variables(q, global_size, nelr, variables, ff_variable);
	k_initialize_variables(q, global_size, nelr, old_variables, ff_variable);
	k_initialize_variables(q, global_size, nelr, fluxes, ff_variable);
	std::cout << "About to memcopy" << std::endl;
	std::vector<float> temp(nelr, 0.0f);
	upload<float>(q, step_factors, temp.data(), nelr);

	// these need to be computed the first time in order to compute time step
	std::cout << "Starting..." << std::endl;

	/* common timer: repetitions around [mesh upload, variable initialisation, the solver iterations, download of the
	 * solution]; JIT/build is above, writing the output files is below. */
	dwarf_init("cfd", "unstructured-grids", "sycl");
	dwarf_set_size(nel, nelr);
	float *h_result = (float*) malloc(sizeof(float) * nelr * NVAR);
	for(int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	dwarf_region_begin();
	upload<float>(q, areas, g_areas, nelr);
	upload<int>(q, elements_surrounding_elements, g_ese, nelr*NNB);
	upload<float>(q, normals, g_normals, nelr*NDIM*NNB);
	{
		float *targets[3] = { variables, old_variables, fluxes };
		for(int t = 0; t < 3; t++)
			k_initialize_variables(q, global_size, nelr, targets[t], ff_variable);
		upload<float>(q, step_factors, temp.data(), nelr);   // step_factors = 0
	}

	// Begin iterations
	for(int i = 0; i < iterations; i++)
	{
		copy<float>(q, old_variables, variables, nelr*NVAR);

		// for the first iteration we compute the time step
		k_compute_step_factor(q, global_size, nelr, variables, areas, step_factors);
		for(int j = 0; j < RK; j++)
		{
			k_compute_flux_contributions(q, global_size, nelr, variables, fc_momentum_x, fc_momentum_y, fc_momentum_z, fc_density_energy);
			k_compute_flux(q, global_size, nelr, elements_surrounding_elements, normals, variables, fc_momentum_x, fc_momentum_y, fc_momentum_z,
					fc_density_energy, fluxes, ff_variable, ff_fc_momentum_x, ff_fc_momentum_y, ff_fc_momentum_z, ff_fc_density_energy);
			k_time_step(q, global_size, j, nelr, old_variables, variables, step_factors, fluxes);
		}
	}

	download<float>(q, h_result, variables, nelr*NVAR);
	dwarf_region_end(rep);
	}
	std::cout << "Finished" << std::endl;
	std::cout << "Saving solution..." << std::endl;
	dump(q, variables, nel, nelr);
	std::cout << "Saved solution..." << std::endl;
	std::cout << "Cleaning up..." << std::endl;

	dealloc<float>(q, areas);
	dealloc<int>(q, elements_surrounding_elements);
	dealloc<float>(q, normals);
	dealloc<float>(q, variables);
	dealloc<float>(q, old_variables);
	dealloc<float>(q, fluxes);
	dealloc<float>(q, step_factors);
	dealloc<float>(q, fc_momentum_x);
	dealloc<float>(q, fc_momentum_y);
	dealloc<float>(q, fc_momentum_z);
	dealloc<float>(q, fc_density_energy);

	std::cout << "Done..." << std::endl;
	/* [FIX-OCL-11] validity check (none before): the solution must be finite with positive density */
	{
		long bad = 0;
		for(int i = 0; i < nel; i++) { float d = h_result[i + VAR_DENSITY*nelr]; if(!(std::isfinite(d) && d > 0.0f)) bad++; }
		printf("%s: %ld cells with non-finite or non-positive density\n", bad ? "FAILURE" : "SUCCESS", bad);
		return dwarf_report(bad == 0);
	}
}
