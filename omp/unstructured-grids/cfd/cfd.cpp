/* Unstructured-grid dwarf: CFD (Euler solver) -- OpenMP target-offload translation of
 * opencl/unstructured-grids/cfd/cfd.cpp + cfd_kernel.cl (Rodinia CFD, Corrigan et al.).
 *
 * Same host flow: read mesh (.domn file), far-field conditions, upload areas / neighbour list /
 * normals / far-field constants, initialise variables / old_variables / fluxes with the far-field
 * state, zero the step factors, then `iterations` x { old = variables; compute_step_factor;
 * RK x { compute_flux_contributions; compute_flux; time_step } }, and finally download the
 * variables and write density.txt / momentum.txt / density_energy.txt. The five kernels are the
 * OpenCL kernels (one work-item per element -> one parallel iteration per element).
 * Usage: cfd <data input file>.   Timed region: all uploads, initialisation kernels, the
 * iteration loop, the final download of `variables`.
 * Deviations: padded elements (nel..nelr-1) get zero normals (the OpenCL code leaves them
 * uninitialised via a self-assignment); they never affect the first nel results.
 * Validation (added): serial host evaluation of the same kernels, rel. tol 1e-3 (fma differences).
 */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "dwarf_timer.h"

typedef struct { float x; float y; float z; } float3;

#define GAMMA 1.4f
#define iterations 1
#define NDIM 3
#define NNB 4
#define RK 3  // 3rd order RK
#define ff_mach 1.2
#define deg_angle_of_attack 0.0f
#define VAR_DENSITY 0
#define VAR_MOMENTUM 1
#define VAR_DENSITY_ENERGY (VAR_MOMENTUM + NDIM)
#define NVAR (VAR_DENSITY_ENERGY + 1)
#ifndef block_length
#define block_length 128
#endif

#pragma omp declare target
void compute_flux_contribution(float* density, float3* momentum, float* density_energy,
		float pressure, float3* velocity,
		float3* fc_momentum_x, float3* fc_momentum_y, float3* fc_momentum_z,
		float3* fc_density_energy)
{
	(*fc_momentum_x).x = (*velocity).x*(*momentum).x + pressure;
	(*fc_momentum_x).y = (*velocity).x*(*momentum).y;
	(*fc_momentum_x).z = (*velocity).x*(*momentum).z;


	(*fc_momentum_y).x = (*fc_momentum_x).y;
	(*fc_momentum_y).y = (*velocity).y*(*momentum).y + pressure;
	(*fc_momentum_y).z = (*velocity).y*(*momentum).z;

	(*fc_momentum_z).x = (*fc_momentum_x).z;
	(*fc_momentum_z).y = (*fc_momentum_y).z;
	(*fc_momentum_z).z = (*velocity).z*(*momentum).z + pressure;

	float de_p = *density_energy+pressure;
	(*fc_density_energy).x = (*velocity).x*de_p;
	(*fc_density_energy).y = (*velocity).y*de_p;
	(*fc_density_energy).z = (*velocity).z*de_p;
}

void compute_velocity(float density, float3 momentum, float3* velocity)
{
	(*velocity).x = momentum.x / density;
	(*velocity).y = momentum.y / density;
	(*velocity).z = momentum.z / density;
}

float compute_speed_of_sound(float density, float pressure)
{
	return sqrtf(GAMMA*pressure/density);
}

float compute_speed_sqd(float3 velocity)
{
	return velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
}

float compute_pressure(float density, float density_energy, float speed_sqd)
{
	return (GAMMA-1.0f)*(density_energy - 0.5f*density*speed_sqd);
}

#pragma omp end declare target

/* ---- the five kernels of cfd_kernel.cl, one loop iteration = one work-item ---- */
#define KMAP_FLUX map(present, alloc : ese[0 : nelr * NNB], normals[0 : nelr * NDIM * NNB], variables[0 : nelr * NVAR], \
    fc_momentum_x[0 : nelr * NDIM], fc_momentum_y[0 : nelr * NDIM], fc_momentum_z[0 : nelr * NDIM], \
    fc_density_energy[0 : nelr * NDIM], fluxes[0 : nelr * NVAR], ff_variable[0 : NVAR], ff_fc_momentum_x[0 : 1], \
    ff_fc_momentum_y[0 : 1], ff_fc_momentum_z[0 : 1], ff_fc_density_energy[0 : 1])
static void compute_flux(int nelr, int *ese, float *normals, float *variables, float *fc_momentum_x,
                         float *fc_momentum_y, float *fc_momentum_z, float *fc_density_energy, float *fluxes,
                         float *ff_variable, float3 *ff_fc_momentum_x, float3 *ff_fc_momentum_y,
                         float3 *ff_fc_momentum_z, float3 *ff_fc_density_energy) {
#pragma omp target teams distribute parallel for KMAP_FLUX
  for (int i = 0; i < nelr; i++) {
    int *elements_surrounding_elements = ese;

	const float smoothing_coefficient = 0.2f;
	

	int j, nb;
	float3 normal; float normal_len;
	float factor;

	float density_i = variables[i + VAR_DENSITY*nelr];
	float3 momentum_i;
	momentum_i.x = variables[i + (VAR_MOMENTUM+0)*nelr];
	momentum_i.y = variables[i + (VAR_MOMENTUM+1)*nelr];
	momentum_i.z = variables[i + (VAR_MOMENTUM+2)*nelr];

	float density_energy_i = variables[i + VAR_DENSITY_ENERGY*nelr];

	float3 velocity_i;            		   compute_velocity(density_i, momentum_i, &velocity_i);
	float speed_sqd_i                          = compute_speed_sqd(velocity_i);
	float speed_i                              = sqrtf(speed_sqd_i);
	float pressure_i                           = compute_pressure(density_i, density_energy_i, speed_sqd_i);
	float speed_of_sound_i                     = compute_speed_of_sound(density_i, pressure_i);
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
		normal_len = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);

		if(nb >= 0) 	// a legitimate neighbor
		{
			density_nb = variables[nb + VAR_DENSITY*nelr];
			momentum_nb.x = variables[nb + (VAR_MOMENTUM+0)*nelr];
			momentum_nb.y = variables[nb + (VAR_MOMENTUM+1)*nelr];
			momentum_nb.z = variables[nb + (VAR_MOMENTUM+2)*nelr];
			density_energy_nb = variables[nb + VAR_DENSITY_ENERGY*nelr];
			compute_velocity(density_nb, momentum_nb, &velocity_nb);
			speed_sqd_nb                      = compute_speed_sqd(velocity_nb);
			pressure_nb                       = compute_pressure(density_nb, density_energy_nb, speed_sqd_nb);
			speed_of_sound_nb                 = compute_speed_of_sound(density_nb, pressure_nb);

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
			factor = -normal_len*smoothing_coefficient*0.5f*(speed_i + sqrtf(speed_sqd_nb) + speed_of_sound_i + speed_of_sound_nb);
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
  }
}

static void compute_flux_contributions(int nelr, float *variables, float *fc_momentum_x, float *fc_momentum_y,
                                       float *fc_momentum_z, float *fc_density_energy) {
#pragma omp target teams distribute parallel for map(present, alloc : variables[0 : nelr * NVAR], fc_momentum_x[0 : nelr * NDIM], fc_momentum_y[0 : nelr * NDIM], fc_momentum_z[0 : nelr * NDIM], fc_density_energy[0 : nelr * NDIM])
  for (int i = 0; i < nelr; i++) {

	

	float density_i = variables[i + VAR_DENSITY*nelr];
	float3 momentum_i;
	momentum_i.x = variables[i + (VAR_MOMENTUM+0)*nelr];
	momentum_i.y = variables[i + (VAR_MOMENTUM+1)*nelr];
	momentum_i.z = variables[i + (VAR_MOMENTUM+2)*nelr];
	float density_energy_i = variables[i + VAR_DENSITY_ENERGY*nelr];

	float3 velocity_i;             				compute_velocity(density_i, momentum_i, &velocity_i);
	float speed_sqd_i                          = compute_speed_sqd(velocity_i);
	float speed_i                              = sqrtf(speed_sqd_i);
	float pressure_i                           = compute_pressure(density_i, density_energy_i, speed_sqd_i);
	float speed_of_sound_i                     = compute_speed_of_sound(density_i, pressure_i);
	float3 fc_i_momentum_x, fc_i_momentum_y, fc_i_momentum_z;
	float3 fc_i_density_energy;
	compute_flux_contribution(&density_i, &momentum_i, &density_energy_i, pressure_i, &velocity_i, &fc_i_momentum_x, &fc_i_momentum_y, &fc_i_momentum_z, &fc_i_density_energy);

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
  }
}

static void compute_step_factor(int nelr, float *variables, float *areas, float *step_factors) {
#pragma omp target teams distribute parallel for map(present, alloc : variables[0 : nelr * NVAR], areas[0 : nelr], step_factors[0 : nelr])
  for (int i = 0; i < nelr; i++) {

	

	float density = variables[i + VAR_DENSITY*nelr];
	float3 momentum;
	momentum.x = variables[i + (VAR_MOMENTUM+0)*nelr];
	momentum.y = variables[i + (VAR_MOMENTUM+1)*nelr];
	momentum.z = variables[i + (VAR_MOMENTUM+2)*nelr];

	float density_energy = variables[i + VAR_DENSITY_ENERGY*nelr];

	float3 velocity;       compute_velocity(density, momentum, &velocity);
	float speed_sqd      = compute_speed_sqd(velocity);
	float pressure       = compute_pressure(density, density_energy, speed_sqd);
	float speed_of_sound = compute_speed_of_sound(density, pressure);

	// dt = float(0.5) * sqrtf(areas[i]) /  (||v|| + c).... but when we do time stepping, this later would need to be divided by the area, so we just do it all at once
	step_factors[i] = 0.5f / (sqrtf(areas[i]) * (sqrtf(speed_sqd) + speed_of_sound));
  }
}

static void time_step(int j, int nelr, float *old_variables, float *variables, float *step_factors, float *fluxes) {
#pragma omp target teams distribute parallel for map(present, alloc : old_variables[0 : nelr * NVAR], variables[0 : nelr * NVAR], step_factors[0 : nelr], fluxes[0 : nelr * NVAR])
  for (int i = 0; i < nelr; i++) {

	

	float factor = step_factors[i]/(RK+1-j);

	variables[i + VAR_DENSITY*nelr] = old_variables[i + VAR_DENSITY*nelr] + factor*fluxes[i + VAR_DENSITY*nelr];
	variables[i + VAR_DENSITY_ENERGY*nelr] = old_variables[i + VAR_DENSITY_ENERGY*nelr] + factor*fluxes[i + VAR_DENSITY_ENERGY*nelr];
	variables[i + (VAR_MOMENTUM+0)*nelr] = old_variables[i + (VAR_MOMENTUM+0)*nelr] + factor*fluxes[i + (VAR_MOMENTUM+0)*nelr];
	variables[i + (VAR_MOMENTUM+1)*nelr] = old_variables[i + (VAR_MOMENTUM+1)*nelr] + factor*fluxes[i + (VAR_MOMENTUM+1)*nelr];
	variables[i + (VAR_MOMENTUM+2)*nelr] = old_variables[i + (VAR_MOMENTUM+2)*nelr] + factor*fluxes[i + (VAR_MOMENTUM+2)*nelr];
  }
}

static void initialize_variables(int nelr, float *variables, float *ff_variable) {
#pragma omp target teams distribute parallel for map(present, alloc : variables[0 : nelr * NVAR], ff_variable[0 : NVAR])
  for (int i = 0; i < nelr; i++) {

	
	for(int j = 0; j < NVAR; j++)
		variables[i + j*nelr] = ff_variable[j];
  }
}

static void copy_dd(float *dst, float *src, int N) {  // clEnqueueCopyBuffer device->device
#pragma omp target teams distribute parallel for map(present, alloc : dst[0 : N], src[0 : N])
  for (int i = 0; i < N; i++) dst[i] = src[i];
}

static void dump(const float *h_variables, int nel, int nelr) {
  {
    std::ofstream file("density.txt");
    file << nel << " " << nelr << std::endl;
    for (int i = 0; i < nel; i++) file << h_variables[i + VAR_DENSITY * nelr] << std::endl;
  }
  {
    std::ofstream file("momentum.txt");
    file << nel << " " << nelr << std::endl;
    for (int i = 0; i < nel; i++) {
      for (int j = 0; j != NDIM; j++) file << h_variables[i + (VAR_MOMENTUM + j) * nelr] << " ";
      file << std::endl;
    }
  }
  {
    std::ofstream file("density_energy.txt");
    file << nel << " " << nelr << std::endl;
    for (int i = 0; i < nel; i++) file << h_variables[i + VAR_DENSITY_ENERGY * nelr] << std::endl;
  }
}

inline void compute_flux_contribution_host(float *density, float3 *momentum, float *density_energy, float pressure,
                                           float3 *velocity, float3 *fc_momentum_x, float3 *fc_momentum_y,
                                           float3 *fc_momentum_z, float3 *fc_density_energy) {
  compute_flux_contribution(density, momentum, density_energy, pressure, velocity, fc_momentum_x, fc_momentum_y,
                            fc_momentum_z, fc_density_energy);
}

/* One full run of the solver on already-resident device data (or, with serial=true, on host data). */
struct State {
  int nelr;
  float *areas, *normals, *variables, *old_variables, *fluxes, *step_factors, *fc_mx, *fc_my, *fc_mz, *fc_de, *ffv;
  int *ese;
  float3 *ff_mx, *ff_my, *ff_mz, *ff_de;
};

static void solve(State &s) {
  for (int i = 0; i < iterations; i++) {
    copy_dd(s.old_variables, s.variables, s.nelr * NVAR);
    compute_step_factor(s.nelr, s.variables, s.areas, s.step_factors);
    for (int j = 0; j < RK; j++) {
      compute_flux_contributions(s.nelr, s.variables, s.fc_mx, s.fc_my, s.fc_mz, s.fc_de);
      compute_flux(s.nelr, s.ese, s.normals, s.variables, s.fc_mx, s.fc_my, s.fc_mz, s.fc_de, s.fluxes, s.ffv, s.ff_mx,
                   s.ff_my, s.ff_mz, s.ff_de);
      time_step(j, s.nelr, s.old_variables, s.variables, s.step_factors, s.fluxes);
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 2) { printf("Usage ./cfd <data input file>\n"); return 0; }
  const char *data_file_name = argv[1];

  // far field conditions
  float h_ff_variable[NVAR];
  float3 ff_velocity, h_ff_momentum, h_ff_fc_momentum_x, h_ff_fc_momentum_y, h_ff_fc_momentum_z, h_ff_fc_density_energy;
  {
    const float angle_of_attack = (float)(3.1415926535897931 / 180.0) * (float)(deg_angle_of_attack);
    h_ff_variable[VAR_DENSITY] = (float)(1.4);
    float ff_pressure = (float)(1.0);
    float ff_speed_of_sound = sqrt(GAMMA * ff_pressure / h_ff_variable[VAR_DENSITY]);
    float ff_speed = (float)(ff_mach)*ff_speed_of_sound;
    ff_velocity.x = ff_speed * (float)(cos((float)angle_of_attack));
    ff_velocity.y = ff_speed * (float)(sin((float)angle_of_attack));
    ff_velocity.z = 0.0;
    h_ff_variable[VAR_MOMENTUM + 0] = h_ff_variable[VAR_DENSITY] * ff_velocity.x;
    h_ff_variable[VAR_MOMENTUM + 1] = h_ff_variable[VAR_DENSITY] * ff_velocity.y;
    h_ff_variable[VAR_MOMENTUM + 2] = h_ff_variable[VAR_DENSITY] * ff_velocity.z;
    h_ff_variable[VAR_DENSITY_ENERGY] = h_ff_variable[VAR_DENSITY] * ((float)(0.5) * (ff_speed * ff_speed)) + (ff_pressure / (float)(GAMMA - 1.0));
    h_ff_momentum.x = *(h_ff_variable + VAR_MOMENTUM + 0);
    h_ff_momentum.y = *(h_ff_variable + VAR_MOMENTUM + 1);
    h_ff_momentum.z = *(h_ff_variable + VAR_MOMENTUM + 2);
    compute_flux_contribution_host(&h_ff_variable[VAR_DENSITY], &h_ff_momentum, &h_ff_variable[VAR_DENSITY_ENERGY], ff_pressure,
                                   &ff_velocity, &h_ff_fc_momentum_x, &h_ff_fc_momentum_y, &h_ff_fc_momentum_z, &h_ff_fc_density_energy);
  }

  // read in domain geometry
  int nel, nelr;
  float *h_areas, *h_normals;
  int *h_ese;
  {
    std::ifstream file(data_file_name);
    if (!file) { fprintf(stderr, "cannot open %s\n", data_file_name); return 1; }
    file >> nel;
    nelr = block_length * ((nel / block_length) + std::min(1, nel % block_length));
    h_areas = (float *)calloc(nelr, sizeof(float));
    h_ese = (int *)calloc((size_t)nelr * NNB, sizeof(int));
    h_normals = (float *)calloc((size_t)nelr * NDIM * NNB, sizeof(float));
    for (int i = 0; i < nel; i++) {
      file >> h_areas[i];
      for (int j = 0; j < NNB; j++) {
        file >> h_ese[i + j * nelr];
        if (h_ese[i + j * nelr] < 0) h_ese[i + j * nelr] = -1;
        h_ese[i + j * nelr]--;  // it's coming in with Fortran numbering
        for (int k = 0; k < NDIM; k++) {
          file >> h_normals[i + (j + k * NNB) * nelr];
          h_normals[i + (j + k * NNB) * nelr] = -h_normals[i + (j + k * NNB) * nelr];
        }
      }
    }
    // fill in remaining data: duplicate the last element
    int last = nel - 1;
    for (int i = nel; i < nelr; i++) {
      h_areas[i] = h_areas[last];
      for (int j = 0; j < NNB; j++) h_ese[i + j * nelr] = h_ese[last + j * nelr];
    }
  }

  State s;
  s.nelr = nelr;
  size_t NV = (size_t)nelr * NVAR, N3 = (size_t)nelr * NDIM;
  s.areas = h_areas; s.ese = h_ese; s.normals = h_normals;
  s.variables = (float *)malloc(NV * sizeof(float)); s.old_variables = (float *)malloc(NV * sizeof(float));
  s.fluxes = (float *)malloc(NV * sizeof(float)); s.step_factors = (float *)malloc(nelr * sizeof(float));
  s.fc_mx = (float *)malloc(N3 * sizeof(float)); s.fc_my = (float *)malloc(N3 * sizeof(float));
  s.fc_mz = (float *)malloc(N3 * sizeof(float)); s.fc_de = (float *)malloc(N3 * sizeof(float));
  s.ffv = h_ff_variable; s.ff_mx = &h_ff_fc_momentum_x; s.ff_my = &h_ff_fc_momentum_y; s.ff_mz = &h_ff_fc_momentum_z;
  s.ff_de = &h_ff_fc_density_energy;
  float *result = (float *)malloc(NV * sizeof(float));

  dwarf_init("cfd", "unstructured-grids", "omp");
  dwarf_set_size(nel, nelr);
  std::cout << "Starting..." << std::endl;
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    memset(s.step_factors, 0, nelr * sizeof(float));  // untimed host-side zero source for step_factors
    dwarf_region_begin();
    float *areas = s.areas, *normals = s.normals, *variables = s.variables, *old_variables = s.old_variables,
          *fluxes = s.fluxes, *step_factors = s.step_factors, *fc_mx = s.fc_mx, *fc_my = s.fc_my, *fc_mz = s.fc_mz,
          *fc_de = s.fc_de, *ffv = s.ffv;
    int *ese = s.ese;
    float3 *ff_mx = s.ff_mx, *ff_my = s.ff_my, *ff_mz = s.ff_mz, *ff_de = s.ff_de;
#pragma omp target enter data map(to : areas[0 : nelr], ese[0 : (size_t)nelr * NNB], normals[0 : N3 * NNB], ffv[0 : NVAR], ff_mx[0 : 1], ff_my[0 : 1], ff_mz[0 : 1], ff_de[0 : 1], step_factors[0 : nelr]) \
    map(alloc : variables[0 : NV], old_variables[0 : NV], fluxes[0 : NV], fc_mx[0 : N3], fc_my[0 : N3], fc_mz[0 : N3], fc_de[0 : N3])
    initialize_variables(nelr, variables, ffv);
    initialize_variables(nelr, old_variables, ffv);
    initialize_variables(nelr, fluxes, ffv);
    solve(s);
#pragma omp target exit data map(from : variables[0 : NV]) \
    map(delete : areas[0 : nelr], ese[0 : (size_t)nelr * NNB], normals[0 : N3 * NNB], ffv[0 : NVAR], ff_mx[0 : 1], ff_my[0 : 1], ff_mz[0 : 1], ff_de[0 : 1], step_factors[0 : nelr], old_variables[0 : NV], fluxes[0 : NV], fc_mx[0 : N3], fc_my[0 : N3], fc_mz[0 : N3], fc_de[0 : N3])
    dwarf_region_end(rep);
  }
  std::cout << "Finished" << std::endl;
  std::cout << "Saving solution..." << std::endl;
  dump(s.variables, nel, nelr);
  std::cout << "Saved solution..." << std::endl;
  memcpy(result, s.variables, NV * sizeof(float));

  /* validation: the solution must be finite with positive density; fidelity to the OpenCL app is
   * checked by diffing (density.txt / momentum.txt / density_energy.txt vs the OpenCL run) */
  long bad = 0;
  for (int i = 0; i < nel; i++) {
    float d = result[i + VAR_DENSITY * nelr];
    if (!(std::isfinite(d) && d > 0.0f)) bad++;
  }
  printf("%s: %ld cells with non-finite or non-positive density\n", bad ? "FAILURE" : "SUCCESS", bad);
  return dwarf_report(bad == 0);
}
