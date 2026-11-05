//////////////////////////////////////////////////////////////////////////////
// 
// esize.cc Calculates the electrostatic size, A, of a molecule using 
// atomic coordinates and radii  from its PQR (or simply xyz) file.
// Author: Grigori Sigalov <sigalov@vt.edu>, <greg_sigalov@yahoo.com>
//
// Algorithm Details are in: 
//
// "Analytical Linearized Poisson--Boltzmann Approach 
//   for Molecular Dynamics Applications", by Sigalov, Fenley and Onufriev
// submitted to  J. Chem. Phys., 2006
//
////////////////////////////////////////////////////////////////////////////


 
// Brief algorithm: First, the moments of inertia of the molecule are
// calculated. If the program is run without second parameter or, which is the
// same, with second parameter being "-det", the electrostatic size is found
// from the third invariant (determinant) of the tensor of inertia.  Otherwise
// some more work has to be done. The _principal_ moments of inertia are found
// by solving the cubic equation det(I - lambda E) = 0. The semiaxes of the
// ellipsoid that has same principal moments of inertia and mass as the
// original molecule are found. Then, there are two possibilities: option
// "-ell" leads to calculation of the exact electrostatic size, A_ell, of the
// effective ellipoid just found, which involves numerical summation by
// Simpson method of an elliptic intergal of the first kind (to handle the
// infinite limit, the integral is transformed to map infinity to zero). With
// option "-elf", the value A_elf, which is an approximation to A_ell
// expressed using elementary functions only, is calculated. Option "-abc"
// prints the semiaxes of the effective ellipsoid.  Option "-tab" prints all
// of the above into a table; option -hea prints table with a header.
// Finally, option -deb prints additionally some extra information. Any of
// these options, if used at all, MUST follow the input file name. 

// By default, the input file is in PQR format. If you have only PDB file
// available, convert it to a plain XYZ format and use an extra option,
// "-xyz". In this case all atomic sizes are set to the same value,
// DEFAULTATOMSIZE.

// Shortly, the following hierarchy can be built:

// 1. -det (default) gives a reasonable estimate of the electrostatic size. In
// most cases, you will want nothing more. Moreover, A_det has simple
// analytical derivatives, so if you are going to use derivatives of the
// electrostatic size, you should choose A_det. This option involves
// calculation of the moments of inertia of the molecule, which is a very
// straighforward procedure.  

// 2. -elf adds some accuracy if the molecule's shape is close to an
// ellipsoid, which rarely happens, at the price of a cubic equation to solve.

// 3. -ell gives an exact solution for an ellipsoidal molecule, and involves
// numerical integration in addition to the cubic equation. Options -abc,
// -tab, -hea, and -deb include -ell and therefore take the same amount of
// calculation. 


#include <stdio.h>  
#include <math.h>   
#include <string>   
#include <iostream>
#include "structures.h"

#define DEFAULTATOMSIZE 1.5
#define ACCURACY 1.E-09 // convergence of Simpson's method, stopping condition  

#define TOL 1.E-12		// values < TOL are treated as zeros

using namespace std; 

int debug = 0;	// this setting (default) suppresses all warnings; call with key -deb if unsure about your input file

extern "C" double estimate_A (residue *residues, int nresidues)
{

   double  x, y, z, r2, r3; /* r^3 (an estimate of atomic mass) */

	long long i, j;

	double xav = 0., yav = 0., zav = 0., // center-of-mass coordinates
		molecule_mass = 0., atom_mass;

    for (i = 0; i < nresidues; i++)
    {
       for (j = 0; j<residues[i].natoms; j++)
       {
         r3 = residues[i].atoms[j].radius;
         r3 *= r3*r3;

		   xav += residues[i].atoms[j].x * r3; 	// atom's "mass" is simply r^3
		   yav += residues[i].atoms[j].y * r3; 
		   zav += residues[i].atoms[j].z * r3;
		   molecule_mass += r3;
	    }
    }
	
	xav /= molecule_mass;	// now it's the center
	yav /= molecule_mass;
	zav /= molecule_mass;
	
	double x2, y2, z2, atoms_MI;
	double I11 = 0., I12 = 0., I13 = 0., I22 = 0., I23 = 0., I33 = 0.;

   for (i = 0; i < nresidues; i++)
   {
	   for (j=0; j<residues[i].natoms; j++)	// calculating the moments of inertia
	   {
         r2 = residues[i].atoms[j].radius;
         r2 *= r2;
         r3 = r2 * residues[i].atoms[j].radius;

		   atom_mass = r3;

		   atoms_MI = r2 / 5.;	// half atom's moment of ineria about diameter (will be added twice!)

         x = residues[i].atoms[j].x - xav;
         y = residues[i].atoms[j].y - yav;
         z = residues[i].atoms[j].z - zav;
		   
		   x2 = x * x + atoms_MI;
		   y2 = y * y + atoms_MI;
		   z2 = z * z + atoms_MI;
		   
		   I11 += atom_mass * (y2 + z2);
		   I22 += atom_mass * (z2 + x2);
		   I33 += atom_mass * (x2 + y2);
		   
		   I12 -= atom_mass * x * y;		// atoms' moments do not contribute 
		   I13 -= atom_mass * x * z;		// because of symmetry
		   I23 -= atom_mass * y * z;
	   }
   }
	
	if (debug)
	{
		printf("Original tensor of inertia: %12.4g %12.4g %12.4g\n", I11, I12, I13);
		printf("                            %12.4g %12.4g %12.4g\n", I12, I22, I23);
		printf("                            %12.4g %12.4g %12.4g\n", I13, I23, I33);
	}
	
   // a determinant code
	double det_I = I11 * I22 * I33 + 2. * I12 * I23 * I13 - I11 * I23 * I23  - I22 * I13 * I13 - I33 * I12 * I12;
	if (det_I <= 0.)
	{
	 	cout << "Problem with the determinant of the tensor of intertia: " << endl;
		printf("Det I = %g\n", det_I);
	}
	return sqrt( 2.5 / molecule_mass ) * pow( det_I, 1./6. );
}
