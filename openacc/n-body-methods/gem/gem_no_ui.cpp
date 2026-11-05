// Stripped GEM

#include <string>
#include <cstring>
#include <iostream>
#include <stdlib.h>

#include <stdio.h>
#include "file_io.h"
#include "defines.h"
#include "calculations.h"
#include "partitioned_open.h"
#include "gem.h"
// #ifdef PGI
// #include <omp.h>
// #include <accel.h>
// #endif


using namespace std;


int cpp_main (int argc, char **argv)
{
   // local variables
   partitioned_open_struct open_dat;
   double A(0),
         proj_len(0.0),
         diel_ext(80.0),
         diel_int(1.0),
         salinity(0.1);

   bool   flags[NUM_FLAGS];
   string options[NUM_FILE_TYPES];

   memset(flags, 0, sizeof(flags));
   memset(&open_dat, 0, sizeof(partitioned_open_struct));

   // parse command line
   if ((argc < 4) || (!check_cmd_line(argc, argv, options, flags)))
   {
        cerr<<"Format for command line options is \'gem [molname diel_ext diel_int salinity] <options>\'"<<endl
            <<"molname, diel_ext, diel_int, and salinity are all required together if present"<<endl
            <<"molname is the name of the molecule"<<endl
            <<"diel_ext is the exterior dielectric charge"<<endl
            <<"diel_int is the interior dielectric charge"<<endl
            <<"salinity is the salinity of the environment"<<endl
            <<"and <options> can be any combination of the following:"<<endl
            <<"\t-projLen       <length>"<<endl
            <<"\t-triDens       <triangles>"<<endl
            <<"\t-probeRad      <radius>"<<endl
            <<"\t-fixA          <electrostatic radius>"<<endl
            <<"\t-dumpVertices  <fileName>"<<endl;

       return 1;
   }

   // open_dat members needed
   open_dat.molname = argv[1];
   open_dat.triDens = 3.0;
   open_dat.probeRadius = 1.5;

   /* get all the dielectric info from command line */
   diel_ext = atof(argv[2]);
   diel_int = atof(argv[3]);
   salinity = atof(argv[4]);

   if (flags[TRI_SPEC])
   {
      open_dat.triDens = atof((char *)options[TRI_SPEC_VAL].c_str());
   }

   if (flags[PROBE_SPEC])
   {
      open_dat.probeRadius = atof((char *)options[PROBE_RADIUS_VAL].c_str());
   }

   if (flags[A_SPEC])
   {
      A = atof((char *)options[A_SPEC_VAL].c_str());
   }

   if (flags[PROJ_SPEC])
   {
      proj_len = atof((char *)options[PROJ_LENGTH_VAL].c_str());
   }

   // get the molecule and the surface
   if (!open_pqr_run_msms(&open_dat))
   {
      // prints more detailed errors in the function call
      return 1;
   }

   if (A < 1) // minimum known atom size is ~1.2 Angstroms
   {
      A = estimate_A (open_dat.residues, open_dat.nresidues);
      cout<<"The estimated electrostatic radius (A) of this molecule is: "<<A<<endl;
   }

   // explicitly set i to 0 (though it was memset to 0 anyway)
   open_dat.i = 0;

   // a good candidate for parallelization
   calc_potential_single_step
   (
      open_dat.residues,
      open_dat.nresidues,
      open_dat.vert,
      open_dat.nvert,
      A,
      proj_len,
      diel_int,
      diel_ext,
      salinity,
      2.0, // ion exclusion radius
      TOTAL_POTENTIAL,
      &(open_dat.i), // first vertex to use
      open_dat.nvert // number of vertices to process
   );

   float  x[3],
          y[3],
          z[3],
         vp[3],
         vc[3];

   extract_statistics
      (
         open_dat.residues, 
         open_dat.nresidues, 
         open_dat.vert, 
         open_dat.nvert, 
         x, y, z, 
         vp, 
         vc
      );

   cout<<"Potential minimum is "<<vp[MIN]<<", maximum is "<<vp[MAX]<<endl;

   if (flags[DUMP_VERTICES])
   {
      cout<<"Dumping vertices to "<<options[VERTEX_FILE_NAME]<<endl;
      dump_vertices
          (
             (char *)options[VERTEX_FILE_NAME].c_str(),
             open_dat.nvert,
             open_dat.vert
          );
   }

   return 0;
}
