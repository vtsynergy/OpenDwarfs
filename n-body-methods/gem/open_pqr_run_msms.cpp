/****************************************************************************
 * GEM -- electrostatics calculations and visualization                     *
 * Copyright (C) 2006  John C. Gordon                                       *
 *                                                                          *
 * This program is free software; you can redistribute it and/or modify     *
 * it under the terms of the GNU General Public License as published by     *
 * the Free Software Foundation; either version 2 of the License, or        *
 * (at your option) any later version.                                      *
 *                                                                          *
 * This program is distributed in the hope that it will be useful,          *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 * GNU General Public License for more details.                             *
 *                                                                          *
 * You should have received a copy of the GNU General Public License along  *
 * with this program; if not, write to the Free Software Foundation, Inc.,  *
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.              *
 ****************************************************************************/

#include "partitioned_open.h"
#include <iostream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "file_io.h"
#ifndef NO_UI
#include "visualize.h"
#endif
#include "calculations.h"

using namespace std;

/* this function opens a pqr file and runs msms on it then extracts
 * the surface from the .vert and .tri files --jcg */
extern "C" int open_pqr_run_msms(partitioned_open_struct *open_dat)
{
   /* local variables */
   /*******************/

   /* various filenames */
   string tempstr1,
          tempstr2;

   /* sanity check */
   if (open_dat == NULL) return 0;

   /* initialize open_dat to sanely return errors */
   /***********************************************/
   open_dat->residues = NULL;
   open_dat->nresidues = 0;
   open_dat->tri = NULL;
   open_dat->ntri = 0;
   open_dat->vert = NULL;
   open_dat->nvert = 0;

   /* read in atom locations from the pqr file */
   /********************************************/
   tempstr1 = open_dat->molname;
   open_dat->nresidues = read_pqr((char *)tempstr1.c_str(), &(open_dat->residues));

   if (open_dat->nresidues == 0)
   {
      tempstr1 += ".pqr";
      open_dat->nresidues = read_pqr((char *)tempstr1.c_str(), &(open_dat->residues));

      if (open_dat->nresidues == 0)
      {
         cout << "Error, no pqr file named either "
              << open_dat->molname<<" or "<<tempstr1<<endl;
         return 0;
      }
   }
      
   fflush(stderr);
   fflush(stdout);

   /* come up with the surface */
   /****************************/
   tempstr1 = open_dat->molname;

   if (!runMsms
        (
          (char *)tempstr1.c_str(),
          open_dat->residues,
          open_dat->nresidues,
          open_dat->triDens,
          open_dat->probeRadius
        ))
     {
        cout<<"Warning: GEM needs the MSMS_EXECUTABLE_PATH environment"
            <<"variable to be set in order to run msms directly."
            <<"searching for precomputed surfaces..."<<endl;
      
         /* read the surface */
         /********************/
         read_msms_output
            (
               (char *)tempstr1.c_str(),
               &(open_dat->vert),
               &(open_dat->tri),
               &(open_dat->nvert),
               &(open_dat->ntri)
            );

         if ((open_dat->nvert > 0) || (open_dat->vert == NULL))
         {
            cout<<"Found precomputed surfaces, continuing."<<endl;
         }
         else
         {
            cerr<<"ERROR, could not find precomputed surfaces."<<endl;
            return 0;
         }
      }
      else
      {
         read_msms_output
           (
              (char *)tempstr1.c_str(),
              &(open_dat->vert),
              &(open_dat->tri),
              &(open_dat->nvert),
              &(open_dat->ntri)
           );

      }
      
      return (((open_dat->nvert < 1) || (open_dat->vert == NULL))?0:1);
}
