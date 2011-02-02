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

#include "file_io.h"
#ifndef NO_UI
#include "tellUser.h"
#else
#define tellUser(type, msg) printf("%s %s", type, msg)
#endif
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/****************************************************************************
 * FUNCTION: runMsms  -- runs the program msms                              *
 *                                                                          *
 * INPUTS: molname    -- the molecule name to be run                        *
 *         residues   -- the residues to be input to the msms algorithm     *
 *         nres       -- the number of residues being input                 *
 *         tri_per_A  -- surface resolution of the desired output           *
 *         probe_size -- the radius of the probe used to represent H20      *
 *                                                                          *
 * OUTPUTS: none                                                            *
 *                                                                          *
 * RETURNS: 1 if MSMS_EXECUTABLE_PATH is set, 0 otherwise                   *
 *                                                                          *
 ****************************************************************************/
int runMsms (char *molname, residue *residues, int nres, double tri_per_A, double probe_size)
{
   /* local variables */
   int  size, len;
   char *fname;
   char *run;
   char *msms_exec;
   char density[64];
   char probe_radius[64];

   /* allocates memory */
   msms_exec = getenv("MSMS_EXECUTABLE_PATH");

   if (msms_exec == NULL)
   {
      return 0;
   }

   /* new memory */
   len = strlen(molname);
   fname = (char *)calloc(len + 6, sizeof(char));

   sprintf(fname, "%s.xyzr", molname);
   sprintf(density, "%2.2f", tri_per_A);
   sprintf(probe_radius, "%2.2f", probe_size);

   write_xyzr(fname, residues, nres);

   /* new memory */
   size = strlen(msms_exec) + 56 + 2*len + strlen(density) + strlen(probe_radius);

   run = (char *)calloc(size, sizeof(char));

   sprintf
     (
        run,
        "%s -if %s -of %s -de %s -probe_radius %s -no_area >/dev/null",
        msms_exec,
        fname,
        molname,
        density,
        probe_radius
     );
   printf("running MSMS with command: %s\n", run);
   fflush(stdout);

   if (system(run) < 0)
       return 0;

   /* free up new memory */
   free(fname);
   free(run);

   return 1;
}
