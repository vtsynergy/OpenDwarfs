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
#include "calculations.h"
#include "structures.h"
#include <math.h>
#include <stdio.h>

/****************************************************************************
 * FUNCTION: extract_statistics  -- one loop to extract statistics for all  *
 *                                     data contained within vertices/atoms *
 *                                                                          *
 * INPUTS:  res  -- the residues making up the molecule in question         *
 *          nres -- number of residues in the molecule                      *
 *          vert -- vertices composing the surface of the molecule          *
 *          nvert -- number of vertices on the surface of the molecule      *
 *                                                                          *
 * OUTPUTS: x, y, z -- spacial min/max/extents of the molecule min bounding *
 *                     box                                                  *
 *          vp      -- vertex potential min/max/extents for the molecule    *
 *          vc      -- vertex charge min/max/extents for the molecule       *
 *                                                                          *
 ****************************************************************************/
int extract_statistics(residue *res, int nres, vertx *vert, int nvert, 
                        float x[3], float y[3], float z[3], float vp[3], float vc[3])
{
  /* local variables */
  int i; 
  double potential_sum;

  if ((res == NULL) || (vert == NULL) || (nvert < 3))
  {
     /* have to initialize xyz bounds for scales to not be 0 */
     x[MIN] = y[MIN] = z[MIN] = -1;
     x[MAX] = y[MAX] = z[MAX] =  1;
     return 0;
  }

  /* else */

  /* reinitialize to the zeroth element */
  vp[MIN] = vp[MAX] = vert[0].potential;
  vc[MIN] = vc[MAX] = vert[0].surface_charge;
  x[MIN]  = x[MAX]  = vert[0].x;
  y[MIN]  = y[MAX]  = vert[0].y;
  z[MIN]  = z[MAX]  = vert[0].z;

  potential_sum = vert[0].potential*vert[0].potential;
 
  /* seek through and grab our values */
  for (i = 1; i < nvert; i++)
  {
     potential_sum += vert[i].potential*vert[i].potential;

     /* total potential bounds */
     if (vp[MIN] > vert[i].potential)
        vp[MIN] = vert[i].potential;
     else if (vp[MAX] < vert[i].potential)
        vp[MAX] = vert[i].potential;

     /* total charge bounds */
     if (vc[MIN] > vert[i].surface_charge)
        vc[MIN] = vert[i].surface_charge;
     else if (vc[MAX] < vert[i].surface_charge)
        vc[MAX] = vert[i].surface_charge;
   
     /* check x bounds */
     if (x[MIN] > vert[i].x)
        x[MIN] = vert[i].x;
     else if (x[MAX] < vert[i].x)
        x[MAX] = vert[i].x;
   
     /* check y bounds */
     if (y[MIN] > vert[i].y)
        y[MIN] = vert[i].y;
     else if (y[MAX] < vert[i].y)
        y[MAX] = vert[i].y;
   
     /* check z bounds */
     if (z[MIN] > vert[i].z)
        z[MIN] = vert[i].z;
     else if (z[MAX] < vert[i].z)
        z[MAX] = vert[i].z;
   
   } /* end for (i < nvert) */

   potential_sum /= nvert;

   printf("RMSD of total phi values is %f\n", sqrt(potential_sum));
   
   return 1;
} /* end extract statistics function */
