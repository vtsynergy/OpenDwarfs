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
#include <math.h>

/***************************************************************************
 * FUNCTION:  dot_prod  --returns the dot product of 2 vectors in 3 space  *
 *                                                                         *
 * INPUTS:    v1, v2    -- 3 space vectors {x, y, z}                       *
 *                                                                         *
 * OUTPUTS:   v1 * v2  <-- * is dot                                        *
 *                                                                         *
 ***************************************************************************/
double dot_prod (double v1[3], double v2[3])
{
   return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

/***************************************************************************
 * FUNCTION:  magnitude  -- returns the magnitude of the vector v          *
 *                                                                         *
 * INPUTS:    v   -- 3 space vector {x, y, z}                              *
 *                                                                         *
 * OUTPUTS:   |v|                                                          *
 *                                                                         *
 ***************************************************************************/
double magnitude (double v[3])
{
    return dist(0.,0.,0., v[0], v[1], v[2]);
}

/***************************************************************************
 * FUNCTION:  cross_prod  -- returns the cross product of v1 and v2        *
 *                                                                         *
 * INPUTS:    v1, v2   -- 3 space vectors {x, y, z}                        *
 *                                                                         *
 * OUTPUTS:   cross = v1 x v2                                              *
 *                                                                         *
 ***************************************************************************/
void cross_prod (double v1[3], double v2[3], double cross[3])
{

   /* local variables */
   const int x = 0,
             y = 1,
             z = 2;

   cross[x] = v1[y]*v2[z] - v1[z]*v2[y];
   cross[y] = v1[z]*v2[x] - v1[x]*v2[z];
   cross[z] = v1[x]*v2[y] - v1[y]*v2[x];
}

