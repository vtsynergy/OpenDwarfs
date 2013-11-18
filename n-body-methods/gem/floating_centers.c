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
#include "defines.h"
#include "structures.h"
#include <stdio.h>
#include <math.h>
/*
 * Ok, this one is hard to explain... but I will do my best.
 *      By using A alone, we estimate the center of the molecule (for any
 *  given "center") as some position A distance along the normal toward the
 *  interior of the molecule.  The A alone method is technically the
 * "floating center".  It can result in wrong results for outcroppings
 *  (by assuming the center is somewhere not within the molecule) and it
 * technically smooths over small bumps and pockets. This algorithm makes
 * an attempt to estimate the geometric center within some radius of the
 * point for which we are calculating EP.  This method generally
 * outperforms the standard method but fails miserably on molecules for
 * which the geometric center is very near the surface somewhere (B_DNA).
 *                                                     --jcg
 */


#if USE_FLOATING_CENTERS
/***************************************************************************
 * FUNCTION: floating_center (well not really)                             *
 *                                                                         *
 * INPUTS:   A (estimated electrostatic radius)                            *
 *           center (location we are looking at)                           *
 *           residue (residues making up the molecule)                     *
 *           nres    (number of residues in the molecule)                  *
 *                                                                         *
 * RETURNS:  The estimated center of the molecule from this point in space *
 ***************************************************************************/
vertx floating_center (double A, vertx center, residue *residues, int nres)
{
	/* local variables */
	int i, j, navg;
	double dist, avg_x=0, avg_y=0, avg_z=0;
	double xdiff, ydiff, zdiff;
	vertx to_return;

	double comp_dist = max_R_mult2*A*A;

	for (navg = i = 0; i < nres; i++)
	{

		/* for each residue, find out if the bounding cube falls
		   within the radius of search before searching the residue. */
		xdiff = min(fabs(center.x - residues[i].x[0]),
				fabs(center.x - residues[i].x[1]));

		ydiff = min(fabs(center.x - residues[i].x[0]),
				fabs(center.x - residues[i].x[1]));

		zdiff = min(fabs(center.x - residues[i].x[0]),
				fabs(center.x - residues[i].x[1]));

		dist = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;

		if (dist < comp_dist)
		{
			for (j = 0; j < residues[i].natoms; j++)
			{
				xdiff = residues[i].atoms[j].x - center.x;
				ydiff = residues[i].atoms[j].y - center.y;
				zdiff = residues[i].atoms[j].z - center.z;

				dist = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;

				if (dist <= comp_dist)
				{
					avg_x += residues[i].atoms[j].x;
					avg_y += residues[i].atoms[j].y;
					avg_z += residues[i].atoms[j].z;
					navg ++;
				}
			} /* end for atoms in the residue */
		} /* end if (cube is within radius) */
	}/* end for (natoms) */

	/* average them */
	if (navg > 0)
	{
		avg_x /= (double)navg;
		avg_y /= (double)navg;
		avg_z /= (double)navg;
	}
	else
	{
		avg_x = 0.;
		avg_y = 0.;
		avg_z = 0.;
	}

	/* set the vertex to return */
	to_return.x = avg_x;
	to_return.y = avg_y;
	to_return.z = avg_z;

	/* return the answer */
	return to_return;

} /* end function (floating_center) */
#endif
