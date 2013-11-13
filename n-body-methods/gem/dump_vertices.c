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
#include <stdio.h>

/****************************************************************************
 * FUNCTION: dump_vertices --very simple function to just dump out ver info *
 *                                                                          *
 * INPUTS:   fname   --file name to write                                   *
 *           nvert   --number of vertices to dump                           *
 *            vert   --vertices to write out                                *
 *                                                                          *
 * OUTPUTS: none                                                            *
 *                                                                          *
 * RETURNS: nothing                                                         *
 ****************************************************************************/
void dump_vertices(char *fname, int nvert, vertx *vert)
{
/* local variables */
int i;
FILE *fp;

   fp = fopen(fname, "w");

   if (fp == NULL)
   {
      fprintf(stderr, "Error opening %s to write\n", fname);
      return;
   }

   fprintf(fp, "X    Y    Z     potential\n");

   for (i = 0; i < nvert; i++)
   {
       fprintf(fp, "%lf %lf %lf    %lf\n", vert[i].x, vert[i].y, vert[i].z, vert[i].potential);
   }

   fclose(fp);
}
