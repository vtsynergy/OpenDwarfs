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
 * FUNCTION: write_xyzr  --writes an xyzr file for msms to read             *
 *                                                                          *
 * INPUTS: fname  -- name of the file to write                              *
 *         residues -- the molecule data to write to file                   *
 *         nres     -- the number of residues in the molecule to write      *
 *                                                                          *
 * OUTPUTS: none                                                            *
 *                                                                          *
 * RETURNS: nothing                                                         *
 *                                                                          *
 ****************************************************************************/
void write_xyzr(char *fname, residue *residues, int nres)
{
   /* local variables */
   int i, j;
   FILE *fp;

   fp = fopen(fname, "w");

   for (i = 0; i < nres; i++)
   {
      for (j = 0; j < residues[i].natoms; j++)
      {
         fprintf(fp, "%f\t%f\t%f\t%f\n",
                     residues[i].atoms[j].x,
                     residues[i].atoms[j].y,
                     residues[i].atoms[j].z,
                     residues[i].atoms[j].radius);
      }
   }

   fclose(fp);

}
