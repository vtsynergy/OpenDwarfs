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

/*******************************************************************
 * Given the radius and position of each atom, we are extrapolating*
 * bonds as hybridizations of orbitals (such that the distance     *
 * between the two atoms is less than their radii                  *
 * If you know a better way, please god make it so <><><>   --jcg  *
 ******************************************************************/
#include "calculations.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

/********************************************************************
 * FUNCTION: extrapolate_bonds (tries to find bonds between atoms)  *
 *                                                                  *
 * INPUTS:   residues  (residues making up the molecule)            *
 *           nres      (number of residues in molecule)             *
 *                                                                  *
 * CHANGES:  bnds      (bond structures found)                      *
 *           CALLER RESPONSIBLE FOR MEMORY IN BNDS RETURNED         *
 *                                                                  *
 * RETURNS:  the number of bonds in the array pointed to by bnds    *
 ********************************************************************/ 
int extrapolate_bonds (residue *residues, int nres, bond **bnds)
{
   /* local variables */
   bond *bonds;
   int nbonds, est;
   int r1, r2, a1, a2;

   /* performing boxing on atoms to check interresidue bonds */
   double nonbondlen,
          res_rad,
          atom_rad,
          distance,
          res_x,
          res_y,
          res_z;

   double d2;

   /* estimate the number of bonds as 3*natoms to start */
   /* estimate number of atoms per residue to be 20 */
   est = 3*nres*20;
   bonds = (bond *)calloc(est, sizeof(bond));

   if (bonds == NULL)
   {
      fprintf(stderr, "failed to allocate space for %i bonds\n", est);
      return 0;
   }

   /* step through all residues */
   for (nbonds = r1 = 0; r1 < nres; r1++)
   {

      /* for each atom in a residue */
      for (a1 = 0; a1 < residues[r1].natoms; a1++)
      {
         /* build a cube that barely contains the atom */
         atom_rad = residues[r1].atoms[a1].radius;

         /* for each noncomputed residue, including this one */
         for (r2=r1; r2 < nres; r2++)
         {
            /* center of the residue */
            res_x = residues[r2].x;
            res_y = residues[r2].y;
            res_z = residues[r2].z;

            /* estimate the residue radius as that of a circumscribing sphere */
            res_rad = residues[r2].rad;

            distance = sqrt(dist2(
                             residues[r1].atoms[a1].x,
                             residues[r1].atoms[a1].y,
                             residues[r1].atoms[a1].z,
                             res_x, res_y, res_z));

            /* if the atom falls within the residue bounds... */
            if (distance < atom_rad + res_rad)
            {
               /* check for bonds within atoms in that residue */             
               for (a2 = 0; a2 < residues[r2].natoms; a2++)
               {

                  /* skip self-self */
                  if ((a2 == a1) && (r2 == r1)) continue;

                  /* calculate the distance squared */
                  d2 = dist2
                          (
                             residues[r1].atoms[a1].x,
                             residues[r1].atoms[a1].y,
                             residues[r1].atoms[a1].z,
                             residues[r2].atoms[a2].x,
                             residues[r2].atoms[a2].y,
                             residues[r2].atoms[a2].z
                          );
   
                  nonbondlen = max(atom_rad, residues[r2].atoms[a2].radius);

                  if (d2 < nonbondlen * nonbondlen)
                  {
                     /* re-estimate if we overstep */
                     if (nbonds == est)
                     {
                        /* add 3 * the remaining number of atoms */
                        est += 3*(nres - r1)*20;
                        if (nres == 0)
                        {
                           est += 3*residues[r1].natoms;
                        }
      
                        bonds = (bond *)realloc(bonds, est * sizeof(bond));
      
                        if (bonds == 0)
                        {
                           fprintf(stderr, "failed to allocate space for %i bonds\n", est);
                           return 0;
                        }
                     }

                     /* set the bond up */
                     bonds[nbonds].a1 = a1;
                     bonds[nbonds].r1 = r1;
                     bonds[nbonds].a2 = a2;
                     bonds[nbonds].r2 = r2;
                     bonds[nbonds].length = sqrt(d2);
   
                     nbonds++;
                  } /* if the distance squared is less than 1 of the radii squared */
               } /* for all atoms within the residue */
            } /* if the atom overlaps the residue */
         } /* for residues >= r1 */
      } /* for atoms in this bond */
   } /* for nbonds (big loop) */

   /* free up overestimation memory */
   bonds = (bond *) realloc (bonds, nbonds*sizeof(bond));

   /* return what we found */
   *bnds = bonds;

   printf("extrapolate bonds found %i bonds\n", nbonds);
   fflush(stdout);

   return nbonds;
}
