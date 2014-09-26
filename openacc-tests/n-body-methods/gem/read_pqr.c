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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef NO_UI
#include "tellUser.h"
#else
#define tellUser(type, msg) printf("%s %s", type, msg)
#endif
#include "structures.h"
#include "calculations.h"
#include <math.h>
#include "defines.h"

#define PQR_LINE_SIZE 512
#define PQR_CHARGE_WARN_TOL .01

/**************************************************************************
 * FUNCTION: gather_mol_info --gathers information about spatial          *
 *                             requirements to store the molecule         *
 *                                                                        *
 * INPUTS:  fp  -- the file pointer (assumedly at beginning of file)      *
 *                                                                        *
 * OUTPUTS: res -- pointer to residues to be allocated and filled         *
 *                 such that the number of atoms are known                *
 *                 CALLER IS RESPONSIBLE FOR MEMORY                       *
 *                                                                        *
 * RETURNS: number of residues read in                                    *
 *                                                                        *
 **************************************************************************/
int gather_mol_info (residue **res, FILE *fp)
{
   /* local variables */
   int nres      = 0;
   residue *residues = NULL;
   int natoms  = 0;
   char last_res[5];
   int  last_res_num = -99999;
   char res_name[5];
   int  res_num, i, total_atoms, linecount;
   char    *fgetsret;
   char line[PQR_LINE_SIZE];

   /* now we actually get to read in the data */

   /* first find out how many residues there are */
   fgetsret = (char *)1;

   total_atoms = natoms = linecount = 0;
   for (i = 0; fgetsret && fp; i++)
   {
      /* read in the line as a string */
      fgetsret = fgets(line, PQR_LINE_SIZE, fp);

      linecount++;

      if ((strstr(line, "ATOM") == NULL) && (strstr(line, "HETATM") == NULL))
         continue;

      /* check to see if there is a chain id here, if so, ignore it */
      if (line[21] != ' ')
      {
            /* all we care about is the residue name and number */
            sscanf
            (
               line,
               "%*s %*d %*s %s %*s %d %*f %*f %*f %*f %*f",
               res_name,
               &res_num
            );
      }
      else
      {
            /* all we care about is the residue name and number */
            sscanf
            (
               line,
               "%*s %*d %*s %s %d %*f %*f %*f %*f %*f",
               res_name,
               &res_num
            );
      }

      if ((res_num != last_res_num) || (strcmp(res_name, last_res) != 0))
      {
         strcpy(last_res, res_name);
         last_res_num = res_num;
         nres++;
      }


   }/* end for i < natoms && fp && fgetsret */

   /* rewind the file */
   fseek(fp, 0, SEEK_SET);

   printf("This file is %i lines long.\n", linecount);

   *res = residues = (residue *)calloc(nres, sizeof(residue));

   last_res_num = -9999;

   /* now we actually get to read in the data about the data */
   fgetsret = (char *)1;
   for (natoms = i = 0;fgetsret && fp && i <= nres;) /* read all atoms in the last res.. */
   {
      /* read in the line as a string */
      fgetsret = fgets(line, PQR_LINE_SIZE, fp);

      if ((!fgetsret)||((strstr(line, "ATOM") == NULL) && (strstr(line, "HETATM") == NULL)))
         continue;

      if (line[21] != ' ')
      {
         /* all we care about is the last 5 fields */
         sscanf
            (
               line,
               "%*s %*d %*s %s %*s %d %*f %*f %*f %*f %*f",
               res_name,
               &res_num
            );
      }
      else
      {
         /* all we care about is the last 5 fields */
         sscanf
            (
               line,
               "%*s %*d %*s %s %d %*f %*f %*f %*f %*f",
               res_name,
               &res_num
            );
      }

      if ((res_num != last_res_num) || (strcmp(res_name, last_res) != 0))
      {
         if (i > 0)
         {
           strcpy(residues[i-1].name, last_res);
           residues[i-1].res_num = last_res_num;
           residues[i-1].natoms = natoms;
           if (natoms != 0) /* unthinkable, but... */
            residues[i-1].atoms = (atom *)calloc(natoms, sizeof(atom));

         }
         strcpy(last_res, res_name);
         last_res_num = res_num;

         total_atoms += natoms;
         natoms = 0;

         i++;
      }

      natoms++;

   } /* end for i < nres && fp && fgetsret */

   if (i > 0)
   {
      strcpy(residues[i-1].name, last_res);
      residues[i-1].res_num = last_res_num;

      residues[i-1].natoms = natoms;
      residues[i-1].atoms = (atom *)calloc(natoms, sizeof(atom));
      strcpy(residues[i-1].name, last_res);
   }
   total_atoms += natoms;

   printf("This molecule contains %i residues and %i atoms.\n", nres, total_atoms);

   /* rewind the file (again) */
   fseek(fp, 0, SEEK_SET);

   return nres;
}

/**************************************************************************
 * FUNCTION: read_pqr --reads in a pqr file and populates info about the  *
 *                                                       molecule         *
 *                                                                        *
 * INPUTS:  fname  -- the file name to read                               *
 *                                                                        *
 * OUTPUTS: residues -- pointer to residues to be allocated and filled    *
 *                 CALLER IS RESPONSIBLE FOR MEMORY                       *
 *                                                                        *
 * RETURNS: number of residues read in                                    *
 *                                                                        *
 **************************************************************************/
int read_pqr(char *fname, residue **residues)
{
   /* local variables */
   FILE    *fp;
   int     nres=0,
           i,
           j;
   float   mol_charge = 0,
           temp;
   char    line[PQR_LINE_SIZE];
   residue *res;
   char    *fgetsret;
   double  x[2], y[2], z[2];
   double  xvec, yvec, zvec;
   double  spacial_bound;


   *residues = NULL;

   fp = fopen(fname, "r");
   if (!fp) 
      return 0;


   /* gather the residue profile and allocate space */
   nres = gather_mol_info(residues, fp);
   res = *residues;

   /* now we actually get to read in the data */
   fgetsret = (char *)1;
   for (i = 0;fgetsret && fp && (i < nres);i++)
   {

      /* initialize bounding box */
      x[0] = 9999.;
      x[1] = -9999.;
      y[0] = 9999.;
      y[1] = -9999.;
      z[0] = 9999.;
      z[1] = -9999.;

      for (j = 0; (j < res[i].natoms); j++)
      {
         /* read in the line as a string */
         fgetsret = fgets(line, PQR_LINE_SIZE, fp);
         if ((strstr(line, "ATOM") == NULL) && (strstr(line, "HETATM") == NULL))
         {
            j--;
            continue;
         }

         if (line[21] != ' ')
         {
            /* all we care about is the last 5 fields */
            sscanf
               (
                  line,
                  "%*s %d %s %s %*s %d %f %f %f %f %f",
                  &res[i].atoms[j].index,
                  res[i].atoms[j].name,
                  res[i].name,
                  &res[i].res_num,
                  &res[i].atoms[j].x,
                  &res[i].atoms[j].y,
                  &res[i].atoms[j].z,
                  &res[i].atoms[j].charge,
                  &res[i].atoms[j].radius
               );
         }
         else
         {
            /* all we care about is the last 5 fields */
            sscanf
               (
                  line,
                  "%*s %d %s %s %d %f %f %f %f %f",
                  &res[i].atoms[j].index,
                  res[i].atoms[j].name,
                  res[i].name,
                  &res[i].res_num,
                  &res[i].atoms[j].x,
                  &res[i].atoms[j].y,
                  &res[i].atoms[j].z,
                  &res[i].atoms[j].charge,
                  &res[i].atoms[j].radius
               );
         }

         /* check the bounds of the residue cube */
         /****************************************/
             res[i].x += res[i].atoms[j].x;
             res[i].y += res[i].atoms[j].y;
             res[i].z += res[i].atoms[j].z;

             /* X min */
             spacial_bound = res[i].atoms[j].x - res[i].atoms[j].radius;
             if (spacial_bound < x[0])
                x[0] = spacial_bound;
             /* X max */
             spacial_bound = res[i].atoms[j].x + res[i].atoms[j].radius;
             if (spacial_bound > x[1])
                x[1] = spacial_bound;

             /* Y min */
             spacial_bound = res[i].atoms[j].y - res[i].atoms[j].radius;
             if (spacial_bound < y[0])
                y[0] = spacial_bound;
             /* Y max */
             spacial_bound = res[i].atoms[j].y + res[i].atoms[j].radius;
             if (spacial_bound > y[1])
                y[1] = spacial_bound;

             /* Z min */
             spacial_bound = res[i].atoms[j].z - res[i].atoms[j].radius;
             if (spacial_bound < z[0])
                z[0] = spacial_bound;
             /* Z max */
             spacial_bound = res[i].atoms[j].z + res[i].atoms[j].radius;
             if (spacial_bound > z[1])
                z[1] = spacial_bound;

         /********************/
         /* Done with bounds */
         res[i].atoms[j].charge *= ELECTROSTATIC_CONVERSION_FACTOR;

         /* keep a running total of the charges on the residue */
         res[i].total_charge += res[i].atoms[j].charge;

      } /* end for (j < natoms) */

      res[i].x /= (float) res[i].natoms;
      res[i].y /= (float) res[i].natoms;
      res[i].z /= (float) res[i].natoms;

      xvec = max((res[i].x - x[0]), (x[1] - res[i].x));
      yvec = max((res[i].y - y[0]), (y[1] - res[i].y));
      zvec = max((res[i].z - z[0]), (z[1] - res[i].z));

      res[i].rad = sqrt(xvec*xvec + yvec*yvec + zvec*zvec);
      
      mol_charge += res[i].total_charge;

      temp = res[i].total_charge / ELECTROSTATIC_CONVERSION_FACTOR;

#ifdef __ANNOYING_TELLUSER_THAT_SOMETIMES_DUMPS__
      /* warn if residue has non-integer charge */
      if ((temp - (int)temp) > PQR_CHARGE_WARN_TOL)
      {
         sprintf(line, "Residue %s %i: total charge is non-integer (%f)\n",
                       res[i].name, res[i].res_num, temp);
         tellUser("Warning", line);
      }
#endif

   }/* end for i < nres && fp && fgetsret */

   printf("Total molecular charge: %2.2f\n", mol_charge/ELECTROSTATIC_CONVERSION_FACTOR);

   fflush(stdout);

   return nres;
} /* end read_pqr */
