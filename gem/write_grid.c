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
/**************************************************************************
 * This file contains a function which writes out a grid in DelPhi format *
 * for use with GRASP/DelPhi related visualization routines. --jcg        *
 *************************************************************************/
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "defines.h"
#include "file_io.h"

/****************************************************************************
 * FUNCTION:   write_grid -- generic write grid, determines grid type from  *
 *                           filename or defaults to delphi                 *
 *                                                                          *
 * INPUTS:     filename   -- name of the file                               *
 *             dim        -- dimension of the cubic map on one edge         *
 *             grid       -- grid data to output                            *
 *             cen_x, cen_y, cen_z -- center of the grid in space           *
 *             wid                 -- dimension of the grid in space        *
 *                                                                          *
 * RETURNS:    -1 == bad arguments                                          *
 *              0 == problems with file                                     *
 *              1 == success                                                *
 ***************************************************************************/
int write_grid (char *filename, int dim, float *grid, float cen_x, float cen_y, float cen_z, float wid)
 {
 /* local variables */
 /*******************/
 const int DELPHI_FORMAT = 0;
 const int MEAD_FORMAT   = 1;
 const int BAD_ARGS      = -1;
 const int FILE_PROBLEMS = 0;
 const int SUCCESS       = 1;
 int format, length, ret_val;
 char *finalname=NULL;

 float min_x, min_y, min_z,
       max_x, max_y, max_z;

    if (filename == NULL)
    {
       fprintf(stderr, "error in write_grid, null filename\n");
       return BAD_ARGS;
    }

    /* else */
    length = strlen(filename);

    if ((length < 5) || (filename[length-1] == '/'))
    {
      fprintf(stderr, "Error in write_grid, grid is not a discernable file type or it is a directory\n");
      return BAD_ARGS;
    }

    /* else */
    if (filename[length-4] == '.')
    {
       if (strcmp(&filename[length-4], ".fld") == 0)
       {
          finalname = (char *)malloc(length*sizeof(char));
          strcpy(finalname, filename);

          format = MEAD_FORMAT;
       }
       else if (strcmp(&filename[length-4], ".phi") == 0)
       {
          finalname = (char *)malloc(length*sizeof(char));
          strcpy(finalname, filename);

          format = DELPHI_FORMAT;
       }
       else
       {
          finalname = (char *)malloc((length+4) * sizeof(char));
          sprintf(finalname, "%s.phi", filename);
   
          format = DELPHI_FORMAT;
       }
    }
    else
    {
       finalname = (char *)malloc((length+4) * sizeof(char));
       sprintf(finalname, "%s.phi", filename);

       format = DELPHI_FORMAT;
    }

    printf("Saving phimap file >%s< in %s format\n",
            finalname, ((format == DELPHI_FORMAT)? "DelPhi" : "MEAD"));

    if (format == DELPHI_FORMAT)
    {
       ret_val = write_phi_grid
                  (
                     finalname,
                     dim,
                     grid,
                     (wid)/(float)(dim-1),
                     cen_x,
                     cen_y,
                     cen_z
                  );

       if (!ret_val) return FILE_PROBLEMS;
    }
    else
    {
       min_x = cen_x - (wid/2.);
       min_y = cen_y - (wid/2.);
       min_z = cen_z - (wid/2.);
       max_x = cen_x + (wid/2.);
       max_y = cen_y + (wid/2.);
       max_z = cen_z + (wid/2.);

       ret_val = write_avs
                 (
                    finalname,
                    grid,
                    dim,
                    min_x, min_y, min_z,
                    max_x, max_y, max_z
                 );

       if (!ret_val) return FILE_PROBLEMS;
    }

    return SUCCESS;
 }




/****************************************************************************
 * According to documentation, the PhiMap .grd file is formatted as follows *
 * character*20 uplbl                                                       *
 * character*10 nxtlbl, character*60 toplbl                                 *
 * real*4 phi(n,n,n) {n=65 typically}                                       *
 * character*16 botlbl                                                      *
 * real*4 scale, oldmid(3)                                                  *
 *                                                                          *
 * scale is the reciprocal grid spacing, mid is the midpoint,               *
 * uplbl, nxtlbl, and toplbl are all just character strings,                *
 * and phi is the potential map data.                                       *
 *                                                            --jcg         *
 ***************************************************************************/

/* one thing about fortran unformatted output is that it appends and
 * concatenates a header to all output, so when emulating fortran output
 * we have to output headers for each line....
 * so writing x is now sizeof(x) x sizeof(x)
 */

/* NOTE:  This routine assumes internal representation of potential (kCAL_PER_MOL) --
 * so it converts back to kT, if you are writing kT data then you need to convert
 * to internal representation first unfortunately --jcg
 */
int write_phi_grid (char *fname, int szg, float *gridData, float scale, float midx, float midy, float midz)
{
   /* local variables */
   /*******************/
   FILE *fp;
   int sz, i;
   int szgcub = szg*szg*szg; /* grid size squared gets reused a lot */

   /* try to open the file */
   fp = fopen (fname, "w");   /* write only */

   /* check for errors */
   if (fp == NULL)
   {
      fprintf(stderr, "Error writing grid file\n");
      return 0;
   }

   /* write uplbl */
   /***************/
   sz = 20;
   /* always starts with "now starting phimap " */

   fwrite(&sz, 4, 1, fp); /* header */
   fwrite("now starting phimap ", 1, 20, fp); /* fixed */
   fwrite(&sz, 4, 1, fp); /* footer */

   /* write nextlbl, toplbl */
   /*************************/
   sz=70;
   fwrite(&sz, 4, 1, fp); /* header */
   fwrite("01234567890123456789012345678901234567890123456789012234567890123456789",sizeof(char),70,fp);
   fwrite(&sz, 4, 1, fp); /* footer */

   /* write the grid */
   /******************/
   sz = szgcub;

   for (i = 0; i < sz; i++)
   {
     gridData[i] *= kCAL_TO_kT_PER_MOL;
   }

   sz*=sizeof(float);

   fwrite(&sz, 4, 1, fp); /* grid header */

   fwrite(gridData, sizeof(float), szgcub, fp);

   fwrite(&sz, 4, 1, fp); /* grid footer */

   /* revert the grid so it doesn't get corrupted over time */
   sz = szgcub;
   for (i = 0; i < sz; i++)
   {
     gridData[i] *= kT_TO_kCAL_PER_MOL;
   }


   /* write the bottom label */
   /**************************/
   sz=16;
   fwrite(&sz, 4, 1, fp); /* header */
   fwrite("scale, midpoint ", 1, 16, fp);
   fwrite(&sz, 4, 1, fp); /* footer */

   /* write the scale and midpoint */
   /********************************/

   sz=16; /* sizeof(float) * 4 */
   fwrite(&sz, 4, 1, fp);    /* header */

   fwrite(&scale, 4, 1, fp); /* scale  */
   fwrite(&midx, 4, 1, fp);  /* midx   */
   fwrite(&midy, 4, 1, fp);  /* midy   */
   fwrite(&midz, 4, 1, fp);  /* midz   */

   fwrite(&sz, 4, 1, fp);    /* footer */

   /* close the file */
   fclose(fp);

   return 1;

}/* end function write_phi_grid*/
