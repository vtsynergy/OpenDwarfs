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
#include <math.h>
#include "file_io.h"
#include "calculations.h"
#include <stdlib.h> /* for exit */
#include <string.h> /* for strcat */

#define MSMS_LINE_SIZE 1024

/**************************************************************************
 * FUNCTION: read_msms_output  --reads in files output by an msms run     *
 *                                                                        *
 * INPUTS:   molname   -- the name of the molecule to run on              *
 *                                                                        *
 * OUTPUTS:  vertices  -- the vertices read from the .vert file           *
 *           triangles -- the triangles read from the .face file          *
 *           nvert     -- the number of vertices read in                  *
 *           nface     -- the number of triangles read in                 *
 *                                                                        *
 *                                                                        *
 * RETURNS: nothing                                                       *
 **************************************************************************/
void read_msms_output(char *molname, vertx **vertices, triangle **triangles, int *nvert, int *nface)
{
   /* local variables */
   FILE *fp; /* the file pointer to the vertex file */

   double mag;
   char line[MSMS_LINE_SIZE]; /* a line in the input file to read */

   int  i,         /* counter for iterations */
        ierr,      /* the number of things read in a sscanf or scanf        */
        nv,
        nt;

   vertx  *vertex;
   triangle *tri;

   *vertices = NULL;
   *triangles = NULL;
   *nvert = *nface = 0;

   /* Open the vert file */
   strcpy(line, molname);
   strcat(line, ".vert");
   fp = fopen(line,"r");
   if (!fp) 
   {
      fprintf(stderr, "Error: could not open file %s\n", line);
      return;
   }

   /* Have a look at the header records in the vertex file.          */
   /* I suspect that the 2nd line should be parsed to determine the  */
   /* actual file contents, but the current documentation doesn't    */
   /* indicate what other items might be on the line.  So for now    */
   /* I'll assume that it must contain the number of vertices as the */
   /* first item on the line.                                        */
   if ((fgets( line, MSMS_LINE_SIZE, fp ) == NULL) ||
       (fgets( line, MSMS_LINE_SIZE, fp ) == NULL) ||
       (fgets( line, MSMS_LINE_SIZE, fp ) == NULL)) {
       fprintf(stderr, "Error: unable to read from file\n");
       return;
   }
   sscanf( line, "%d %d", &nv, &ierr);

   /* allocate the vertex array */
   vertex = (vertx *)calloc(nv, sizeof(vertx));

   /* check for errors */
   if(vertex == NULL)
   {
      fprintf(stderr, "Error, too many vertices to process on this machine\n");
      return;
   }

   for (i=0; i< nv; i++)
   {
       if (fgets( line, MSMS_LINE_SIZE, fp ) == NULL) {
           fprintf(stderr, "Error: unable to read from file\n");
           return;
       }
      ierr = sscanf
             (
                 line,
                 "%f %f %f %f %f %f %*d %d %*f", /* this is the full format of the line */
                 &vertex[i].x,
                 &vertex[i].y,
                 &vertex[i].z,
                 &vertex[i].xNorm,
                 &vertex[i].yNorm,
                 &vertex[i].zNorm,
                 &vertex[i].nearest_atom
             );

      /* make sure normals are unit normals */
      mag = dist(vertex[i].xNorm, vertex[i].yNorm, vertex[i].zNorm, 0., 0., 0.);

      vertex[i].xNorm /= mag;
      vertex[i].yNorm /= mag;
      vertex[i].zNorm /= mag;

      /* MSMS starts with 1 instead of 0 */
      vertex[i].nearest_atom -= 1;

      if (ierr != 7)
      {
         fprintf(stderr, "Error reading vertex %d", i);
         free(vertex);
         return;
      }

    }/* end for (nv) loop */

    /* close the vertex file */
    fclose(fp);

    strcpy(line, molname);
    strcat(line, ".face");
    fp = fopen(line,"r");

    /* read header information */
    if ((fgets(line, MSMS_LINE_SIZE, fp) == NULL) ||
        (fgets(line, MSMS_LINE_SIZE, fp) == NULL) ||
        (fgets(line, MSMS_LINE_SIZE, fp) == NULL)) {
        fprintf(stderr, "Error: unable to read from file\n");
        return;
    }
    sscanf(line, "%d", &nt);

    /* allocate our triangles here */
    tri = (triangle *)calloc(nt, sizeof(triangle));

    /* now start reading triangles */
    for (i = 0; i < nt; i++)
    {
        if (fgets(line, MSMS_LINE_SIZE, fp) == NULL) {
            fprintf(stderr, "Error: unable to read from file\n");
            return;
        }
       ierr = sscanf(line,"%d %d %d %*d %*d", &tri[i].v[0], &tri[i].v[1], &tri[i].v[2]);

       tri[i].v[0]--;
       tri[i].v[1]--;
       tri[i].v[2]--;

       if (ierr != 3) /* check for errors while reading */
       {
          fprintf(stderr, "Error reading from face file %s.face\n", molname);
          free(vertex);
          free(tri);
          return;
       }

       /* assign each triangle a nearest atom in a logical way */
       tri[i].nearest_atom = vertex[tri[i].v[0]].nearest_atom;
    
    } /* end for (i < nt) */

    /* set output pointers to data retrieved */
    *nvert = nv;
    *nface = nt;
    *triangles = tri;
    *vertices = vertex;

    printf("Molecular surface consists of %i vertices and %i triangles.\n", *nvert, *nface);
    fflush(stdout);

    /* close the face file */
    fclose(fp);

}/* end read_msms function */
