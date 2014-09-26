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
#ifndef __GEM_STRUCTURES_H__
#define __GEM_STRUCTURES_H__

/* definitions for array indices corresponding to specific bounding
 * values */
#define MIN    0
#define MAX    1
#define EXTENT 2
#define CENTER 3

#define POTENTIAL   0
#define ATM_CHARGE  1

#include <sys/types.h>
//typedef unsigned short ushort;
//typedef unsigned long  ulong;
typedef unsigned char  uchar;

typedef struct atom
   {
      float  x, y, z,        /* position in space, x, y, z                */
             radius,         /* radius of the atom                        */
             charge;         /* charge of the atom                        */

      float pot_est;        /* an estimate of the potential of the atom  */

      float  r,g,b;          /* red, green, and blue color components     */

      char   name[5];        /* the name of this atom                     */
      int    index;          /* atom index from the input file            */

   } atom;

typedef struct residue
   {

     char name[5];  /* name of the residue         */

     int res_num,   /* this residue number         */

          natoms;   /* number of atoms in the residue */
     atom *atoms;   /* pointer to array of atoms   */

     float x, y, z, rad;

     float total_charge; /* total charge of the residue */

   } residue;

typedef struct vertx
   {
      ushort zdepth;  /* the z depth of this vertex */

      float x, y, z,             /* position in space, x, y, z            */
            xNorm, yNorm, zNorm, /* the normal at this point              */
            r,g,b;               /* red, green, and blue color components */

      float potential;           /* potential at the vertex               */
      float surface_charge;      /* surface charge at the vertex          */

      int nearest_atom;          /* the atom nearest to this vertex       */

   } vertx;

typedef struct triangle
   {
      ushort zdepth;  /* the z depth of this triangle (average) */
      int v[3]; /*the three vertices making up a triangle */
      int nearest_atom;  /* nearest atom by consensus */
   } triangle;

/* a bond exists between two atoms, so indexes for each atom in the bond */
typedef struct bond
   {
      int a1, /* atom 1    */
          r1, /* residue 1 */
          a2, /* atom 2    */
          r2; /* residue 2 */

      float length;
   } bond;
#endif
