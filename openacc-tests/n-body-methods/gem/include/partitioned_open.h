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
#ifndef __PARTITIONED_OPEN_H__
#define __PARTITIONED_OPEN_H__

#ifndef NO_UI
#include "visualize.h"
#else
#define Widget void* 
#define vis_data_struct void* 
#endif

#include "structures.h"

#define POPEN_CALC_NONE     0 
#define POPEN_CALC_DIRECTLY 1
#define POPEN_READ          2

typedef struct
{
    int  calctype, /* do we wish to calculate this directly?    */
         compare_phi,   /* do we wish to compare values to a phimap? */
         phiType;       /* type of potential to calculate            */

    double diel_int,    /* internal dielectric constant */
           diel_ext,    /* external dielectric constant */
           ion_exc_rad, /* ion exclusion radius         */
           sal,         /* salinity */
           proj_len,    /* projection length */
           probeRadius, /* probe radius for MSMS     */
           triDens;     /* triangle density for MSMS */

    float  vp[3];       /* used to keep a hold on the potential bounds
                           for comparison mapping where values are only
                           intermediate for the vertices */

    vis_data_struct *vis; /* visual parameters and information */

    /* mirrored temporary pointers to prevent corruption during process */
    /********************************************************************/
    int nvert,       /* number of vertices on the surface  */
        nresidues,   /* number of residues in the molecule */
        i,           /* current progress counter           */
        ntri;        /* number of faces on the surface     */

    residue   *residues;     /* residues making up the molecule          */
    vertx   *vert;           /* vertices of the triangles on the surface */
    triangle *tri;           /* faces of the triangles on the surface    */

    float *grid;             /* grid values read from a phimap file      */

    double A;                /* Radius of molecule if it were a sphere   */

    char *molname;           /* molecule name                  */
    char *auxStr,            /* auxiliary file information     */
         *compare_phi_file;  /* phimap file to compare data to */

} partitioned_open_struct;


#ifdef __cplusplus
extern "C" {
#endif

/* encapsulates the process of opening via a dialog */
void popen_file (Widget, char*, int, char*, double, double, double, double, double, double, double, int, double, vis_data_struct *); 

/* a single step of the partioning process */
int partitioned_open (partitioned_open_struct *, int, int);

/* what to do when a partitioned open finishes */
void partitioned_open_finalize (partitioned_open_struct *);

/* function just to open a pqr file and run msms */
int open_pqr_run_msms(partitioned_open_struct *);

/* function to open a decoration molecule */
int open_decoration_molecule(char *, char *, vis_data_struct *);

#ifdef __cplusplus
}
#endif

#endif
