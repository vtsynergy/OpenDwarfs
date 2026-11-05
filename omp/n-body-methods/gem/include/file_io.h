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
#ifndef __READWRITE_H__
#define __READWRITE_H__

/*********************************************************************
 * Header file for i/o routines and their associated data structures *
 * --jcg                                                             *
 ********************************************************************/

#include "structures.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* writes out an image */
/***********************/
void SaveScreenShot(void *, char *,int,int);

/* writes out tga image header */
void TGA_header(FILE *, int, int);

/* write bytes to a tga image captured from ogl */
void TGA_bytes(FILE *, int, int, unsigned char *);

/* reads in atom locations and returns number of atoms read */
int read_pqr(char*, residue**);

/* reads in msms surface files */
void read_msms_output (char*, vertx**, triangle**, int*, int*);

/* run msms */
int runMsms(char*, residue*, int, double, double);

/* writes a phimap of calculations */
int write_phi_grid (char*, int, float*, float, float, float, float);

/*reads in a phimap of precalculated potentials */
int read_phi_grid (char*, float**,float*,float*, float*, float*,float*,float*);

/* reads in an AVS grid file as a potential map */
int read_avs(char*,float**,int*,int*,int*,float*,float*,float*,float*,float*,float*);

/* writes out an AVS grid file as a potential map */
int write_avs(char*, float*, int, float, float, float, float, float, float);

/* writes out a grid and determines type from filename */
int write_grid (char *, int, float *, float, float, float, float);

/* writes an xyzr file for MSMS */
void write_xyzr(char*, residue*, int);

/* writes out a text file containing vertex information */
void dump_vertices(char *, int, vertx *);

/* batch execution of the partitioned_open function, populates the vis_data_struct and runs all the usual computations --jcg */
int openFile(char*,int,char*,double,double,double,double,double,double,double,int,char*,void*);

/* byte swaps an array of data (4 bytes long) of size (size) */
void swap4(void *data, int size);

#ifdef __cplusplus
}
#endif

#endif
