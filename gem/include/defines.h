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
#ifndef __defines_h__
#define __defines_h__

/**************************************************************************
 * This file contains all the defined constants to be used by calculation *
 * methods, if you want to optimize a constant or experiment with         *
 * constant values, this is the place.   --jcg                            *
 **************************************************************************/

/* various parameters for reading in options */
/* flag definitions */
#define NUM_FLAGS       14
#define READ_PHI        0
#define WRITE_PHI       1
#define VISUALIZE       2
#define PROJ_SPEC       3
#define TRI_SPEC        4
#define A_SPEC          5
#define DUMP_IMAGE      6
#define COMPARE_PHI     7
#define DUMP_VERTICES   8
#define PROBE_SPEC      9
#define FIX_COLOR_SCALE 10
#define DECOR_MOLECULE  11
#define DECOR_TYPE      12
#define BG_COLOR        13

/* parameter data definitions */
#define NUM_FILE_TYPES        15
#define PHI_IN_FILE           0
#define PHI_OUT_FILE          1
#define PHI_EDGE_VAL          2
#define PHI_WID_VAL           3
#define PROJ_LENGTH_VAL       4
#define TRI_SPEC_VAL          5
#define A_SPEC_VAL            6
#define IMAGE_VALS            7
#define PHI_COMPARE_FILE      8
#define VERTEX_FILE_NAME      9
#define PROBE_RADIUS_VAL      10
#define COLOR_SCALE_VAL       11
#define DECOR_MOL_NAME        12
#define DECOR_MOL_LOOK        13
#define BG_COLOR_RGB          14

/* Makes it easier to understand compilation flags */
#define YES 1
#define NO  0

#define NO_POTENTIAL 0
#define COULOMB_POTENTIAL 1
#define REACTION_POTENTIAL 2
#define TOTAL_POTENTIAL 3

/* Definitions for Onufriev-Fenley method */
#define ALPHA_OF .571412f /* "Incorporating variable dielectric environments into the generalized Born model" Sigalov, Scheffe, Onufriev J. Chem. Phys. */
#define min_A_mult .75
#define max_R_mult 1.5
#define max_R_mult2 max_R_mult * max_R_mult

/* convenient definitions */
#define SQRT_TWO 1.4142136 

#define ELECTROSTATIC_CONVERSION_FACTOR 332.

#define kT_TO_kCAL_PER_MOL .592
#define kCAL_TO_kT_PER_MOL 1.68918919

/* methodological switches */
#define USE_FLOATING_CENTERS NO /* yes or no answers here all caps */

/* yes causes the system to use an elliptic integral estimation method
   for estimating A, NO causes the system to use the "simple" method that
   is rumored to be more accurate
 */
#define USE_ELLIPTIC_INTEGRAL NO

/* This is the MCLA coarse-graining threshold. As long as it is larger */
/* Than the system size, it is effectively OFF */ 
#define GRANULAR_CUTOFF 1000.0
#define GRANULAR_CUTOFF2 1000000.0    /* cutoff**2 to reduce computations */

#endif /* __defines_h__ */
