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
#ifndef __GEM_MAIN_H__
#define __GEM_MAIN_H__

#include "structures.h"
#ifndef NO_UI
#include "visualize.h"
#endif

#ifdef __cplusplus /* it wont link with C, so go away C (string class) */

#include <string>

using namespace std;

/* checks command line arguments */
bool check_cmd_line (int, char**, string*, bool*);

/* checks to see if a string is a number */
int stringIsNum(string);

extern "C" {
#endif

    int cpp_main(int, char**);
/* place cross C/C++ functions here */


#ifdef __cplusplus
}
#endif

#endif
