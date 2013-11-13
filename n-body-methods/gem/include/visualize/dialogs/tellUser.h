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
#ifndef __TELLUSER_H__
#define __TELLUSER_H__
/***************************************************************
 * This is the header to use the tellUser function to tell the *
 * user about warnings or errors that may occur during regular *
 * operation.                                  --jcg           *
***************************************************************/
#include <Xm/Xm.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TELLUSER_WIDTH  575
#define TELLUSER_HEIGHT 60

/* global values */
extern Widget mainForm;  /* created in  init.c */

/* general function to print important information to screen */
void tellUser(char *windowTitle, char *thing);

#ifdef __cplusplus
}
#endif

#endif
