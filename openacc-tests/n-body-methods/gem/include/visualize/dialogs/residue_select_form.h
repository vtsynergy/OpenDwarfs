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

#ifndef __RESIDUE_SELECT_WINDOW_H__
#define __RESIDUE_SELECT_WINDOW_H__

#include <Xm/Xm.h> /* for widget type */

typedef struct {

   Widget rowCol;        /* the row column parent widget     */

   int selected_residue; /* the selected residue index       */

   void *vis_ptr;        /* the visualization struct pointer */

   Pixel select_bg_color, /* keep track of pixel values we will be using */
         normal_bg_color; 

} residue_selection_struct;

#ifdef __cplusplus
extern "C" {
#endif

/* creates the initial (empty) residue selection form */
residue_selection_struct *create_residue_select_form (Widget, void *);

/* fills the buttons in the residue select form */
void fill_residue_select_form (residue_selection_struct *);

/* clears the buttons in the residue select form */
void clear_residue_select_form (residue_selection_struct *);

/* sets the selected residue in the residue selection struct */
void residue_select_form_set_selected (residue_selection_struct *, int, int);

#ifdef __cplusplus
}
#endif

#endif
