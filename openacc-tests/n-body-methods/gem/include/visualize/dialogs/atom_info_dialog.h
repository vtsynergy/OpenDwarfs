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

#ifndef __ATOM_INFO_DIALOG_H__
#define __ATOM_INFO_DIALOG_H__

/* the dialog box that has the atom info in it */
typedef struct
{
   Widget topLevel,         /* toplevel widget in the dialog       */
          residueLabel,     /* label for the name of the residue   */
          atomRowCol,       /* the row column widget for the atoms */
          totalChargeLabel; /* displays the total charge of the residue */

   int popped_up,         /* whether or not the interface is up */
       current_res;       /* residue being displayed right now  */

   void *vis_ptr;

}atom_info_dialog_struct;

#ifdef __cplusplus
extern "C" {
#endif

/* associated with the atom info box */
/*************************************/
/* creates the box */
atom_info_dialog_struct *create_atom_info_dialog(Widget, void *vis);
/* pops up the box */
void popup_atom_info_dialog(atom_info_dialog_struct *);
/* pops down the box */
void popdown_atom_info_dialog(atom_info_dialog_struct *);
/* updates the fields in the box */
void update_atom_info_dialog(atom_info_dialog_struct *);

#ifdef __cplusplus
}
#endif

#endif
