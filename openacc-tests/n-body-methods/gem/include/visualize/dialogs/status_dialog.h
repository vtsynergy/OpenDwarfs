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
#ifndef __STATUS_DIALOG_H__
#define __STATUS_DIALOG_H__

#include <Xm/Xm.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
   Widget toplevel,  /* the toplevel widget                     */
          statusBar, /* the drawing area for the status bar     */
          progressW; /* the widget with the progress percentage */

   char   *label;    /* the label for this run      */
   int    progress;  /* progress of the run         */

   void *userData;   /* data needed to do the run   */

   int timer,        /* the time between executions of the job      */
       lastUp;       /* last time we updated the bar and text field */

   int  (*func)  (void *, int); /* the function to run */
   void (*clean) (void *);      /* the cleanup function for userData */

   GLXContext con;   /* context for the status bar  */

   XtWorkProcId      workId; /* the id for this run */

} status_dialog_struct;

/* associated with status dialogs */
/**********************************/
/* creates the status dialog */
status_dialog_struct *create_status_dialog(Widget, void *,int (*f)(void *, int), char *, XVisualInfo *);

/* starts a run through the status dialog */
void status_start_run (XtAppContext, status_dialog_struct *);

/* pauses a run */
void status_pause_run (status_dialog_struct *);

/* destroys the dialog */
void destroy_status_dialog(status_dialog_struct *);

/* assigns a finalization routine */
void status_set_finalize(status_dialog_struct *, void (*c)(void *));

#ifdef __cplusplus
}
#endif

#endif
