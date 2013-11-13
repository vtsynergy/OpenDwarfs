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

#ifndef __VISUALIZATION_H__
#define __VISUALIZATION_H__

/* killing myself to get the GLhandleARB type */
#include "GLee.h"
#include <GL/glx.h>

#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include <Xm/Xm.h>
#include "status_dialog.h"
#include "structures.h"
#include "residue_select_form.h"
#include "atom_info_dialog.h"

#define READ_ONLY       0
#define GEM_ONLY        1
#define READ_DIFF_READ  2
#define GEM_DIFF_READ   3
#define DIFF_FLAG       2
#define CALC_FLAG       1

#define SURFACE_CM 0
#define ATOM_CM    1

/* definitions defining the current type of mouse interraction */
#define M_NONE   0
#define M_MOVE   1
#define M_ZOOM   2
#define M_ROTATE 3
#define M_SELECT 4

/* definitions for defining various selections for drawing */
#define DRAW_NONE          0
#define S_DRAW             1
#define A_DRAW_STICKS      1
#define S_DRAW_WIREFRAME   2
#define A_DRAW_BALL_STICK  2
#define A_DRAW_SPACEFILL   3
#define A_DRAW_BACKBONE    4

typedef struct
{
   residue *residues; /* residues making up the molecule    */
   int  nresidues; /* number of residues in the molecule */

   bond *bonds; /* the bonds in the molecule */
   int  nbonds; /* number of bonds in the molecule */

   vertx *vert;  /* points to a dynamic array of vertices */
   int  nvert;  /* number of vertices in the array       */

   triangle *tri;  /* points to dynamic array of triangles   */
   int  ntri;  /* number of triangles making up surface  */

   int atomDrawType;
   int surfaceDrawType;

}decor_mol_data_struct;


typedef struct
{
   char *molname; /* the molecule name we are viewing */

   /* this chunk supports grid writing/subsampling */
   /************************************************/
   void  *grid;       /* Grid containing either floats or doubles
                         representing potential in space, this is
                         effectively always a float, though the 
                         AVS format supports doubles.  Will need a
                         flag and some ifs if MEAD ever uses this
                         capacity within AVS */

   int   grid_x_dim,  /* number of cells in X dimension of grid */
         grid_y_dim,  /* number of cells in Y dimension of grid */
         grid_z_dim;  /* number of cells in Z dimension of grid */

   /*******************/
   /* grid chunk over */

   float grid_world_dims[6]; /* xmin, ymin, zmin, xmax, ymax, zmax */

   int phiType;

   double diel_int,  /* internal dielectric constant */
          diel_ext,  /* external dielectric constant */
          ion_exc_rad, /* ion exclusion radius about the molecule */
          A,           /* electrostatic radius of the molecule */
          sal;       /* salinity */

   residue *residues; /* residues making up the molecule    */
   int  nresidues; /* number of residues in the molecule */

   bond *bonds; /* the bonds in the molecule */
   int  nbonds; /* number of bonds in the molecule */

   vertx *vert;  /* points to a dynamic array of vertices */
   int  nvert;  /* number of vertices in the array       */

   triangle *tri;  /* points to dynamic array of triangles   */
   int  ntri;  /* number of triangles making up surface  */

   int calc_type;  /* type of calculation to be done */
                   /* either GEM_ONLY      */
                   /*        GEM_DIFF_READ */
                   /* or     READ_ONLY     */

} mol_data_struct;

typedef struct
{
   GLuint vertex_shader,    /* vertex shader id   */
          fragment_shader,  /* fragment shader id */
          program;          /* program id         */

}GL_20_shader_def_struct;

typedef struct
{
   GLhandleARB vertex_shader,    /* vertex shader id   */
               fragment_shader,  /* fragment shader id */
               program;          /* program id         */

}GL_ARB_shader_def_struct;

typedef struct
{
   int type;  /* 2 is gl2, 1 is arb */

   union
   {
       GL_20_shader_def_struct  gl2;
       GL_ARB_shader_def_struct arb;
   }in;

}GL_shader_def_struct;

typedef struct
{
   int  selected_res; /* the selected residue in the structure */
   /* vertex potential bounds */
   float  vp[4], /* min, max, extent, and center of potentials */
          vc[4], /* min, max, extent, and center of charges */

   /* vertex surface x, y, z location bounds */
           x[4], /* min, max, extent, and center of x bounds    */
           y[4], /* min, max, extent, and center of y bounds    */
           z[4]; /* min, max, extent, and center of z bounds    */

   Widget     toplevel, /* toplevel widget for the whole system */
              drawTop,  /* toplevel for the drawing area */
              drawW,   /* drawing area widget         */
              colrW;   /* the color scale widget      */

   double xshift, yshift,   /* shifts applied    */
          trans,            /* transparency      */
          scale;            /* scaling applied   */

    int     surfDrawMode,   /* how do we draw the surface? */
            atomDrawMode,   /* how do we draw the atoms?   */
            colorType;      /* color by charge or potential? Charge = 1, pot = 0 */

    int     colorMapType,   /* either SURFACE_CM, or ATOM_CM */
            start_x,        /* the mouse x position when button was pressed */
            start_y,        /* the mouse y position when button was pressed */
            click_x,        /* the x location where the mouse was clicked */
            click_y,        /* the y location where the mouse was clicked */
            mouse_mode,     /* either MOVE, ZOOM, or ROTATE */
            last_mode,      /* support for meta-clicking    */
            draw_wid,       /* dimensions of the main draw widget */
            draw_hit,
            rotation_dirty; /* tells us if we have rotated since we last sorted */

    atom_info_dialog_struct *atom_info; /* the atom info dialog (always up) */

    double rot_mat[16],           /* the rotation matrix to use      */
           rot_store_mat[16];     /* the stored last rotation matrix */

    char up; /* tells us whether the interface is up or not */

    XtAppContext app;  /* app context for looping in the xtAppMainEventLoop etc   */
    XVisualInfo *vi;   /* visual info struct, defines capabilities of the display */

    GLXContext context_D,  /* context for the drawing area */
               context_S;  /* context for the scale        */

    /* display lists for atoms and surface */
    GLuint atomList,     /* atom display list for main molecule */
           surfaceList,  /* surf display list for main molecule */
           decorList;    /* decoration molecule display list    */
           
    GL_shader_def_struct phong_shader;  /* shaders if enabled */

    XFontStruct *font; /* contains our fontinfo from XLoadQueryFont */

    residue_selection_struct *res_sel; /* always up */

    float bg_color[3]; // r = 0, g = 1, b = 2

} vis_param_struct;


/* this structure just encapsulates the parameters used by
   visualize and its subfunctions */
typedef struct 
{
   mol_data_struct mol_data;/* data regarding the molecule being viewed */

   decor_mol_data_struct decor_mol; /* molecule to be shown for decoration only */

   vis_param_struct params; /* data regarding view parameters, etc.  */

} vis_data_struct;


#ifdef __cplusplus /* make it a c function so c can access it */
extern "C" {
#endif

/* gets an initialized visual structure */
vis_data_struct *get_vis_struct(residue*, vertx*, triangle*, int, int, int);

/* frees a visual structure */
void free_vis_struct(vis_data_struct *);
                                                                                
/* creates the menubar */
Widget createMenubar (Widget, vis_data_struct *);

/* creates all the windows and starts the visual interface */
void visualize (vis_data_struct *vis, int argc, char **argv);

/* calculates information like min, max, extent of everything */
void populate_statistics(vis_data_struct *);

/* draws to a pixmap and then dumps a screenshot to file */
int screenshot(vis_data_struct *vis, int width, int height, char *fname);

/* initializes the fonts for the colorbar */
XFontStruct *initialize_gl_fonts(Display *);

/* generates color information for everything at the beginning */
void generate_colors(vis_data_struct *);

/* updates all zdepths on the surface */
void update_zdepths (vertx *vert, int nvert, triangle *tri, int ntri);

/* opens a decor file */
int open_decoration_molecule(char *, char *, vis_data_struct*);

/* regenerates main molecule display lists to trim that code out */
void gen_main_mol_lists (vis_data_struct *, char, char, double);

/* regenerates decor molecule display list to trim that code out */
void gen_decor_mol_list (vis_data_struct *);

/* utility routines used throughout */
/************************************/
void rotate_to_vector (vis_data_struct *, float, float, float);
void updateTranslations (vis_data_struct *, int);
void set_phong_shader (GL_shader_def_struct *shader);
void write_grid_from_vis (char *, int, float, vis_data_struct *);

/* draw routines for various models */
void draw(vis_data_struct *vis);
void drawAtoms(residue *, int, bond *, int, int, int);
void drawSelection(residue *, int);
void drawSurface(vertx *, int, triangle *, int, residue *, int, int, int, double);
void drawCylinderVec(GLUquadricObj*,float,float,float,float,float,float);
void drawColrMap(vis_data_struct *vis);         /* draws to whatever is current, swaps buffers */
void just_draw_color_map(vis_data_struct *vis, int, int); /* draws to whatever is current, no swaps, etc */
void normalize(float *, float *);
int  find_nearest_res (int, int, vis_data_struct *);
void rotate (int, int, vis_data_struct *);

/* takes index into array */
void find_res_atm_pair(residue *, int, int, int *, int *);
int  find_linear_index(residue *, int, int);

/* various dialogs */
/*******************/
/* prompts for data, and sets the color scale */
void set_color_scale_dialog(Widget, vis_data_struct*);
/* prompts for pertinent information, then opens a file */
void open_file_dialog(Widget, vis_data_struct*);
/* prompts for grid stuff, then saves the grid */
void save_grid_dialog(Widget, vis_data_struct*);
/* prompts for width, height, filename, then writes the file */
void save_image_dialog(Widget, vis_data_struct*);
/* same as open_file_dialog but for a deocoration molecule */
void open_decor_file_dialog(Widget, vis_data_struct*);

void clear_moldata (mol_data_struct *);
void clear_decor_mol (decor_mol_data_struct *);

#ifdef __cplusplus
}
#endif

#endif
