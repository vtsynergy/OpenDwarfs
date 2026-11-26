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

#ifndef __VERTEX_SHADER_H__
#define __VERTEX_SHADER_H__

/* both shaders from GLSL Tutorial on lighthouse3d
http://www.lighthouse3d.com/opengl/glsl/index.php?ogldir1
*/

#define VERTEX_SHADER_CODE \
" \
varying vec4 diffuse, ambientGlobal, ambient; \
varying vec3 normal, lightDir, halfVector; \
varying float dist; \
\
void main() \
{ \
    vec4 ecpos; \
    vec3 aux; \
\
    normal = normalize(gl_NormalMatrix * gl_Normal); \
    \
    ecPos = gl_ModelViewMatrix * gl_Vertex; \
    aux = vec3(gl_LightSource[0].position-ecPos); \
    lightDir = normalize(aux); \
    dist = length(aux); \
    \
    halfVector = normalize(gl_LightSource[0].halfVector.xyz); \
    diffuse = gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse; \
    \
    ambient = gl_FrontMaterial.ambient * gl_LightSource[0].ambient; \
    ambientGlobal = gl_LightModel.ambient * gl_FrontMaterial.ambient; \
    \
    gl_Position = ftransform(); \
}"

#endif
