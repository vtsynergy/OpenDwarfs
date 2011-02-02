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

#ifndef __FRAGMENT_SHADER_H__
#define __FRAGMENT_SHADER_H__

/* both shaders from GLSL Tutorial on lighthouse3d
http://www.lighthouse3d.com/opengl/glsl/index.php?ogldir1
*/

#define FRAGMENT_SHADER_CODE \
" \
varying vec4 diffuse, ambientGlobal, ambient; \
varying vec3 normal, lightDir, halfVector; \
varying float dist; \
void main() \
{ \
vec3 n, halfV, viewV, ldir; \
float NdotL, NdotHV; \
vec4 color = ambientGlobal; \
float att; \
\
   n = normalize(normal); \
   NdotL = max(dot(n, normalize(lightDir)), 0.0); \
\
   if (NdotL > 0.0) {\
     att = 1.0 / (gl_LightSource[0].constantAttenuation + \
                  gl_LightSource[0].linearAttenuation * dist + \
                  gl_LightSource[0].quadraticAttenuation * dist * dist); \
     color += att * (diffuse * NdotL + ambient); \
\
     halfV = normalize(halfVector); \
     NdotHV = max(dot(n, halfV), 0.0); \
     color += att * gl_FrontMaterial.specular * gl_LightSource[0].specular * \
                    pow(NdotHV, gl_FrontMaterial.shininess);\
   }\
\
gl_FragColor = gl_Color; \
}"
#endif
