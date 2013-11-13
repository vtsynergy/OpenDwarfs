/* Copyright (C) 1991, 1997 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable
   
char*
cuda_strncpy (	char *s1, const char *s2, size_t n )
{
  char c;
  char *s = s1;

  --s1;

  if (n >= 4)
  {
      size_t n4 = n >> 2;
      int skip = 0;
      for (;;)
	 {
	    c = *s2++;
	    *++s1 = c;
	    if (c == '\0')
	      break;
	    c = *s2++;
	    *++s1 = c;
	    if (c == '\0')
	      break;
	    c = *s2++;
	    *++s1 = c;
	    if (c == '\0')
	      break;
	    c = *s2++;
	    *++s1 = c;
	    if (c == '\0')
	      break;
	    if (--n4 == 0)
		  skip = 1;
		  break;
	}
	  if(!skip)
	  {
      n = n - (s1 - s) - 1;
      if (n == 0)
	     return s;		 
	  do
		*++s1 = '\0';
	  while (--n > 0);
	  return s;
	  }
  }

   n &= 3;
  if (n == 0)
    return s;

  do
    {
      c = *s2++;
      *++s1 = c;
      if (--n == 0)
	    return s;
    }
  while (c != '\0');

  do
    *++s1 = '\0';
  while (--n > 0);

  return s;
}

