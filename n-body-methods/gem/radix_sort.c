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
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "calculations.h"

#include <stdio.h>

/* this version of radix sort works on arbitrary structures
 * by sorting by a given key of a given length in the standard
 * manner.  -- Will not operate on C++ objects.
 */

uchar big_endian (void)
{
	union {long l; uchar c[sizeof(long)]; } u;

	u.l = 1;
	return (u.c[sizeof(long) -1]);
}

/***********************************************************************
 * FUNCTION: radix    --sorts a given list of structures by radix sort * 
 *                         on the key at keyloc bytes into the objects *
 *                                                                     *
 * INPUTS:   byte      --the given byte by which we will radix sort    *
 *           keyloc    --the byte offset of the key location in the    *
 *                           individual objects                        *
 *           member_size -- size of each member in the array           *
 *           members     -- the number of members in the array         *
 *           source      -- the data source (unsorted)                 *
 *                                                                     *
 * OUTPUTS:  dest        -- the data destination (will be sorted)      *
 *                                                                     *
 ***********************************************************************/
static void radix (ushort byte, ushort keyloc, ushort member_size,
		ulong members, void *source, void *dest)
{
	/* local variables */
	/* count will be an accumulation array (serving as an index for each list) */
	ulong count[256];

	ulong   sum, /* sum for summing elements into the accumulation buffer */
		c,   /* character value for summing into accumulation buffer  */
		*cp,   /* count position (offset into count array)              */
		i;   /* counter for various loops */

	uchar *sp, /* source pointer      */
	      *dp; /* destination pointer */

	uchar *bp, /* byte pointer for our byte by which we are sorting */
	      *src  = (uchar *) source,
	      *dst  = (uchar *) dest;

	/* fewer increment operations than memset and no function overhead */
	cp = count;
	for (i = 256; i > 0; --i, ++cp)
		*cp = 0; /* initialize count array to 0 */

	/* count occurences of every byte value */
	/****************************************/

	/* sort by the byte at "byte" into the key at "keyloc" into the elem */
	bp = src + keyloc + byte;
	for (i = members; i > 0; --i, bp += member_size)
	{
		++count[*bp];
	}

	/*
	 * transform count into index by summing elements 
	 *********************************************************************
	 * this may seem like a strange way to pass through
	 * an array to some people, but typical array dereferences
	 * cost more time than this becuase they do "base + count*elem_size"
	 * which is significantly more work per index than what we are doing
	 * here.
	 */
	sum = 0; /* initialize sum to 0 */
	cp = count;
	for (i = 0; i < 256; ++i, ++cp)
	{
		c = *cp;
		*cp = sum;
		sum += c;
	}

	/* fill dest with the right values in the right place */
	/******************************************************/

	/* byte pointer to index into the count array */
	bp = src + keyloc + byte;
	sp = src;
	for (i = members; i > 0; --i, bp += member_size, sp += member_size)
	{
		cp = count + *bp;

		dp = dst + ((*cp) * member_size);

		memcpy(dp, sp, member_size);

		++(*cp);
	}
}

/***********************************************************************
 * FUNCTION: radix_sort  --sorts a given list of structures by radix   * 
 *                         on the key at keyloc bytes into the objects *
 *                                                                     *
 * INPUTS:   keylen    --the length of the key in bytes                *
 *           keyloc    --the byte offset of the key location in the    *
 *                           individual objects                        *
 *           m_size    --size of each member in the array              *
 *           members   -- the number of members in the array           *
 *           to_sort   -- the data source (unsorted)                   *
 *                                                                     *
 * OUTPUTS:  to_sort   -- will be sorted                               *
 *                                                                     *
 ***********************************************************************/
void radix_sort (ushort keylen, ushort keyloc, ushort m_size, ulong members, void *to_sort)
{
	/* local variables */
	void *temp = malloc (members * m_size);
	void *src, *swp, *tmp;
	short i, start, end;
	short increment;

	assert (temp != NULL);

	src = to_sort;
	tmp = temp;

	if (big_endian())
	{
		start = keylen-1;
		end = -1;
		increment = -1;
	}
	else
	{
		start = 0;
		end = keylen;
		increment = 1;
	}

	for (i = start; i != end; i+=increment)
	{
		radix (i, keyloc, m_size, members, src, tmp);

		swp = tmp;
		tmp = src;
		src = swp;
	}

	if ((keylen % 2) != 0)
	{
		memcpy (to_sort, temp, members * m_size);
	}

	free (temp);
}
