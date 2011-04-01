#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "utils.h"
#include "blast.h"

void init_lookups(void);

unsigned char htinsert[4][256];
unsigned char htlookup[4][256];

void init_lookups(void)
{
   int a,b,c,d;
   int curval = 0;
   unsigned char htbytes[3][3][3][3];

   for (a = 0; a < 3; a++)
      for (b = 0; b < 3; b++)
         for (c = 0; c < 3; c++)
            for (d = 0; d < 3; d++)
            {
               htlookup[0][curval] = a;
               htlookup[1][curval] = b;
               htlookup[2][curval] = c;
               htlookup[3][curval] = d;
               htinsert[0][curval] = curval;
               htinsert[1][curval] = curval;
               htinsert[2][curval] = curval;
               htinsert[3][curval] = curval;
               htbytes[a][b][c][d] = curval;
               curval++;
            }

   curval = 0;
   for (a = 0; a < 3; a++)
      for (b = 0; b < 3; b++)
         for (c = 0; c < 3; c++)
            for (d = 0; d < 3; d++)
            {
               if (a < 2)
                  htinsert[0][htbytes[a][b][c][d]] = htbytes[a+1][b][c][d];
               if (b < 2)
                  htinsert[1][htbytes[a][b][c][d]] = htbytes[a][b+1][c][d];
               if (c < 2)
                  htinsert[2][htbytes[a][b][c][d]] = htbytes[a][b][c+1][d];
               if (d < 2)
                  htinsert[3][htbytes[a][b][c][d]] = htbytes[a][b][c][d+1];
            }

   return;
}

Hashcounter *hashcounter_new(unsigned long logsize, unsigned long seed)
{
   Hashcounter *new_counter;
   int c;

   init_lookups();
   new_counter = (Hashcounter *) global_malloc(sizeof(Hashcounter));

   for (new_counter->indmask = 0, c = 0; c < logsize; c++)
   {
      new_counter->indmask <<= 1;
      new_counter->indmask++;
   }
   new_counter->indmask <<= 2;
   new_counter->fieldmask = 3;
   
   new_counter->byte_size = 1 << logsize;
   /* We can fit 4 entries in each byte */
   new_counter->size = 4 * new_counter->byte_size;
   new_counter->hashseed = seed;
   new_counter->insertions = 0;
   
   new_counter->hashtable = global_malloc(new_counter->byte_size);

   printf("New hashcounter size=%.1f Mb\n", (float)new_counter->byte_size / 1024.0 / 1024.0);

   hashcounter_reset(new_counter);

   return new_counter;
}

void hashcounter_insert(Hashcounter *counter, uint4 hashValue)
{
   unsigned long ind, field;
   unsigned char cur;

   ind = (hashValue & counter->indmask) >> 2;
   field = hashValue & counter->fieldmask;

   cur = counter->hashtable[ind];
   counter->hashtable[ind] = htinsert[field][cur];

   counter->insertions++;

   return;
}

int hashcounter_multiple(Hashcounter *counter, uint4 hashValue)
{
   unsigned char cur;

   cur = counter->hashtable[(hashValue & counter->indmask) >> 2];
   return htlookup[hashValue & 3][cur] & 2;

}

int hashcounter_single(Hashcounter *counter, uint4 hashValue)
{
   unsigned char cur;

   cur = counter->hashtable[(hashValue & counter->indmask) >> 2];
   return htlookup[hashValue & 3][cur];
   
}

unsigned long hashcounter_insertions(Hashcounter *counter)
{
   return counter->insertions;
}

void hashcounter_write(Hashcounter *counter, FILE *stream)
{
   uint64_t c;

   for (c = 0; c < counter->size / 4; c++)
   {
      fputc(counter->hashtable[c], stream);
   }
}

void hashcounter_count(Hashcounter *counter, FILE *stream)
{
   uint64_t singles = 0, multiples = 0;
   uint64_t c = 0, d = 0;
   unsigned char cur;

   for (c = 0; c < counter->byte_size; c++)
   {
      cur = counter->hashtable[c];
      for (d = 0; d < 4; d++)
      {
         switch( htlookup[d][cur])
         {
            case 1: 
               singles++;
            case 0:
               break;

            default:
               multiples++;
         }
      }
   }

   fprintf(stream, "Hashcounter state:\nSingles: %"PRIu64"\nDuplicates: %"PRIu64"\n", singles, multiples);

   return;
}

void hashcounter_reset(Hashcounter *counter)
{

   counter->insertions = 0;
   memset(counter->hashtable, 0, counter->byte_size);

   return;
}

void hashcounter_free(Hashcounter *counter)
{
   free(counter->hashtable);
   free(counter);

   return;
}

