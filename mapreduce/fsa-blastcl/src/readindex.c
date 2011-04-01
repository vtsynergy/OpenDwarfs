#include "readindex.h"
#include "utils.h"
#include "vbyte.h"
#include "_vec.h"
#include "vec.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

#define GRAMLEN_MAX 256

struct index_scanner *open_index(char *index_prefix)
{
   struct index_scanner *scanner;
   char fnamebuf[FILENAME_MAX + 1];

   fnamebuf[FILENAME_MAX] = '\0';

   scanner = (struct index_scanner *) malloc(sizeof(struct index_scanner));

   snprintf(fnamebuf, FILENAME_MAX, "%s.idx", index_prefix);
   scanner->idx_file = fopen(fnamebuf, "r");

   snprintf(fnamebuf, FILENAME_MAX, "%s.lo", index_prefix);
   scanner->offsets = fopen(fnamebuf, "r");

   return scanner;
}

struct listinfo *get_next_list(struct index_scanner *iscn, struct listinfo *info)
{
   return get_list_at(iscn, info, -1);
}

struct listinfo *get_next_biglist(struct index_scanner *iscn, struct listinfo *info)
{
   struct lo_info *lo_info;

   lo_info = get_next_lo_info(iscn, NULL);

   if (!lo_info)
      return NULL;

   return get_list_at(iscn, info, lo_info->offset);
}

struct lo_info *get_next_lo_info(struct index_scanner *iscn, struct lo_info *lo_info)
{
   FILE *offsetfile = iscn->offsets;

   if (!lo_info)
      lo_info = malloc(sizeof(struct lo_info));

   vbyte_read(offsetfile, &(lo_info->size));
   if (!(vbyte_read(offsetfile, &(lo_info->offset))))
   {
      free(lo_info);
      return NULL;
   }

   return lo_info;
}

struct listinfo *get_list_at(struct index_scanner *iscn, struct listinfo *info, long offset)
{
   struct vec *vector;
   unsigned long num, last_doc = -1;
   int c, d;
   FILE *idx = iscn->idx_file;
   unsigned long veclen,
                 numdocs,
                 occurs,
                 last_offset;

   if (offset >= 0)
      fseek(idx, offset, SEEK_SET);

   /* There is nothing left to read; free structures */
   if (!( vbyte_read(idx, &numdocs)
            && vbyte_read(idx, &occurs)
            && vbyte_read(idx, &veclen) ))
   {
      info_free(info);
      return NULL;
   }

   /* Initialise the info structure if we weren't already given one */
   if (info == NULL)
   {
      info = malloc(sizeof(struct listinfo));
      info->doc_numbers = malloc(numdocs * sizeof(unsigned long));
      info->phrase_frequency = malloc(numdocs * sizeof(unsigned long));
      info->phrase_offsets = malloc(numdocs * sizeof(unsigned long *));
      memset(info->phrase_offsets, 0, numdocs * sizeof(unsigned long *));
      info->size = numdocs;
   }
   
   /* We may need to enlarge various aspects of the structure if it is not
    * big enough to hold the current postings list */
   if (info->size < numdocs)
   {
      info->doc_numbers = realloc(info->doc_numbers, numdocs * sizeof(unsigned long));
      info->phrase_frequency = realloc(info->phrase_frequency, numdocs * sizeof(unsigned long));
      info->phrase_offsets = realloc(info->phrase_offsets, numdocs * sizeof(unsigned long *));
      memset(info->phrase_offsets + info->size, 0, (numdocs - info->size) * sizeof(unsigned long *));
      info->size = numdocs;
   }

   info->doc_count = numdocs;

   vector = vec_init(veclen);
   vector->len = veclen;

   /* Populate the vector's block of memory from the appropriate location
    * in the vector file */
   fread(vector->vector, veclen, 1, idx);

   /* Read in the stats for each document in the postings list */
   for (c = 0; c < numdocs; c++)
   {
      /* Get the document number (stored as a d-gap) */
      vec_getvbyte(vector, &num);
      info->doc_numbers[c] = last_doc + num + 1;
      last_doc += num + 1;

      /* Read the number of times the phrase occurs in this document */
      vec_getvbyte(vector, &info->phrase_frequency[c]);
      if (info->phrase_offsets[c] != NULL)
         free(info->phrase_offsets[c]);

      /* Prepare the memory for the appropriate number of offsets */
      info->phrase_offsets[c] = malloc(info->phrase_frequency[c] * sizeof(unsigned long));
      memset(info->phrase_offsets[c], 0, info->phrase_frequency[c] * sizeof(unsigned long));
      
      /* Read in all the document offsets */
      last_offset = -1;
      for (d = 0; d < info->phrase_frequency[c]; d++)
      {
         vec_getvbyte(vector, info->phrase_offsets[c] + d);
         last_offset += (info->phrase_offsets[c])[d] + 1;
         (info->phrase_offsets[c])[d] = last_offset;
      }
   }

   vec_free(vector);

   return info;
}

void info_free(struct listinfo *info)
{
   free(info->doc_numbers);
   free(info->phrase_frequency);

   free(info);
}

struct index_scanner *reset_index(struct index_scanner *scn)
{
   fseek(scn->idx_file, 0, SEEK_SET);

   return scn;
}

struct index_scanner *close_index(struct index_scanner *scn)
{
   fclose(scn->idx_file);

   free(scn);

   return scn;
}

