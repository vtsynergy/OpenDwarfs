#ifndef CRC_FORMATS_H
#define CRC_FORMATS_H

#include<stdio.h>
#include<stdlib.h>

#include "../../../include/common_args.h"

unsigned int* read_crc(unsigned int* num_pages,unsigned int* page_size,const char* file_path);
void write_crc(const unsigned int** pages, const unsigned int num_pages, const unsigned int page_size,const char* file_path);
unsigned int** rand_crc(const unsigned int num_pages,const unsigned int page_size,const unsigned int seed);
void free_crc(unsigned int** pages, const unsigned int num_pages);


#endif
