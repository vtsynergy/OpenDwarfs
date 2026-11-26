#include "../inc/crc_formats.h"

unsigned int* read_crc(unsigned int* num_pages,unsigned int* page_size,const char* file_path)
{
	FILE* fp;
	unsigned int i,j,read_count,num_words;
	unsigned int* page;

	fp = fopen(file_path,"r");
	check(fp != NULL,"crc_formats.read_crc() - Cannot Open File");
	fscanf(fp,"%u\n",num_pages);
	fscanf(fp,"%u\n\n",page_size);
	num_words = *page_size / 4;

	page = int_new_array(sizeof(int)*(*num_pages)*num_words,"crc_formats.read_crc() - Heap Overflow! Cannot allocate space for page");

	for(j=0; j<*num_pages; j++)
	{
		read_count = 0;
		for(i=0; i<num_words; i++)
			read_count += fscanf(fp,"%u ",&page[j*num_words+i]);
		check(read_count == num_words,"crc_formats.read_crc() - Input file corrupted! Read count differs from page size");
		fscanf(fp,"\n");
	}

	fclose(fp);
	return page;
}

void write_crc(const unsigned int** pages, const unsigned int num_pages, const unsigned int page_size,const char* file_path)
{
	FILE* fp;
	unsigned int i,j,num_words;

	num_words = page_size / 4;
	fp = fopen(file_path,"w");
	check(fp != NULL,"crc_formats.write_crc() - Cannot Open File");
	fprintf(fp,"%u\n",num_pages);
	fprintf(fp,"%u\n\n",page_size);

	for(j=0; j<num_pages; j++)
	{
		for(i=0; i<num_words; i++)
			fprintf(fp,"%u ",pages[j][i]);
		fprintf(fp,"\n");
	}

	fclose(fp);
}

unsigned int** rand_crc(const unsigned int num_pages,const unsigned int page_size,const unsigned int seed)
{
	unsigned int i,j,num_words;
	unsigned int** page;

	srand(seed);
	page = malloc(sizeof(void*)*num_pages);
	check(page != NULL,"crc_formats.rand_crc() - Heap Overflow! Cannot allocate space for page*");
	num_words = page_size / 4;
	for(i=0; i<num_pages; i++)
	{
		page[i] = int_new_array(num_words,"crc_formats.rand_crc() - Heap Overflow! Cannot allocate space for pages");
		for(j=0; j<num_words; j++)
			page[i][j] = rand();
	}
	return page;
}

void free_crc(unsigned int** pages, const unsigned int num_pages)
{
	unsigned int j;
	for(j=0; j<num_pages; j++)
	{
		free(pages[j]);
	}
	free(pages);
}
