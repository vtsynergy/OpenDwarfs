#include "../../../include/common_args.h"
#include "../inc/crc_formats.h"

#include<stdio.h>
#include<stdlib.h>
#include<sys/types.h>
#include<unistd.h>
#include<getopt.h>
#include<time.h>

#define CRC_NAME_MAX_LENGTH 256

static struct option long_options[] = {
	/* name, has_arg, flag, val */
	{"print", 0, NULL, 'p'},
	{"num-pages",1,NULL, 'n'},
	{"size",1,NULL, 's'},
	{"crc-file",1,NULL,'f'},
	{"no-rand",0,NULL,'R'},
	{"no-save",0,NULL,'S'},
	{0,0,0,0}
};

int main(int argc, char** argv)
{
	unsigned int page_size = 1024,num_pages=1,i;
	unsigned long seed=10000;
	char *file_path=NULL,do_print=0,do_rand=1,do_save=1,free_file=0;
	time_t t;
	struct tm tm;
	int opt,option_index=0;

	unsigned int** pages;

	const char* usage = "Usage: %s [-n <num_pages>] [-s <page_size>] [-p] [-f <file_path>] [-R] [-S]\n\n \
			     -n: Create <num_pages> pages - Default is 1\n \
			     -s: Set # of bytes with each page to <page_size> - Default is 1024\n \
			     -f: Save CRC pages to file <file_path> rather than the default: ../test/combinational-logic/crc/crcfile_N<num_pages>_S<page_size>_<year>-<month>-<day>-<hour>-<min>\n \
			     -p: Print data to stdout in decimal format\n \
			     -R: Do NOT seed random number generator in order to produce repeatable results\n \
			     -S: Do NOT save the CRC messages to a file\n\n";

	while ((opt = getopt_long(argc, argv, "n:s:f:pRS", long_options, &option_index)) != -1 )
	{
		switch(opt)
		{
			case 'n':
				if(optarg != NULL)
					num_pages = atoi(optarg);
				else
					num_pages = atoi(argv[optind]);
				break;
			case 's':
				if(optarg != NULL)
					page_size = atoi(optarg);
				else
					page_size = atoi(argv[optind]);
				break;
			case 'f':
				if(optarg != NULL)
					file_path = optarg;
				else
					file_path = argv[optind];
				break;
			case 'p':
				do_print = 1;
				break;
			case 'R':
				do_rand = 0;
				break;
			case 'S':
				do_save = 0;
				break;
			default:
				fprintf(stderr, usage,argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	if(do_rand) seed = (unsigned long) getpid();

	printf("Generating Messages for CRC...\n\n");
	pages = rand_crc(num_pages,page_size,seed);

	//if(do_print) ; //TODO: print data here
	//else
	printf("Data Generated.\n");

	if(do_save)
	{
		if(!file_path)
		{
			file_path = malloc(sizeof(char)*CRC_NAME_MAX_LENGTH);
			check(file_path != NULL,"createcrc.main() - Heap Overflow! Cannot allocate space for 'file_path'");
			t = time(NULL);
			tm = *localtime(&t);
			snprintf(file_path,(sizeof(char)*CRC_NAME_MAX_LENGTH),"../test/combinational-logic/crc/crcfile_N%u_S%u_%d-%d-%d-%d-%d",num_pages,page_size,(tm.tm_year-100),tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min);
			free_file=1;
		}
		printf("Saving CRC Messages to File '%s'...\n\n",file_path);
		write_crc((const unsigned int**)pages,num_pages,page_size,file_path);
	}
	if(free_file) free(file_path);
	free_crc(pages,num_pages);
	return 0;
}
