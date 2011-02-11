// printAccession.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Given a FASTA file and description location (dloc) prints the
// description line including the accession number

#include "blast.h"

int main(int argc, char* argv[])
{
	char* description;
	char file[1000];

    if (argc < 3)
	{
		printf("Useage: printDescription <database> <dloc>\n");
		exit(-1);
	}

	sprintf(file, "%s.descriptions", argv[1]);

	descriptions_open(file);
	description = descriptions_getDescription(atoi(argv[2]), 20);
	printf("%s\n", description);
    descriptions_close();
}
