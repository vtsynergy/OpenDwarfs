// parameters.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Parameters used by BLAST and code for processing command line arguments
// to setting default and non-default values for parameters

#include <errno.h>
#include "blast.h"

#define DEFAULT -1
//Shucai
int parameters_blockNum = 56;
int parameters_threadNum = 128;

char* parameters_queryFile = "";
char* parameters_subjectDatabaseFile = "";
char* parameters_hspindexFilename = "";
int4 parameters_numDisplayAlignments = 500;
int4 parameters_numDisplayTracebacks = 250;
int4 parameters_databaseSize = 0;
int4 parameters_numberOfSequences = 0;
char parameters_outputType = 0;

char parameters_scoringSystem = 0;
char parameters_semiGappedScoring = 1;
char parameters_restrictedInsertionScoring = 1;
char parameters_bytepackedScoring = 1;
char parameters_tableScoring = 0;
char parameters_oneHitTrigger = 0;
char parameters_getDescriptions = 1;
char parameters_ssearch = 0;
char parameters_strands = 3;
char parameters_useIndex = 1;
char parameters_filterEnabled = 1;
char parameters_allClusterMembers = 1;

char* parameters_scoringMatrix = "BLOSUM62";
char* parameters_scoringMatrixPath;
int4 parameters_mismatchScore = -3;
int4 parameters_matchScore = 1;

char parameters_wordSize = DEFAULT;
int4 parameters_T = 11;
int4 parameters_A = 40;
char parameters_overlap = 3;
int2 parameters_startGap = DEFAULT;
int2 parameters_extendGap = DEFAULT;
int2 parameters_openGap = DEFAULT;

int4 parameters_wordTableBytes;
int4 parameters_wordTableLetters;
int4 parameters_wordExtraLetters;
int4 parameters_wordExtraBytes;

double parameters_cutoff = 10;

float parameters_ungappedNormalizedTrigger = DEFAULT;
float parameters_ungappedNormalizedDropoff = DEFAULT;
float parameters_gappedNormalizedDropoff = DEFAULT;
float parameters_gappedFinalNormalizedDropoff = DEFAULT;

float parameters_semiGappedR1 = DEFAULT;
float parameters_semiGappedR2 = DEFAULT;
int4 parameters_semiGappedExtensionN = 10;
int2 parameters_semiGappedStartGap = DEFAULT;
int2 parameters_semiGappedExtendGap = DEFAULT;
int2 parameters_semiGappedOpenGap = DEFAULT;
int4 parameters_semiGappedDropoffIncrease = 0;

int2 parameters_bytepackStartGap = DEFAULT;
int2 parameters_bytepackOpenGap = DEFAULT;
int2 parameters_bytepackExtendGap = DEFAULT;
int2 parameters_bytepackOpen4Gap = DEFAULT;
int2 parameters_bytepackExtend4Gap = DEFAULT;
int4 parameters_bytepackDropoffDecrease = DEFAULT;

int4 parameters_verboseDloc = 0;

struct argument
{
	char flag;
    char type;
    void* address;
    char* description;
    char* defaultText;
};

struct defaultPenalties
{
	char* scoringMatrix;
    int2 startGap;
    int2 extendGap;
    int2 semiGappedDifference;
};

const char parameters_STRING = 0;
const char parameters_INT = 1;
const char parameters_DOUBLE = 2;
const char parameters_FLOAT = 3;
const char parameters_BYTE = 4;
const char parameters_SHORTINT = 5;

struct defaultPenalties parameters_getDefaultPenalties(char* scoringMatrix);
void parameters_printArguments(struct argument* parameters_arguments);

// Given an argument such as -i returns the i. If not a flag arguements
// returns zero. If - is passed, returns space character.
char getFlag(char* argument)
{
	if (argument[0] == '-')
	{
    	if (strlen(argument) == 2)
        {
			return argument[1];
        }
        else
        {
        	return ' ';
        }
	}
	else
	{
		return 0;
	}
}

// Process command line arguments and set parameters value
void parameters_processArguments(int4 argc, char* argv[])
{
	int4 count = 1, value, argumentCount;
	float floatValue;
    double doubleValue;
	char *commandLineArgument, *nextArgument;
	char flag;
    struct argument argument;

    #define NUM_ARGUMENTS 26
    struct argument parameters_arguments[NUM_ARGUMENTS] =
    {
	 {'B', parameters_INT, &parameters_blockNum,
	  "Number of blocks in the kernel", "30"},
	 {'T', parameters_INT, &parameters_threadNum,
	  "Thread number per block", "256"},
     {'d', parameters_STRING, &parameters_subjectDatabaseFile,
      "Database", ""},
     {'i', parameters_STRING, &parameters_queryFile,
      "Query File", ""},
     {'A', parameters_INT, &parameters_A,
      "Multiple Hits window size (protein only)", "40"},
     {'f', parameters_INT, &parameters_T,
      "Threshold for extending hits (protein only)", "11"},
     {'e', parameters_DOUBLE, &parameters_cutoff,
      "Expectation value (E)", "10.0"},
     {'m', parameters_BYTE, &parameters_outputType,
      "alignment view options:\n0 = pairwise,\n7 = XML Blast output,\n8 = tabular", "0"},
     {'y', parameters_FLOAT, &parameters_ungappedNormalizedDropoff,
      "Dropoff (X) for blast extensions in bits", "7.0 for protein, 20.0 for nucleotide"},
     {'P', parameters_BYTE, &parameters_oneHitTrigger,
      "0 for multiple hit, 1 for single hit", "0"},
     {'F', parameters_BYTE, &parameters_filterEnabled,
      "Filter query sequence using DUST/SEG (1=True or 0=False)", "1"},
     {'G', parameters_SHORTINT, &parameters_startGap,
      "Cost to open a gap (Matrix dependant)", "11 for protein, 5 for nucleotide"},
     {'E', parameters_SHORTINT, &parameters_extendGap,
      "Cost to extend a gap (Matrix dependant)", "1 for protein, 2 for nucleotide"},
     {'O', parameters_SHORTINT, &parameters_semiGappedStartGap,
      "Cost to open a gap (Semi-gapped alignment, Matrix dependant)", "7"},
     {'X', parameters_FLOAT, &parameters_gappedNormalizedDropoff,
      "X dropoff value for gapped alignment (in bits)", "15.0 for protein, 30.0 for nucleotide"},
     {'N', parameters_FLOAT, &parameters_ungappedNormalizedTrigger,
      "Number of bits to trigger gapping", "22.0 for protein, 25.0 for nucleotide"},
     {'Z', parameters_FLOAT, &parameters_gappedFinalNormalizedDropoff,
      "X dropoff value for final gapped alignment (in bits)",
      "25.0 for protein, 50.0 for nucleotide"},
     {'M', parameters_STRING, &parameters_scoringMatrix,
      "Matrix", "BLOSUM62"},
     {'v', parameters_INT, &parameters_numDisplayAlignments,
      "Number of database sequences to show one-line descriptions for (V)", "500"},
     {'b', parameters_INT, &parameters_numDisplayTracebacks,
      "Number of database sequences to show alignments for (B)", "250"},
     {'W', parameters_BYTE, &parameters_wordSize,
      "Word size", "3 for protein, 11 for nucleotide"},
     {'z', parameters_INT, &parameters_databaseSize,
      "Effective length of the database (use zero for the real size)", "0"},
     {'q', parameters_INT, &parameters_mismatchScore,
      "Penalty for a nucleotide mismatch", "-3"},
     {'r', parameters_INT, &parameters_matchScore,
      "Reward for a nucleotide match", "1"},
     {'S', parameters_BYTE, &parameters_strands,
      "Query strands to search against database (nucleotide only)\n      3 is both, 1 is top, 2 is bottom", "3"},
     {'a', parameters_BYTE, &parameters_allClusterMembers,
      "Display alignments for all members of a cluster  (1=True or 0=False)", "1"},

/*     {'I', parameters_BYTE, &parameters_useIndex,
      "Use inverted index to perform search if available (1 or 0)", "1"},
     {'h', parameters_STRING, &parameters_hspindexFilename,
      "HSP Index File", ""},
     {'s', parameters_BYTE, &parameters_scoringSystem,
      "Scoring System (0 = Bytepacked + Semi + RI  1 = RI only  2 = Semi only  3 = Original  4 = Semi + RI  5 = Table + RI)", "0"},
     {'R', parameters_FLOAT, &parameters_semiGappedR1,
      "Ratio between semi-gapped/byte-packed and gapped cutoffs", "0.68 for semi-gapped, 0.85 for byte-packed"},
     {'J', parameters_INT, &parameters_semiGappedExtensionN,
      "Semi-gapped alignment insertion frequency", "10"},
     {'D', parameters_BYTE, &parameters_getDescriptions,
      "Get sequence descriptions when printing results (1 or 0)", "1"},
     {'V', parameters_BYTE, &parameters_overlap,
      "Maximum distance between two hits for them to be regarded as overlapping", "3"},
     {'l', parameters_INT, &parameters_verboseDloc,
      "Verbose description location (for debugging purposes)", "0"},
     {'B', parameters_INT, &parameters_bytepackDropoffDecrease,
      "Bytepacked alignment dropoff decrease", "0"},
     {'1', parameters_SHORTINT, &parameters_bytepackStartGap,
      "Bytepacked alignment cost to open a gap", "0"},
     {'2', parameters_SHORTINT, &parameters_bytepackExtendGap,
      "Bytepacked/Table-driven alignment cost to extend a gap", "2"}*/
    };

    while (count < argc)
	{
		commandLineArgument = argv[count];
		flag = getFlag(commandLineArgument);
		if (flag)
		{
			// Get (possible) next argument
			nextArgument = NULL;
			count++;
			if (count < argc)
			{
				nextArgument = argv[count];
			}

            // For each possible inbuild argument, check if flag matches it
            argumentCount = 0;
            while (argumentCount < NUM_ARGUMENTS)
            {
				argument = parameters_arguments[argumentCount];
                if (argument.flag == flag)
                {
                	// Check the next commandline argument is not also a flag
                    if (nextArgument == NULL)// || getFlag(nextArgument))
                    {
                        parameters_printArguments(parameters_arguments);
                        fprintf(stderr, "ERROR: No argument given for %s\n", argument.description);
                        exit(-1);
                    }
                    // If argument type is a string, get pointer
                    if (argument.type == parameters_STRING)
                    {
                    	*((char**)argument.address) = nextArgument;
                    }
                    // If it is an integer, get value
                    if (argument.type == parameters_INT)
                    {
                        value = strtol(nextArgument, &nextArgument, 10);
                        if (nextArgument[0] != '\0')
                        {
                            parameters_printArguments(parameters_arguments);
                            fprintf(stderr, "ERROR: Invalid argument given for %s\n",
                                    argument.description);
                            exit(-1);
                        }
                        *((int4*)argument.address) = value;
                    }
                    // If it is a int2eger, get value
                    if (argument.type == parameters_SHORTINT)
                    {
                        value = strtol(nextArgument, &nextArgument, 10);
                        if (nextArgument[0] != '\0')
                        {
                            parameters_printArguments(parameters_arguments);
                            fprintf(stderr, "ERROR: Invalid argument given for %s\n",
                                    argument.description);
                            exit(-1);
                        }
                        *((int2*)argument.address) = value;
                    }
                    // If it is a byte, get value
                    if (argument.type == parameters_BYTE)
                    {
                        value = strtol(nextArgument, &nextArgument, 10);
                        if (nextArgument[0] != '\0')
                        {
                            parameters_printArguments(parameters_arguments);
                            fprintf(stderr, "ERROR: Invalid argument given for %s\n",
                                    argument.description);
                            exit(-1);
                        }
                        *((char*)argument.address) = value;
                    }
                    // If it is a double, get value
                    if (argument.type == parameters_DOUBLE)
                    {
                        doubleValue = atof(nextArgument);
                        if (doubleValue == 0)
                        {
                           parameters_printArguments(parameters_arguments);
                           fprintf(stderr, "ERROR: Invalid argument given for %s\n",
                                    argument.description);
                            exit(-1);
                        }
                        else
                        {
                            *((double*)argument.address) = doubleValue;
                        }
                    }
                    // If it is a float, get value
                    if (argument.type == parameters_FLOAT)
                    {
                        floatValue = atof(nextArgument);
                        if (floatValue == 0)
                        {
                            parameters_printArguments(parameters_arguments);
                            fprintf(stderr, "ERROR: Invalid argument given for %s\n",
                                    argument.description);
                            exit(-1);
                        }
                        else
                        {
                            *((float*)argument.address) = floatValue;
                        }
                    }

                    break;
                }

                argumentCount++;
            }

            // If we reached of list of arguments and didn't find a match
            // then print error
            if (argumentCount == NUM_ARGUMENTS)
            {
				parameters_printArguments(parameters_arguments);
            	if (flag != ' ')
                {
                	fprintf(stderr, "ERROR: Unrecognised parameter -%c\n", flag);
                }
                exit(-1);
            }
		}
		count++;
	}

	// Check at least query and subject database were specified
	if (parameters_queryFile[0] == '\0')
	{
        parameters_printArguments(parameters_arguments);
		fprintf(stderr, "ERROR: Query File not specified\n");
		exit(-1);
	}
	if (parameters_subjectDatabaseFile[0] == '\0')
	{
        parameters_printArguments(parameters_arguments);
		fprintf(stderr, "ERROR: Subject Database File not specified\n");
		exit(-1);
	}

	if (parameters_numDisplayAlignments < parameters_numDisplayTracebacks)
		parameters_numDisplayAlignments = parameters_numDisplayTracebacks;

    if (parameters_mismatchScore > 0)
    	parameters_mismatchScore = -parameters_mismatchScore;
}

#define parameters_setDefault(parameter, value) \
        if (parameter == DEFAULT) parameter = value;

// Given the alphabet type of the query, use the relevant parameter defaults
void parameters_loadDefaults(unsigned char alphabetType)
{
	struct defaultPenalties defaultPenalties;

    // Decode scoring system parameter
    if (alphabetType == encoding_nucleotide)
    {
    	// Nulceotide case
        if (parameters_scoringSystem == 0)
        {
    		parameters_bytepackedScoring = 1;
            parameters_semiGappedScoring = 0;
            parameters_restrictedInsertionScoring = 1;
        }
        else if (parameters_scoringSystem == 1)
        {
    		parameters_bytepackedScoring = 0;
            parameters_semiGappedScoring = 0;
            parameters_restrictedInsertionScoring = 1;
        }
        else if (parameters_scoringSystem == 2)
        {
    		parameters_bytepackedScoring = 0;
            parameters_semiGappedScoring = 1;
            parameters_restrictedInsertionScoring = 0;
        }
        else if (parameters_scoringSystem == 3)
        {
    		parameters_bytepackedScoring = 0;
            parameters_semiGappedScoring = 0;
            parameters_restrictedInsertionScoring = 0;
        }
        else if (parameters_scoringSystem == 4)
        {
            parameters_bytepackedScoring = 0;
            parameters_semiGappedScoring = 1;
            parameters_restrictedInsertionScoring = 1;
        }
        else if (parameters_scoringSystem == 5)
        {
        	parameters_tableScoring = 1;
    		parameters_bytepackedScoring = 0;
            parameters_semiGappedScoring = 0;
            parameters_restrictedInsertionScoring = 1;
        }
	}
    else
    {
    	// Protein case
        parameters_bytepackedScoring = 0;
        if (parameters_scoringSystem == 0)
        {
            parameters_semiGappedScoring = 1;
            parameters_restrictedInsertionScoring = 1;
        }
        else if (parameters_scoringSystem == 1)
        {
            parameters_semiGappedScoring = 0;
            parameters_restrictedInsertionScoring = 1;
        }
        else if (parameters_scoringSystem == 2)
        {
            parameters_semiGappedScoring = 1;
            parameters_restrictedInsertionScoring = 0;
        }
        else if (parameters_scoringSystem == 3)
        {
            parameters_semiGappedScoring = 0;
            parameters_restrictedInsertionScoring = 0;
        }
        else if (parameters_scoringSystem == 4)
        {
            parameters_semiGappedScoring = 1;
            parameters_restrictedInsertionScoring = 1;
        }
        else if (parameters_scoringSystem == 5)
        {
            parameters_semiGappedScoring = 0;
            parameters_restrictedInsertionScoring = 1;
        }
    }

    // If a nucleotide alphabet
    if (alphabetType == encoding_nucleotide)
    {
    	// Use default nucleotide gap penalties
		defaultPenalties.startGap = 5;
		defaultPenalties.extendGap = 2;
		defaultPenalties.semiGappedDifference = 2;

        // Bytepacked alignment parameters
        parameters_setDefault(parameters_bytepackStartGap, 0);
        parameters_setDefault(parameters_bytepackExtendGap, defaultPenalties.extendGap);
        parameters_bytepackOpenGap = parameters_bytepackStartGap + parameters_bytepackExtendGap;

        parameters_bytepackExtend4Gap = 4 * parameters_bytepackExtendGap;
        parameters_bytepackOpen4Gap = parameters_bytepackStartGap + parameters_bytepackExtend4Gap;

        parameters_setDefault(parameters_bytepackDropoffDecrease, 0);

        // Set default word size
        parameters_setDefault(parameters_wordSize, 11);

        // Set trigger and dropoffs
		parameters_setDefault(parameters_ungappedNormalizedTrigger, 25.0);
        parameters_setDefault(parameters_ungappedNormalizedDropoff, 20.0);
		parameters_setDefault(parameters_gappedNormalizedDropoff, 30.0);
		parameters_setDefault(parameters_gappedFinalNormalizedDropoff, 50.0);

        // Set matrix name for nucleotide search
		parameters_scoringMatrix = (char*)global_malloc(sizeof(char) * 30);
        sprintf(parameters_scoringMatrix, "blastn matrix:%d %d",
                parameters_matchScore, parameters_mismatchScore);

        // Set values for performing hit detection
        if (parameters_wordSize >= 11)
        {
            parameters_wordTableBytes = 2;
            parameters_wordExtraBytes = (parameters_wordSize - 11) / 4;
            parameters_wordExtraLetters = ((parameters_wordSize - 11) % 4) + 3;
        }
       	// Current disabled word size < 11
/*        else if (parameters_wordSize >= 7)
        {
            parameters_wordTableBytes = 1;
            parameters_wordExtraBytes = 0;
            parameters_wordExtraLetters = (parameters_wordSize - 4);
        }*/
        else
        {
            fprintf(stderr, "Error: Word size W=%d is too small\n", parameters_wordSize);
            fflush(stderr);
            exit(-1);
        }

        parameters_wordTableLetters = parameters_wordTableBytes * 4;

    }
    // If a protein alphabet
	else
    {
        // Get default penalties for selected protein scoring matrix
        defaultPenalties = parameters_getDefaultPenalties(parameters_scoringMatrix);

        // Set default word size
        parameters_setDefault(parameters_wordSize, 3);

        parameters_setDefault(parameters_ungappedNormalizedTrigger, 22.0);
	    parameters_setDefault(parameters_ungappedNormalizedDropoff, 7.0);
        parameters_setDefault(parameters_gappedNormalizedDropoff, 15.0);
        parameters_setDefault(parameters_gappedFinalNormalizedDropoff, 25.0);
    }

    // If using byte-packed scoring instead of semi-gapped
    if (parameters_bytepackedScoring)
    {
    	// Set byte-packed default R value
    	parameters_setDefault(parameters_semiGappedR1, 0.85);
    	parameters_setDefault(parameters_semiGappedR2, 1.2);
	}
    else if (parameters_tableScoring)
    {
    	// Otherwise if using table driven scoring
    	parameters_setDefault(parameters_semiGappedR1, 1.0);
    	parameters_setDefault(parameters_semiGappedR2, 1.5);
    }
    else if (parameters_semiGappedScoring)
    {
    	// Otherwise set semi-gapped default R value
    	parameters_setDefault(parameters_semiGappedR1, 0.68);
    	parameters_setDefault(parameters_semiGappedR2, 1.2);
    }
    else
    {
    	// Regular gapped scoring
    	parameters_setDefault(parameters_semiGappedR1, 1.0);
    	parameters_setDefault(parameters_semiGappedR2, 1.0);
    }

    // If no open gap value given at the command line, use default
    parameters_setDefault(parameters_startGap, defaultPenalties.startGap);

    // If no extend gap value given at the command line, use default
    parameters_setDefault(parameters_extendGap, defaultPenalties.extendGap);

    // Invert sign of gap existance penalty and extend gap penalty if either is negative
    if (parameters_startGap < 0)
    	parameters_startGap = -parameters_startGap;
    if (parameters_extendGap < 0)
    	parameters_extendGap = -parameters_extendGap;

    // Calculate openGap from startGap
    parameters_openGap = parameters_startGap + parameters_extendGap;

    // If no semi-gapped open gap value given at the command line, use default difference
    parameters_setDefault(parameters_semiGappedStartGap,
    	parameters_startGap - defaultPenalties.semiGappedDifference);

    // Use the same gap extend value for semi-gapped alignment
    parameters_semiGappedExtendGap = parameters_extendGap;

    // Invert sign of semi-gapped gap existance penalty if negative
    if (parameters_semiGappedStartGap < 0)
    	parameters_semiGappedStartGap = -parameters_semiGappedStartGap;

    // Calculate semiGappedOpenGap from semiGappedStartGap
    parameters_semiGappedOpenGap = parameters_semiGappedStartGap + parameters_semiGappedExtendGap;

    // If not CAFE nor nucleotide search
    #ifndef CAFEMODE
    if (alphabetType != encoding_nucleotide)
    {
        parameters_findScoringMatrix();
	}
    #endif
}

// Load a BLOSUM or PAM scoring matrix
void parameters_findScoringMatrix()
{
    FILE* matrixFile, *ncbircFile;
    char* homeDirectory, *ncbircFilename;

    // Get home user's directory
    homeDirectory = getenv("HOME");

    // Construct name of .ncbirc file
    ncbircFilename = (char*)global_malloc(sizeof(char) * (strlen(homeDirectory) + 9));
    sprintf(ncbircFilename, "%s/.ncbirc", homeDirectory);

    // Check for existence of NCBI file
	if ((ncbircFile = fopen(ncbircFilename, "r")) != NULL)
    {
        // Determine the location of the scoring matrix file by consulting the .ncbirc file
        parameters_scoringMatrixPath = (char*)global_malloc(sizeof(char) * 1024);

        if (!(fscanf(ncbircFile, "[NCBI]\nData=%s", parameters_scoringMatrixPath)))
        {
            fprintf(stderr, "Error reading scoring matrix path from %s file\n", ncbircFilename);
            fprintf(stderr, "BLAST requires the file .ncbirc in the user's home directory\n");
            fprintf(stderr, "containing the text:\n\n");
            fprintf(stderr, "[NCBI]\n");
            fprintf(stderr, "Data=/home/user/fsa-blast/data\n\n");
            fprintf(stderr, "Where the path specified contains scoring matrix files (ie. BLOSUM62)\n");
            exit(-1);
        }
    	fclose(ncbircFile);
    }
    else
    {
    	// If not available guess location of scoring matrix files
        parameters_scoringMatrixPath = (char*)global_malloc(sizeof(char) * 1024);
        strcpy(parameters_scoringMatrixPath, "data");
    }

    // Append scoring matrix filename then check for existance
    sprintf(parameters_scoringMatrixPath, "%s/%s", parameters_scoringMatrixPath,
            parameters_scoringMatrix);
    matrixFile = fopen(parameters_scoringMatrixPath, "r");
    if (matrixFile == NULL)
    {
        fprintf(stderr, "%s\n", strerror(errno));
        fprintf(stderr, "Error reading matrix file %s\n", parameters_scoringMatrixPath);
        fprintf(stderr, "BLAST requires the file .ncbirc in the user's home directory\n");
        fprintf(stderr, "containing the text:\n\n");
        fprintf(stderr, "[NCBI]\n");
        fprintf(stderr, "Data=/home/user/fsa-blast/data\n\n");
        fprintf(stderr, "Where the path specified contains scoring matrix files (ie. BLOSUM62)\n");
        exit(-1);
    }
    fclose(matrixFile);
    free(ncbircFilename);
}

void parameters_printArguments(struct argument* parameters_arguments)
{
	int4 count = 0;
	struct argument argument;

    while (count < NUM_ARGUMENTS)
    {
    	argument = parameters_arguments[count];
		printf("  -%c  %s \n", argument.flag, argument.description);

        printf("    default = %s\n", argument.defaultText);

        count++;
    }
    printf("\n");
}


#define NUM_DEFAULT_PENALTIES 6
static const struct defaultPenalties parameters_defaultPenalties[NUM_DEFAULT_PENALTIES] =
{
	{"BLOSUM45", 14, 2, 4},
	{"BLOSUM50", 13, 2, 4},
	{"BLOSUM62", 11, 1, 4},
	{"BLOSUM80", 10, 1, 5},
	{"PAM30", 10, 1, 7},
	{"PAM70", 9, 1, 5},
};

// Lookup default gap penalties for the scoring matrix being used
struct defaultPenalties parameters_getDefaultPenalties(char* scoringMatrix)
{
	int4 count = 0;

    // Lookup up each entry in the table
    while (count < NUM_DEFAULT_PENALTIES)
    {
    	// If a match is found, use it
		if (strcmp(parameters_defaultPenalties[count].scoringMatrix, scoringMatrix) == 0)
        {
        	return parameters_defaultPenalties[count];
        }
    	count++;
    }

    // No matching entries found, error
    fprintf(stderr, "Error: Matrix %s not supported\n", scoringMatrix);
    exit(-1);
}

void parameters_free()
{
	free(parameters_scoringMatrixPath);
	if (encoding_alphabetType == encoding_nucleotide)
		free(parameters_scoringMatrix);
}
