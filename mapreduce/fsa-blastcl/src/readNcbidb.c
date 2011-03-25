#include <blast.h>

int readNcbidb_version, readNcbidb_alphabetType;
uint4 readNcbidb_numberOfSequences, readNcbidb_longestSequenceLength;
uint8 readNcbidb_numberOfLetters;
uint4 *readNcbidb_sequenceOffsets, *readNcbidb_descriptionOffsets;
FILE* readNcbidb_pinFile, *readNcbidb_phrFile, *readNcbidb_psqFile;
char* readNcbidb_pinFilename, *readNcbidb_phrFilename, *readNcbidb_psqFilename;

#define readNcbidb_numCodes 25
unsigned char readNcbidb_codes[readNcbidb_numCodes] =
{31,1,20,18,8,5,14,2,17,10,7,0,16,12,9,13,11,3,6,4,19,22,15,21,23};

// Read an 64-bit integer
uint8 readNcbidb_read64Int(FILE* file, char* filename);
// Read an integer
uint4 readNcbidb_readInt(FILE* file, char* filename);

int main()
{
	unsigned char* sequence, *description;
	uint4 sequenceLength, descriptionLength, sequenceCount = 0;

//	readNcbidb_open("/home/mcam/hadden/data/nr5");
    readNcbidb_open("/home/mcam/hadden/data/oldnrs/nr");

    encoding_initialize(readNcbidb_alphabetType);

    // For each sequence in collection
    while (sequenceCount < readNcbidb_numberOfSequences)
    {
    	// Get sequence, description
    	sequence = readNcbidb_getSequence(sequenceCount, &sequenceLength);
		description = readNcbidb_getDescription(sequenceCount, &descriptionLength);

        // Display in FASTA format
		printf(">%s\n", description);
        print_singleSequence(sequence, sequenceLength); printf("\n");

		free(description);
        free(sequence);

        sequenceCount++;
    }
}

// Open NCBI collection files for reading
void readNcbidb_open(char* filename)
{
    int titleLength, createdLength;
    char* title, *created;
    uint4 sequenceNumber = 0, headerOffset;

    // Open PIN file
    readNcbidb_pinFilename = global_malloc(strlen(filename) + 5);
    sprintf(readNcbidb_pinFilename, "%s.pin", filename);
    if ((readNcbidb_pinFile = fopen(readNcbidb_pinFilename, "r")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for reading\n", readNcbidb_pinFilename);
		exit(-1);
	}

    // Read header information
	readNcbidb_version = readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename);
	readNcbidb_alphabetType = !(readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename));

    // Read title
	titleLength = readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename);
    title = global_malloc(titleLength + 1);
    if (fread(title, sizeof(char), titleLength, readNcbidb_pinFile) < titleLength)
    {
		fprintf(stderr, "Error reading from file %s\n", readNcbidb_pinFilename);
		exit(-1);
    }
    title[titleLength] = '\0';
//	printf("Title = \"%s\"\n", title);

    // Read creation time
	createdLength = readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename);
    created = global_malloc(createdLength + 1);
    if (fread(created, sizeof(char), createdLength, readNcbidb_pinFile) < createdLength)
    {
		fprintf(stderr, "Error reading from file %s\n", readNcbidb_pinFilename);
		exit(-1);
    }
    created[createdLength] = '\0';
//	printf("Created = \"%s\"\n", created);

    // Read more header info
    readNcbidb_numberOfSequences = readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename);
    readNcbidb_numberOfLetters = readNcbidb_read64Int(readNcbidb_pinFile, readNcbidb_pinFilename);
	readNcbidb_longestSequenceLength = readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename);

//    printf("%d,%llu,%d\n", readNcbidb_numberOfSequences, readNcbidb_numberOfLetters, readNcbidb_longestSequenceLength);

    // Read description offsets
    readNcbidb_descriptionOffsets = global_malloc(sizeof(uint4) * (readNcbidb_numberOfSequences + 1));
    sequenceNumber = 0;
    while (sequenceNumber <= readNcbidb_numberOfSequences)
    {
    	readNcbidb_descriptionOffsets[sequenceNumber]
        	= readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename);
        sequenceNumber++;
    }

    // Read sequence offsets
    readNcbidb_sequenceOffsets = global_malloc(sizeof(uint4) * (readNcbidb_numberOfSequences + 1));
    sequenceNumber = 0;
    while (sequenceNumber <= readNcbidb_numberOfSequences)
    {
    	readNcbidb_sequenceOffsets[sequenceNumber]
        	= readNcbidb_readInt(readNcbidb_pinFile, readNcbidb_pinFilename);

//		printf("[%d]", readNcbidb_sequenceOffsets[sequenceNumber]);
        sequenceNumber++;
    }

    fclose(readNcbidb_pinFile);

	// Open PHR and PSQ files for reading
    readNcbidb_phrFilename = global_malloc(strlen(filename) + 5);
    sprintf(readNcbidb_phrFilename, "%s.phr", filename);
    if ((readNcbidb_phrFile = fopen(readNcbidb_phrFilename, "r")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for reading\n", readNcbidb_phrFilename);
		exit(-1);
	}
    readNcbidb_psqFilename = global_malloc(strlen(filename) + 5);
    sprintf(readNcbidb_psqFilename, "%s.psq", filename);
    if ((readNcbidb_psqFile = fopen(readNcbidb_psqFilename, "r")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for reading\n", readNcbidb_psqFilename);
		exit(-1);
	}
}

// Convert a code from NCBI-BLAST database encoding to FSA-BLAST encoding
unsigned char readNcbidb_convertCode(unsigned char code)
{
	if (code >= readNcbidb_numCodes)
    {
    	fprintf(stderr, "Error invalid code %d in NCBI-BLAST database\n", code);
        exit(-1);
    }

	return readNcbidb_codes[code];
}

// Read a sequence from the collection
unsigned char* readNcbidb_getSequence(uint4 sequenceNumber, uint4* length)
{
	char* sequence;
    uint4 offset, position = 0;

    if (sequenceNumber >= readNcbidb_numberOfSequences)
    	return NULL;

	*length = readNcbidb_sequenceOffsets[sequenceNumber + 1]
            - readNcbidb_sequenceOffsets[sequenceNumber] - 1;
	offset = readNcbidb_sequenceOffsets[sequenceNumber];

    sequence = (char*)global_malloc(*length);
    if (fseek(readNcbidb_psqFile, offset, SEEK_SET) == -1 ||
        fread(sequence, sizeof(char), *length, readNcbidb_psqFile) < *length)
    {
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error reading sequence from %s\n", readNcbidb_psqFilename);
		exit(-1);
    }

    while (position < *length)
    {
		sequence[position] = readNcbidb_convertCode(sequence[position]);
    	position++;
    }

    return sequence;
}

// Print the bits in a byte
void readNcbidb_printBits(unsigned char value)
{
	int count = 8;
    while (count > 0)
    {
    	count--;
        printf("%d", (value >> count) & 1);
    }
    printf("\n");
}

// Read sequence length encoding in ASN.1 field length format
unsigned char* readNcbidb_readASN1sequenceLength(unsigned char* data, uint4* length)
{
	int lengthLength;

	if ((*data) & 0x80)
    {
		lengthLength = *data & 0x7F;
    	*length = 0;
        while (lengthLength > 0)
        {
        	lengthLength--;
            data++;
            (*length) <<= 8;
            (*length) |= *data;
        }
    }
    else
    {
    	(*length) = *data & 0x7F;
    }
    data++;

    return data;
}

// Read a description from the collection
unsigned char* readNcbidb_getDescription(uint4 sequenceNumber, uint4* length)
{
	unsigned char* description, *newDescription, *descriptionStart;
    uint4 offset, startPos;

    if (sequenceNumber >= readNcbidb_numberOfSequences)
    	return NULL;

	*length = readNcbidb_descriptionOffsets[sequenceNumber + 1]
            - readNcbidb_descriptionOffsets[sequenceNumber] - 1;
	offset = readNcbidb_descriptionOffsets[sequenceNumber];

//    printf("offset=%d length=%d\n", offset, *length); fflush(stdout);

    description = (char*)global_malloc((*length) + 1);
    if (fseek(readNcbidb_phrFile, offset, SEEK_SET) == -1 ||
        fread(description, sizeof(char), *length, readNcbidb_phrFile) < *length)
    {
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error reading description from %s\n", readNcbidb_phrFilename);
		exit(-1);
    }
    description[*length] = '\0';

	// Read sequence length in ASN.1 format
    descriptionStart = description + 7;
    descriptionStart = readNcbidb_readASN1sequenceLength(descriptionStart, length);

    // Copy description to new buffer
    newDescription = global_malloc((*length) + 1);
    memcpy(newDescription, descriptionStart, (*length) + 1);
    free(description);

    return newDescription;
}

// Read an integer
uint4 readNcbidb_readInt(FILE* file, char* filename)
{
	uint4 value = 0, count = 0;
	unsigned char bytes[4];

    if (fread(bytes, sizeof(char), 4, file) < 4)
    {
		fprintf(stderr, "Error reading from file %s\n", filename);
		exit(-1);
    }

    while (count < 4)
    {
    	value <<= 8;
        value |= bytes[count];
//		printf("[%d:%d]", bytes[count], value);
    	count++;
    }

    return value;
}

// Read an 64-bit integer
uint8 readNcbidb_read64Int(FILE* file, char* filename)
{
	uint8 value = 0;
    int count = 8;
	unsigned char bytes[8];

    if (fread(bytes, sizeof(char), 8, file) < 8)
    {
		fprintf(stderr, "Error reading from file %s\n", filename);
		exit(-1);
    }

    while (count > 0)
    {
    	count--;
    	value <<= 8;
        value |= bytes[count];
//		printf("[%d:%d]", bytes[count], value);
    }

    return value;
}
