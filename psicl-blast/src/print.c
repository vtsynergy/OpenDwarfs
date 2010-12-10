// print.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Functions for displaying final alignments to the user in standard BLAST format

#include "blast.h"

// Prototypes
void print_gappedExtension(struct gappedExtension* gappedExtension, struct PSSMatrix PSSMatrix,
                           char* query, unsigned char* subject);
void print_ungappedExtension(struct ungappedExtension* ungappedExtension, char* query,
                             struct PSSMatrix PSSMatrix, unsigned char* subjectSequence);
char* print_cafeDescription(int4 sequenceNumber, int4 startChar);
// Print a gapped extension using XML output
void print_XMLgappedExtension(struct gappedExtension* gappedExtension, struct PSSMatrix PSSMatrix,
                              char* query, unsigned char* subject, uint4 hspNum);
// Print a gapped extension using tabular output
void print_tabularGappedExtension(struct gappedExtension* gappedExtension, struct PSSMatrix PSSMatrix,
                                  char* query, unsigned char* subject, char* queryDescription,
                                  char* subjectDescription);

// Construct the pairwise alignment
void print_constructAlignment(char* queryLine, char* subjectLine, char* midLine,
                              uint4 *identities, uint4 *positives, uint4* gaps, uint4* gapopens,
                              char reverseComplement, struct trace trace, char* query,
                              unsigned char* subject, struct PSSMatrix PSSMatrix)
{
	int4 count = 0;
	int4 inInsertion = 0, charCount = 0, charPosition;
	int4 queryPosition, subjectPosition;
	unsigned char* traceCodes;

    if (reverseComplement)
	{
        charPosition = trace.length;
        queryPosition = trace.queryStart;
        subjectPosition = trace.subjectStart;
	}
    else
    {
		charPosition = 0;
        queryPosition = trace.queryStart;
        subjectPosition = trace.subjectStart;
	}
	traceCodes = trace.traceCodes;

    // For each position in the trace
	while (count < trace.length)
	{
		// A match
		if (traceCodes[count] == 0)
		{
			inInsertion = 0;

			// Add query and subject characters
			queryLine[charPosition] = query[queryPosition];

            if (reverseComplement)
				subjectLine[charPosition] = encoding_getLetter(encoding_getComplement(subject[subjectPosition]));
            else
				subjectLine[charPosition] = encoding_getLetter(subject[subjectPosition]);

			// If chars match, using char for midline
			if (queryLine[charPosition] == subjectLine[charPosition])
			{
				(*identities)++;
				(*positives)++;
            	if (encoding_alphabetType == encoding_nucleotide)
					midLine[charPosition] = '|';
                else
					midLine[charPosition] = queryLine[charPosition];
			}
			// If positive add a plus symbol
			else if (encoding_alphabetType == encoding_protein &&
                     PSSMatrix.matrix[queryPosition][subject[subjectPosition]] > 0)
			{
				(*positives)++;
                midLine[charPosition] = '+';
			}
			// If 0 or negative, blank space
			else
			{
				midLine[charPosition] = ' ';
			}

            if (reverseComplement)
                queryPosition--;
            else
                queryPosition++;

            subjectPosition++;
        }
		// An insertion wrt. query
		else if (traceCodes[count] == 1)
		{
        	gaps++;
            if (!inInsertion)
            	(*gapopens)++;

			queryLine[charPosition] = '-';
			midLine[charPosition] = ' ';
            if (reverseComplement)
				subjectLine[charPosition] = encoding_getLetter(encoding_getComplement(subject[subjectPosition]));
			else
				subjectLine[charPosition] = encoding_getLetter(subject[subjectPosition]);

			inInsertion = 1;

            subjectPosition++;
		}
		// An insertion wrt. subject
		else
		{
        	gaps++;
            if (!inInsertion)
            	(*gapopens)++;

            subjectLine[charPosition] = '-';
			midLine[charPosition] = ' ';
			queryLine[charPosition] = query[queryPosition];

			inInsertion = 1;

            if (reverseComplement)
                queryPosition--;
            else
                queryPosition++;
		}

        // Display nucleotide sequences in lower case
        if (encoding_alphabetType == encoding_nucleotide)
        {
            subjectLine[charPosition] = tolower(subjectLine[charPosition]);
            queryLine[charPosition] = tolower(queryLine[charPosition]);
        }

        if (reverseComplement)
        	charPosition--;
        else
        	charPosition++;

		count++;
        charCount++;
	}

	// Alignment must end with a match
	queryLine[charPosition] = query[queryPosition];
    if (reverseComplement)
        subjectLine[charPosition] = encoding_getLetter(encoding_getComplement(subject[subjectPosition]));
    else
        subjectLine[charPosition] = encoding_getLetter(subject[subjectPosition]);

    // Display nucleotide sequences in lower case
    if (encoding_alphabetType == encoding_nucleotide)
    {
        subjectLine[charPosition] = tolower(subjectLine[charPosition]);
        queryLine[charPosition] = tolower(queryLine[charPosition]);
    }

    // If chars match, using char for midline
	if (queryLine[charPosition] == subjectLine[charPosition])
	{
		(*identities)++;
		(*positives)++;
        if (encoding_alphabetType == encoding_nucleotide)
            midLine[charPosition] = '|';
        else
            midLine[charPosition] = queryLine[charPosition];
	}
	// If positive add a plus symbol
	else if (PSSMatrix.matrix[queryPosition][subject[subjectPosition]] > 0)
	{
		(*positives)++;
		midLine[charPosition] = '+';
	}
	// If 0 or negative, blank space
	else
	{
		midLine[charPosition] = ' ';
	}
	charCount++;

	// Null-terminate each of the strings
	subjectLine[charCount] = '\0';
	queryLine[charCount] = '\0';
	midLine[charCount] = '\0';
}

// Format a sequence description so it wraps around lines nicely
char* print_formatDescription(char* description, int4 firstLineIndent, int4 remainingLinesIndent,
                              int4 maxLineLength)
{
	int4 length, startWord = 0, endWord = 0, startLine = 0, newDescriptionPos = 0;
    int4 maxLengthAfterIndent, count;
    char* newDescription;

    length = strlen(description);

    // Declare memory for formatted description, which will be at most double its length
    newDescription = (char*)global_malloc(sizeof(char) * length * 2);

    startLine = startWord = endWord = 0;

    while (endWord < length)
	{
        // Find the end of the current word
        while (endWord < length && description[endWord] != ' ' &&
               description[endWord] != '\n' && description[endWord] != '\t')
        {
            endWord++;
        }

        if (startLine == 0)
        	maxLengthAfterIndent = maxLineLength - firstLineIndent;
		else
        	maxLengthAfterIndent = maxLineLength - remainingLinesIndent;

        // If it is too far along
        if (endWord - startLine > maxLengthAfterIndent)
        {
            // If this is the first word of this line, we need to cut the word in half
            if (startLine == startWord)
            {
                // Copy portion of the word
                memcpy(newDescription + newDescriptionPos, description + startLine,
                       sizeof(char) * maxLineLength);
                newDescriptionPos += maxLineLength;

                // Add carriage return
                newDescription[newDescriptionPos] = '\n';
                newDescriptionPos++;

                // Add indent
                count = 0;
                while (count < remainingLinesIndent)
                {
                    newDescription[newDescriptionPos] = ' ';
                    newDescriptionPos++;
                    count++;
                }

                // Set pointers for next line
                startLine += maxLineLength;
                endWord = startWord = startLine;
            }
            // Else exclude this last word and copy the rest of the line over
            else
            {
            	// Copy line over
				memcpy(newDescription + newDescriptionPos, description + startLine,
                       sizeof(char) * (startWord - startLine));
                newDescriptionPos += (startWord - startLine);

                // Add carriage return
                newDescription[newDescriptionPos] = '\n';
                newDescriptionPos++;

                // Add indent
                count = 0;
                while (count < remainingLinesIndent)
                {
                    newDescription[newDescriptionPos] = ' ';
                    newDescriptionPos++;
                    count++;
                }

                // Set pointers for next line
                startLine = endWord = startWord;
            }
        }
        else
        {
        	// Advance past current whitespace
            endWord++;

            // This is the start of the next word
            startWord = endWord;
        }
	}

    // Copy the last bit over
    memcpy(newDescription + newDescriptionPos, description + startLine,
           sizeof(char) * (length - startLine));
    newDescriptionPos += (length - startLine);

    // Null terminate
    newDescription[newDescriptionPos] = '\0';

    // Free the unformatted description
    free(description);

    return newDescription;
}

// Return new version of description that ends at first whitespace
char* print_untilWhitespace(char* description)
{
	char* newDescription;
	uint4 count = 0, length;
    length = strlen(description);

    while (count < length)
    {
    	if (description[count] == ' ')
			break;

        count++;
    }

    newDescription = (char*)global_malloc(count + 1);
    memcpy(newDescription, description, count);
	newDescription[count] = '\0';

    free(description);

    return newDescription;
}

// Format a sequence description so it wraps around lines nicely
char* print_encodeGreaterLessThan(char* description)
{
	uint4 oldPos = 0, newPos = 0, length;
    char* newDescription;
    length = strlen(description);

    // Declare memory for formatted description, which will be at most quadruple its length
    newDescription = (char*)global_malloc(sizeof(char) * length * 4);

    while (oldPos <= length)
    {
    	// Replace > and < symbols with &gt; and &lt;
		if (description[oldPos] == '>')
        {
        	newDescription[newPos] = '&'; newPos++;
        	newDescription[newPos] = 'g'; newPos++;
        	newDescription[newPos] = 't'; newPos++;
        	newDescription[newPos] = ';'; newPos++;
        }
		else if (description[oldPos] == '<')
        {
        	newDescription[newPos] = '&'; newPos++;
        	newDescription[newPos] = 'l'; newPos++;
        	newDescription[newPos] = 't'; newPos++;
        	newDescription[newPos] = ';'; newPos++;
        }
        else
        {
        	newDescription[newPos] = description[oldPos];
        	newPos++;
        }
    	oldPos++;
    }

    // Free the unformatted description
    free(description);

    return newDescription;
}

// Print full gapped alignments
void print_gappedAlignmentsFull(char* query, struct PSSMatrix PSSMatrix)
{
	struct alignment* alignment;
    struct finalAlignment* finalAlignment;
	struct gappedExtension* currentExtension;
	char* description, *queryDescription;
    unsigned char* subject;
	int4 count = 0, descriptionLength, hspNum;
    uint4 clusterSizes[alignments_numClusters];
    uint4 queryDescriptionLength;

    queryDescriptionLength = strlen(blast_queryDescription);
    queryDescription = (char*)global_malloc(queryDescriptionLength + 1);
    strcpy(queryDescription, blast_queryDescription);
    queryDescription = print_untilWhitespace(queryDescription);

    // Clear array of cluster sizes
    count = 0;
    while (count < alignments_numClusters)
    {
    	clusterSizes[count] = 0;
    	count++;
    }

	// For each alignment we are to print the traceback
    count = 0;
	while ((count < parameters_numDisplayTracebacks || parameters_numDisplayTracebacks == 0)
         && count < alignments_finalAlignments->numEntries)
	{
    	finalAlignment = memSingleBlock_getEntry(alignments_finalAlignments, count);
		alignment = finalAlignment->alignment;
		if (alignment->cluster > 0)
        {
			clusterSizes[alignment->cluster]++;
        }
        count++;
	}

    // For each alignment we are to print the traceback
    count = 0;
	while ((count < parameters_numDisplayTracebacks || parameters_numDisplayTracebacks == 0)
         && count < alignments_finalAlignments->numEntries)
	{
    	finalAlignment = memSingleBlock_getEntry(alignments_finalAlignments, count);
		alignment = finalAlignment->alignment;

        // If this alignment is to be displayed
        if (parameters_allClusterMembers || alignment->cluster == 0 ||
            clusterSizes[alignment->cluster] > 0)
		{
            // Print description
            if (parameters_getDescriptions)
            {
                // If we are in cafe mode, use cafe function to get the description
                #ifdef CAFEMODE
                description = print_cafeDescription(alignment->descriptionLocation, 0);
                #else
                descriptionLength = strlen(finalAlignment->description);
                description = (char*)global_malloc(sizeof(char) * (descriptionLength + 1));
                strcpy(description, finalAlignment->description);
                #endif

                if (parameters_outputType == parameters_xml)
                {
					description = print_encodeGreaterLessThan(description);
				}
                else if (parameters_outputType == parameters_tabular)
                {
                	description = print_untilWhitespace(description);
                }
                else
                {
                	description = print_formatDescription(description, 0, 11, 68);
                }
            }
            else
            {
                description = (char*)global_malloc(sizeof(char));
                strcpy(description, "");
            }

            if (parameters_outputType == parameters_xml)
            {
                printf("        <Hit>\n");
                printf("          <Hit_num>%d</Hit_num>\n", count + 1);
                printf("          <Hit_id>gnl|BL_ORD_ID|%d</Hit_id>\n", alignment->descriptionLocation);
                printf("          <Hit_def>%s</Hit_def>\n", description);
                printf("          <Hit_accession>%d</Hit_accession>\n", alignment->descriptionLocation);
                printf("          <Hit_len>%d</Hit_len>\n", alignment->subjectLength);
                printf("          <Hit_hsps>\n");
            }
            else if (parameters_outputType == parameters_tabular)
            {
            }
			else
            {
                printf(">%s\n", description);
                if (!parameters_allClusterMembers && clusterSizes[alignment->cluster] > 1)
                {
                    printf("          [%d near-identical alignment(s) not displayed]\n",
                           clusterSizes[alignment->cluster] - 1);
                    clusterSizes[alignment->cluster] = 0;
                }
//                printf("          Length = %d DescriptionLocation = %d\n\n",
//                       alignment->subjectLength, alignment->descriptionLocation);
        		printf("          Length = %d\n\n", alignment->subjectLength);
			}

            // Get list of gapped extensions
            currentExtension = alignment->gappedExtensions;
			hspNum = 0;

            // For each gapped extension
            while (currentExtension != NULL)
            {
            	hspNum++;

                if (parameters_ssearch)
                    subject = alignment->subject;
                else
                    // Get the unpacked subject
                    subject = unpack_selectRegion(alignment->unpackRegions, alignment->numUnpackRegions,
                                                  currentExtension->subjectEnd)->unpackedSubject;

                // Print the gapped extension
                if (parameters_outputType == parameters_xml)
                    print_XMLgappedExtension(currentExtension, PSSMatrix, query, subject, hspNum);
                else if (parameters_outputType == parameters_tabular)
                    print_tabularGappedExtension(currentExtension, PSSMatrix, query, subject,
                                                 queryDescription, description);
                else if (parameters_outputType != parameters_tabular)
                	print_gappedExtension(currentExtension, PSSMatrix, query, subject);

                currentExtension = currentExtension->next;
            }

            if (parameters_outputType == parameters_xml)
            {
                printf("          </Hit_hsps>\n");
                printf("        </Hit>\n");
			}

            free(description);
        }

		count++;
	}

    free(queryDescription);
}

// Print 1 line description of gapped alignments
void print_gappedAlignmentsBrief()
{
	char* line;
	int4 descriptionLength;
	struct alignment* currentAlignment;
    struct finalAlignment* finalAlignment;
	struct gappedExtension* currentExtension;
	int4 count = 0;
	char* description;
    double normalizedScore, eValue;

	printf("                                                                 Score    E\n");
	printf("Sequences producing significant alignments:                      (bits) Value\n\n");

	// For each alignment
	while (count < alignments_finalAlignments->numEntries
        && count < parameters_numDisplayAlignments)
	{
    	finalAlignment = memSingleBlock_getEntry(alignments_finalAlignments, count);
		currentAlignment = finalAlignment->alignment;
		currentExtension = currentAlignment->gappedExtensions;

        // Get description of subject
		if (parameters_getDescriptions)
        {
            #ifdef CAFEMODE
            description = print_cafeDescription(currentAlignment->descriptionLocation, 1);
            descriptionLength = strlen(description);
            #else
            descriptionLength = strlen(finalAlignment->description);
        	description = (char*)global_malloc(sizeof(char) * (descriptionLength + 1));
            strcpy(description, finalAlignment->description);
			#endif
        }
        else
        {
        	description = (char*)global_malloc(sizeof(char)); strcpy(description, "");
            descriptionLength = 0;
		}

        // Realloc description string to length displayed
        description = (char*)global_realloc(description, sizeof(char) * 67);

        // Construct line with description, normalized score and eValue
        line = (char*)global_malloc(sizeof(char) * 90);

        // If description is more than 63 chars long
        if (descriptionLength > 63)
        {
            // Cut off end and add ...
            description[63] = '\0';
            strcat(description, "...");
        }
        else
        {
            // Otherwise pad out remaining with spaces
            while (descriptionLength < 66)
            {
                description[descriptionLength] = ' ';
                descriptionLength++;
            }
            description[descriptionLength] = '\0';
        }

        // If no gapped extensions have been performed and scored, use the alignment's
        // highest nominal score to calculate normalized score and e-value
		if (currentExtension == NULL)
		{
        	normalizedScore = statistics_gappedNominal2normalized(
                              finalAlignment->highestNominalScore);
            eValue = statistics_gappedCalculateEvalue(normalizedScore);

            sprintf(line, "%s  %4.0f  %s\n", description, normalizedScore,
			        print_eValue2String(eValue));
        }
        else
        {
			sprintf(line, "%s  %4.0f  %s\n", description, currentExtension->normalizedScore,
			        print_eValue2String(currentExtension->eValue));
		}

        printf("%s", line);

        // Free memory used constructing line
        free(line);
        free(description);

		count++;
	}
	printf("\n");
}

// Print a gapped extension
void print_gappedExtension(struct gappedExtension* gappedExtension, struct PSSMatrix PSSMatrix,
                           char* query, unsigned char* subject)
{
	int4 queryPosition, subjectPosition;
	int4 count = 0, charCount = 0, length, lineStart;
	char *queryLine, *subjectLine, *midLine;
	int4 identities = 0, positives = 0, gaps = 0, gapopens = 0;
	int4 numberOfSections;
	char* pairwiseAlignment;
	char* finalText;
	char temp[71], temp2[20];
	struct trace trace;
	int4 longestNumber = 0;
    char reverseComplement = 0;

	trace = gappedExtension->trace;

    if (trace.queryStart >= PSSMatrix.strandLength)
    {
    	trace.queryStart = PSSMatrix.length - trace.queryStart - 1;
    	reverseComplement = 1;
    }

    // Determine the length of the longest query/subject position number
	if ((log(trace.queryStart + 1) / log(10)) > longestNumber)
    {
    	longestNumber = (log(trace.queryStart) / log(10));
    }
	if ((log(trace.subjectStart + 1) / log(10)) > longestNumber)
    {
    	longestNumber = (log(trace.subjectStart) / log(10));
    }
	if ((log(gappedExtension->queryEnd + 1) / log(10)) > longestNumber)
    {
    	longestNumber = (log(gappedExtension->queryEnd) / log(10));
    }
	if ((log(gappedExtension->subjectEnd + 1) / log(10)) > longestNumber)
    {
    	longestNumber = (log(gappedExtension->subjectEnd) / log(10));
    }
    longestNumber++;

    // Declare memory for query, subject and midlines
	queryLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	subjectLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	midLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));

	queryLine[0] = '\0'; subjectLine[0] = '\0'; midLine[0] = '\0';

    print_constructAlignment(queryLine, subjectLine, midLine, &identities, &positives, &gaps,
                             &gapopens, reverseComplement, trace, query, subject, PSSMatrix);

	// Allocate memory for building final pairwise alignment
	length = trace.length + 1;
	numberOfSections = ((length - 1) / 60) + 1;
	pairwiseAlignment = (char*)global_malloc(sizeof(char) * (numberOfSections * 250));
    pairwiseAlignment[0] = '\0';

    if (reverseComplement)
	{
        queryPosition = PSSMatrix.length - gappedExtension->queryEnd - 1;
        subjectPosition = gappedExtension->subjectEnd;
    }
    else
    {
        queryPosition = trace.queryStart;
        subjectPosition = trace.subjectStart;
	}

	// Until we reach end of the alignment
	charCount = 0;
	while (charCount < length)
	{
		// Remember where the line started in the strings
		lineStart = charCount;

		// First build the query line. Example:
		// Query:    4 GDESERIVIN-----HQTYRSTLRTLPGTRLAWLAEPDAHSH-----EDYDPRADEFFFD 58
		strcat(pairwiseAlignment, "Query: ");
        sprintf(temp2, "%%-%dd ", longestNumber);
		sprintf(temp, temp2, queryPosition + 1);
		strcat(pairwiseAlignment, temp);

		count = 0;
		while (charCount < length && count < 60)
		{
			temp[count] = queryLine[charCount];
			if (temp[count] != '-')
			{
                queryPosition++;
			}
			else
			{
				gaps++;
			}

			charCount++;
			count++;
		}

		temp[count] = '\0';
		strcat(pairwiseAlignment, temp);

		sprintf(temp, " %d\n", queryPosition);
		strcat(pairwiseAlignment, temp);

        // Second build the middle line. Example:
		//             G+ S R+++NVGG +H+    TL  +P TRL  L + + H        DY    +E+FFD
		charCount = lineStart;
        sprintf(temp2, "       %%%ds ", longestNumber);
        sprintf(temp, temp2, "");

        strcat(pairwiseAlignment, temp);

		count = 0;
		while (charCount < length && count < 60)
		{
			temp[count] = midLine[charCount];

			charCount++;
			count++;
		}

		temp[count] = '\n';
		temp[count + 1] = '\0';
		strcat(pairwiseAlignment, temp);

		// Finally build the subject line. Example:
		// Sbjct: 2194 GEVSRRVILNVGGVKHEVLWRTLDRVPHTRLGKLKDCNTHDAIVDLCDDYSLAENEYFFD 2784
		charCount = lineStart;
		strcat(pairwiseAlignment, "Sbjct: ");
        sprintf(temp2, "%%-%dd ", longestNumber);
		sprintf(temp, temp2, subjectPosition + 1);
		strcat(pairwiseAlignment, temp);

		count = 0;
		while (charCount < length && count < 60)
		{
			temp[count] = subjectLine[charCount];
			if (temp[count] != '-')
			{
                if (reverseComplement)
                    subjectPosition--;
                else
                    subjectPosition++;
			}
			else
			{
				gaps++;
			}

			charCount++;
			count++;
		}

		temp[count] = '\0';
		strcat(pairwiseAlignment, temp);

        if (reverseComplement)
            sprintf(temp, " %d\n\n", subjectPosition + 2);
		else
            sprintf(temp, " %d\n\n", subjectPosition);

        strcat(pairwiseAlignment, temp);
	}

	// Construct the final text
	finalText = (char*)global_malloc(numberOfSections * 300 + 300);

    sprintf(finalText,
	" Score = %.1f bits (%d), Expect = %s\n Identities = %d/%d (%d%%)",
		gappedExtension->normalizedScore, gappedExtension->nominalScore,
		print_eValue2String(gappedExtension->eValue), identities, length,
		(int4)(identities * 100 / length));

    if (encoding_alphabetType == encoding_protein)
    {
    	sprintf(temp, ", Positives = %d/%d (%d%%)", positives, length, (int4)(positives * 100 / length));
    	strcat(finalText, temp);
    }

    if (gaps > 0)
    {
    	sprintf(temp, ", Gaps = %d/%d (%d%%)", gaps, length, (int4)(gaps * 100 / length));
    	strcat(finalText, temp);
    }

	if (reverseComplement)
    {
    	strcat(finalText, "\n Strand = Plus / Minus\n\n");
    }
    else
    {
    	strcat(finalText, "\n Strand = Plus / Plus\n\n");
    }

    strcat(finalText, pairwiseAlignment);
	free(pairwiseAlignment);

	printf("%s\n", finalText);
	free(finalText);
	free(queryLine);
	free(subjectLine);
	free(midLine);
}

// Convert an e-value to a string representation, using the same
// technique as NCBI-BLAST
char print_eValueBuf[10];
char* print_eValue2String(double evalue)
{
	char* eValueBuffer = print_eValueBuf;

	if (evalue < 1.0e-180) {
        sprintf(print_eValueBuf, "0.0");
    } else if (evalue < 1.0e-99) {
        sprintf(print_eValueBuf, "%2.0le", evalue);
        eValueBuffer++;	/* Knock off digit. */
    } else if (evalue < 0.0009) {
        sprintf(print_eValueBuf, "%3.0le", evalue);
    } else if (evalue < 0.1) {
        sprintf(print_eValueBuf, "%4.3lf", evalue);
    } else if (evalue < 1.0) {
        sprintf(print_eValueBuf, "%3.2lf", evalue);
    } else if (evalue < 10.0) {
        sprintf(print_eValueBuf, "%2.1lf", evalue);
    } else {
        sprintf(print_eValueBuf, "%5.0lf", evalue);
    }

	return eValueBuffer;
}

// Print a single sequence
void print_singleSequence(unsigned char* sequence, int4 length)
{
	int4 count = 0;
    while (count < length)
    {
		printf("%c", encoding_getLetter(sequence[count]));
    	count++;
    }
}

// Header for XML output
void print_XMLheader(char* query, struct PSSMatrix PSSMatrix)
{
	printf("<?xml version=\"1.0\"?>\n");
    printf("<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"NCBI_BlastOutput.dtd\">\n");
    printf("<BlastOutput>\n");
    if (encoding_alphabetType == encoding_protein)
		printf("  <BlastOutput_program>blastp</BlastOutput_program>\n");
	else
		printf("  <BlastOutput_program>blastn</BlastOutput_program>\n");
	printf("  <BlastOutput_version>FSA-BLAST 1.03</BlastOutput_version>\n");
    printf("  <BlastOutput_reference>~Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, ~Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), ~&quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search~programs&quot;,  Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>\n");
	printf("  <BlastOutput_db>%s</BlastOutput_db>\n", parameters_subjectDatabaseFile);
    printf("  <BlastOutput_query-ID>lcl|QUERY</BlastOutput_query-ID>\n");
    printf("  <BlastOutput_query-def>%s</BlastOutput_query-def>\n", blast_queryDescription);
    printf("  <BlastOutput_query-len>%d</BlastOutput_query-len>\n", PSSMatrix.length);
    printf("  <BlastOutput_param>\n");
    printf("    <Parameters>\n");
    printf("      <Parameters_matrix>%s</Parameters_matrix>\n", parameters_scoringMatrix);
    printf("      <Parameters_expect>%g</Parameters_expect>\n", parameters_cutoff);
    printf("      <Parameters_gap-open>%d</Parameters_gap-open>\n", parameters_startGap);
    printf("      <Parameters_gap-extend>%d</Parameters_gap-extend>\n", parameters_extendGap);
	if (parameters_filterEnabled)
    {
    	if (encoding_alphabetType == encoding_protein)
        	printf("      <Parameters_filter>S</Parameters_filter>\n");
		else
        	printf("      <Parameters_filter>D</Parameters_filter>\n");
    }
    else
    printf("      <Parameters_filter>F</Parameters_filter>\n");
    printf("    </Parameters>\n");
    printf("  </BlastOutput_param>\n");
    printf("  <BlastOutput_iterations>\n");
    printf("    <Iteration>\n");
    printf("      <Iteration_iter-num>1</Iteration_iter-num>\n");
    printf("      <Iteration_hits>\n");
}

// Footer of XML output
void print_XMLfooter()
{
    printf("      </Iteration_hits>\n");
    printf("      <Iteration_stat>\n");
    printf("        <Statistics>\n");
    printf("          <Statistics_db-num>%d</Statistics_db-num>\n", readdb_numberOfSequences);
    printf("          <Statistics_db-len>%llu</Statistics_db-len>\n", statistics_databaseSize);
    printf("          <Statistics_hsp-len>%d</Statistics_hsp-len>\n", statistics_lengthAdjust);
    printf("          <Statistics_eff-space>%llu</Statistics_eff-space>\n", statistics_searchSpaceSize);
    printf("          <Statistics_kappa>%.3f</Statistics_kappa>\n", statistics_gappedParams.K);
    printf("          <Statistics_lambda>%.3f</Statistics_lambda>\n", statistics_gappedParams.lambda);
    printf("          <Statistics_entropy>%.3f</Statistics_entropy>\n", statistics_gappedParams.H);
    printf("        </Statistics>\n");
    printf("      </Iteration_stat>\n");
    printf("    </Iteration>\n");
    printf("  </BlastOutput_iterations>\n");
    printf("</BlastOutput>\n");
}

// Print a gapped extension using tabular output
void print_tabularGappedExtension(struct gappedExtension* gappedExtension, struct PSSMatrix PSSMatrix,
                                  char* query, unsigned char* subject, char* queryDescription,
                                  char* subjectDescription)
{
	char *queryLine, *subjectLine, *midLine;
	int4 identities = 0, positives = 0, gaps = 0, length, gapopens = 0;
	struct trace trace;
    char reverseComplement = 0;

	trace = gappedExtension->trace;

    if (trace.queryStart >= PSSMatrix.strandLength)
    {
    	trace.queryStart = PSSMatrix.length - trace.queryStart - 1;
    	reverseComplement = 1;
    }

    // Declare memory for query, subject and midlines
	queryLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	subjectLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	midLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	queryLine[0] = '\0'; subjectLine[0] = '\0'; midLine[0] = '\0';

    length = trace.length + 1;

    print_constructAlignment(queryLine, subjectLine, midLine, &identities, &positives, &gaps,
                             &gapopens, reverseComplement, trace, query, subject, PSSMatrix);

    // Fields: query id, subject ids, % identity, alignment length, mismatches, gap opens,
    //         q. start, q. end, s. start, s. end, evalue, bit score
    printf("%s	%s	%.2f	%d	%d	%d	%d	%d	%d	%d	%s	%.1f\n",
           queryDescription, subjectDescription, 100.0 * (float)identities / (float)length,
           length, length - identities, gapopens, trace.queryStart + 1,
           gappedExtension->queryEnd + 1, trace.subjectStart + 1, gappedExtension->subjectEnd + 1,
           print_eValue2String(gappedExtension->eValue), gappedExtension->normalizedScore);

    free(queryLine);
	free(subjectLine);
	free(midLine);
}

// Print a gapped extension using XML output
void print_XMLgappedExtension(struct gappedExtension* gappedExtension, struct PSSMatrix PSSMatrix,
                              char* query, unsigned char* subject, uint4 hspNum)
{
	char *queryLine, *subjectLine, *midLine;
	int4 identities = 0, positives = 0, gaps = 0, gapopens = 0;
	struct trace trace;
    char reverseComplement = 0;

	trace = gappedExtension->trace;

    if (trace.queryStart >= PSSMatrix.strandLength)
    {
    	trace.queryStart = PSSMatrix.length - trace.queryStart - 1;
    	reverseComplement = 1;
    }

    // Declare memory for query, subject and midlines
	queryLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	subjectLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	midLine = (char*)global_malloc(sizeof(char) * (trace.length + 2));
	queryLine[0] = '\0'; subjectLine[0] = '\0'; midLine[0] = '\0';

    print_constructAlignment(queryLine, subjectLine, midLine, &identities, &positives, &gaps,
                             &gapopens, reverseComplement, trace, query, subject, PSSMatrix);

    printf("            <Hsp>\n");
    printf("              <Hsp_num>%d</Hsp_num>\n", hspNum);
    printf("              <Hsp_bit-score>%.3f</Hsp_bit-score>\n", gappedExtension->normalizedScore);
    printf("              <Hsp_score>%d</Hsp_score>\n", gappedExtension->nominalScore);
    printf("              <Hsp_evalue>%g</Hsp_evalue>\n", gappedExtension->eValue);
    printf("              <Hsp_query-from>%d</Hsp_query-from>\n", trace.queryStart + 1);
    printf("              <Hsp_query-to>%d</Hsp_query-to>\n", gappedExtension->queryEnd + 1);
    printf("              <Hsp_hit-from>%d</Hsp_hit-from>\n", trace.subjectStart + 1);
    printf("              <Hsp_hit-to>%d</Hsp_hit-to>\n", gappedExtension->subjectEnd + 1);
    printf("              <Hsp_query-frame>1</Hsp_query-frame>\n");
    printf("              <Hsp_hit-frame>1</Hsp_hit-frame>\n");
    printf("              <Hsp_identity>%d</Hsp_identity>\n", identities);
    printf("              <Hsp_positive>%d</Hsp_positive>\n", positives);
    printf("              <Hsp_gaps>%d</Hsp_gaps>\n", gaps);
    printf("              <Hsp_align-len>%d</Hsp_align-len>\n", trace.length + 1);
    printf("              <Hsp_qseq>%s</Hsp_qseq>\n", queryLine);
    printf("              <Hsp_hseq>%s</Hsp_hseq>\n", subjectLine);
    printf("              <Hsp_midline>%s</Hsp_midline>\n", midLine);
    printf("            </Hsp>\n");

	free(queryLine);
	free(subjectLine);
	free(midLine);

}

void print_seqData(struct sequenceData *strSequence, uint4 nNum)
{
	int i;
	for (i = 0; i < nNum; i++)
	{
		printf("%d, seqLen = %d, desSt = %d, desLen = %d, codeLen = %d\n", 
				i,
				strSequence[i].sequenceLength,
				strSequence[i].descriptionStart,
				strSequence[i].descriptionLength,
				strSequence[i].encodedLength);

	}
	printf("\n");
}

void print_seqDataFP(struct sequenceDataFP *strSequence, uint4 nNum)
{
	int i;
	for (i = 0; i < nNum; i++)
	{
		printf("FP: %d, seqLen = %d, desSt = %d, desLen = %d, codeLen = %d, offset = %d\n", 
				i,
				strSequence[i].sequenceLength,
				strSequence[i].descriptionStart,
				strSequence[i].descriptionLength,
				strSequence[i].encodedLength,
				strSequence[i].offset);

	}
	printf("\n");
}

