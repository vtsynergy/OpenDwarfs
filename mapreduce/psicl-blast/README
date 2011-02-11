
For the latest information about FSA-BLAST visit http://www.fsa-blast.org/

 About FSA-BLAST

FSA-BLAST is a new version of the popular BLAST (Basic Local Alignment Search Tool) bioinformatics tool, used to search genomic databases containing either protein or nucleotide sequences. FSA stands for Faster Search Algorithm; FSA-BLAST is twice as fast as NCBI-BLAST for nucleotide searches and 30% faster for protein searches with no loss in accuracy. These speed gains are due to a range of improvements to the BLAST algorithm described in detail in recent scientific publications. The software is freely available for download and open source under the BSD license agreement.

The FSA-BLAST software is designed to be as similar as possible in usage to the NCBI-BLAST application. Results are output in an almost identical format. Most command line options are the same, and parameters such as word length, hit threshold, alignment dropoff and gapped alignment trigger are comparable to NCBI-BLAST. FSA-BLAST uses the exact same statistical analysis to evaluate E-values and normalized scores for alignments.

FSA-BLAST is currently under development. Some features are unavailable, including the ability to perform translated searches (BLASTX, TBLASTN, TBLASTX). The following features are currently available:

 Protein vs protein (BLASTP) and nucleotide vs. nucleotide searches (BLASTN)
 Query filtering using DUST and SEG low complexity filters
 Reverse complement nucleotide searches
 Karlin-Altschul alignment statistics
 Control over a range of parameters including; word length, hit threshold, multiple hit window size, ungapped and gapped dropoff, open and extend gap penalties, scoring matrix, match and mismatch penalties, number of reported/displayed alignments, and gapped alignment trigger score.

FSA-BLAST is written and maintained by Michael Cameron. Improvements to the BLAST algorithm are the result of research conducted by Michael Cameron, Hugh E. Williams, Yaniv Bernstein and Adam Cannane at RMIT University, Australia.

 Download and installing

FSA-BLAST is available for download as source code or precompiled binaries for x86 Linux systems:

fsablast101-source.tar.gz
    Source code, unix tar/gzip format
fsablast101-x86binaries.tar.gz
    Intel x86 Linux pre-compiled binaries, unix tar/gzip format
fsablast101-G5binaries.tar.gz
    PowerMac G5 pre-compiled binaries, unix tar/gzip format

All source code is freely available under the BSD license agreement. You can view a copy of the license agreement here.
The software has been tested on a handful of Unix systems. If you have any problems compiling the software under Unix please let us know.

Decompress the above files using the following commands:

gzip -d fsablast101-source.tar.gz 
tar -xf fsablast101-source.tar 

To compile the software type:

make 

Which will generate the following binaries:

formatdb
Converts a FASTA format database into files readable by FSA-BLAST
cluster
Clusters a protein collection for faster BLAST searches (around 22% faster for GenBank NR database)
blast
Performs BLASTP and BLASTN searches
readdb
Outputs the contents of a collection that has been processed using formatdb in FASTA format
ssearch
Performs Smith-Waterman search against a formatted collection using BLAST output, scoring and statistics

BLAST also needs to know the location of scoring matrix files, such as BLOSUM62. BLAST consults the file .ncbirc in the user's home directory to find the location of the scoring files. The .ncbirc file can be created using a text editor and should be formatted as follows:

[NCBI]
Data=/home/user/blast/data

where the directory specified contains the scoring matrix files. In absence of a .ncbirc file, FSA-BLAST will attempt to locate the files in the /data subdirectory of the current working directory.

 Usage instructions 

Before searching a collection, you will first need to format it using the formatdb tool provided with FSA-BLAST. Note that this tool is different from the formatdb application that comes with NCBI-BLAST which uses a different format. The tool will generate three files with extensions .data .sequences and .descriptions in the same directory as the collection.

The following illustrates how to format a collection:

$ ls -al ~/data
total 79256
drwxr-xr-x    2 mcam     mcam         4096 Sep  1 16:33 ./
drwxr-xr-x   24 mcam     mcam         4096 Sep  1 16:32 ../
-rw-r--r--    1 mcam     mcam     81061294 Sep  1 16:33 pdb

$ ./formatdb ~/data/pdb
PROTEIN database detected.
Formatting database..............................done.
277408 sequences processed.
50873918 letters processed.
0 wildcards encoded.
1 volume(s) created.
Longest/shortest sequence was 15281/6 letters

$ ls -al ~/data
total 158828
drwxr-xr-x    2 mcam     mcam         4096 Sep  1 16:36 ./
drwxr-xr-x   24 mcam     mcam         4096 Sep  1 16:32 ../
-rw-r--r--    1 mcam     mcam     81061294 Sep  1 16:33 pdb
-rw-r--r--    1 mcam     mcam      1121748 Sep  1 16:36 pdb.data
-rw-r--r--    1 mcam     mcam     28832910 Sep  1 16:36 pdb.descriptions
-rw-r--r--    1 mcam     mcam     51428734 Sep  1 16:36 pdb.sequences

NOTE: To format a collection for use with FSA-BLAST you will need the database to be in FASTA format. To convert a database from NCBI-BLAST format (with files extensions such as nhr, nin, nsq, phr, pin, and psq) back to FASTA format you will need to use the fastacmd tool that comes with the NCBI toolkit (in the /build directory).

Once the collection has been formatted you can search it with blast using a command such as follows:

$ ./blast -i query -d ~/data/pdb

Which will produce output looking like:

CW-BLAST

Query= gi|38083732|ref|XP_357594.1| similar to KIAA0960 protein [Mus
musculus]
         (60 letters)

Database: /home/mcam/data/pdb
           277,408 sequences; 50,873,918 total letters

Searching....................................................done.



                                                                 Score    E
Sequences producing significant alignments:                      (bits) Value

gi|62777659|gb|AAY03565.1| Sequence 105 from patent US 6861256 ...    34  0.062
gi|40115719|gb|AAR55617.1| Sequence 53 from patent US 6613544 >...    28  5.8

>gi|62777659|gb|AAY03565.1| Sequence 105 from patent US 6861256
           >gi|15113415|gb|AAE68999.1| Sequence 105 from patent US
           6225120
          Length = 383 DescriptionLocation = 17309878

 Score = 34.3 bits (77), Expect = 0.062
 Identities = 17/55 (30%), Positives = 26/55 (47%), Gaps = 5/55 (9%)
 Strand = Plus / Plus

Query: 9   DRQLCRDAIFCDAPCPKDCVL-----WSSCSHTCSGKSILAYAGEEGGIRCPNIS 58
           +R+LC ++  C   CP+ C        + CS  C G  ++   G E  I C N+S
Sbjct: 76  NRRLCWNSKLCQTKCPEKCRNNCIDEHTCCSQDCLGGCVIDKNGNESCISCRNVS 130


>gi|40115719|gb|AAR55617.1| Sequence 53 from patent US 6613544
           >gi|21504139|gb|AAM56900.1| Sequence 139 from patent US
           6369027 >gi|17906295|gb|AAE81374.1| Sequence 52 from
           patent US 6288032
          Length = 191 DescriptionLocation = 23345524

 Score = 27.7 bits (60), Expect = 5.8
 Identities = 17/56 (30%), Positives = 26/56 (46%), Gaps = 5/56 (8%)
 Strand = Plus / Plus

Query: 2   CFPGEEVDRQLCRD---AIFCDAPCPKDCVLWSSCSHTCSGKSILAYAGEEGGIRC 54
           C PG+E+ +Q C+      F D      C  W++CS    G+S+L     E  + C
Sbjct: 105 CRPGQELTKQGCKTCSLGTFNDQNGTGVCRPWTNCS--LDGRSVLKTGTTEKDVVC 158


  Database: /home/mcam/data/pdb  (1 volumes)
  Number of letters in database: 50,873,918
  Number of sequences in database:  277408

Lambda     K      H     (ungapped)
 0.325     0.140  0.498

Lambda     K      H     (gapped)
 0.267     0.041  0.140


Matrix: BLOSUM62
Gap Penalties: Existence: 11, Extension: 1
Semi-Gapped Gap Penalties: Existence: 7, Extension: 1
Number of Hits to DB: 7,326,899
Number of Sequences: 277408
Number of extensions: 295218
Number of successful extensions: 969
Number of sequences with successful extensions: 830
Number of sequences with semi-gapped score above cutoff: 830
Number of sequences better than 10: 2
Number of HSP's that attempted semi-gapping: 963
Number of HSP's that attempted gapping: 961
Number of HSP's contained and not gapped: 6
Number of HSP's succeeded/attempted join: 0/0
Total subject bytes copied/unpacked = 0/0
length of query: 60
length of database: 50,873,918
effective HSP length: 30
effective length of query: 30
effective length of database: 42,551,678
effective search space: 1276550340
effective search space used: 1276550340
T: 11
A: 40
X1: 15
X2: 38
X3: 64
S1: 40
S2: 58
F2: 40

To view a complete list of BLAST parameters simple execute:

$ ./blast

CW-BLAST

  -d  Database
    default =
  -i  Query File
    default =
  -A  Multiple Hits window size (protein only)
    default = 40
  -f  Threshold for extending hits (protein only)
    default = 11
  -e  Expectation value (E)
    default = 10.0
  -y  Dropoff (X) for blast extensions in bits
    default = 7.0 for protein, 20.0 for nucleotide
  -P  0 for multiple hit, 1 for single hit
    default = 0
  -G  Cost to open a gap (Matrix dependant)
    default = 11 for protein, 5 for nucleotide
  -E  Cost to extend a gap (Matrix dependant)
    default = 1 for protein, 2 for nucleotide
  -O  Cost to open a gap (Semi-gapped alignment, Matrix dependant)
    default = 7
  -X  X dropoff value for gapped alignment (in bits)
    default = 15.0 for protein, 30.0 for nucleotide
  -N  Number of bits to trigger gapping
    default = 22.0 for protein, 25.0 for nucleotide
  -Z  X dropoff value for final gapped alignment (in bits)
    default = 25.0 for protein, 50.0 for nucleotide
  -M  Matrix
    default = BLOSUM62
  -v  Number of database sequences to show one-line descriptions for (V)
    default = 500
  -b  Number of database sequences to show alignments for (B)
    default = 250
  -W  Word size
    default = 3 for protein, 11 for nucleotide
  -z  Effective length of the database (use zero for the real size)
    default = 0
  -q  Penalty for a nucleotide mismatch
    default = -3
  -r  Reward for a nucleotide match
    default = 1
  -S  Query strands to search against database (nucleotide only)
      3 is both, 1 is top, 2 is bottom
    default = 3

ERROR: Query File not specified

To perform faster protein BLAST searches, you can cluster the collection using the cluster command:

 ./cluster ~/data/pdb

Number of sequences = 1822
Total number of letters = 574741
Length of longest sequence = 4560
Alphabet type = Protein
1822 sequences read.
Iteration 0 WordLength=12
Hashcounter state:
Singles: 13051
Duplicates: 18264
Iteration time=0.03 secs
Iteration 1 WordLength=21
Hashcounter state:
Singles: 8892
Duplicates: 13927
Iteration time=0.25 secs
Iteration 2 WordLength=30
Hashcounter state:
Singles: 6589
Duplicates: 10720
Iteration time=0.43 secs
Initialized diagonals
Constructing and processing postings lists..................................................done.
Identifying high-scoring sequence pairs...done. (5500 pairs)
Performing hierarchical clustering...done.
Total bytes saved=152948
Writing clusters to disk...done.
Writing remaining sequences to disk...done.

BLAST searches against the clustered database will then be faster. The amount of redundancy in the original collection will affect the speed increase obtained by clustering, although our experiments have shown a 22% speed increase when searching the GenBank NR database.

Also provided is a tool for converting a formatted collection back into FASTA format. The command:

 ./readdb ~/data/pdb

will output the database to stdout in FASTA format.

 Papers

The following papers describe improvements to the BLAST algorithm used by FSA-BLAST to increase search speed without any loss in accuracy:

M. Cameron, H.E. Williams, and A. Cannane, ``Improved Gapped Alignment in BLAST'', IEEE/ACM Transactions on Computational Biology and Bioinformatics, 1(3), 116-129, 2004. Postscript PDF Abstract and source code

M. Cameron, H.E. Williams, and A. Cannane, ``A Deterministic Finite Automaton for Faster Protein Hit Detection in BLAST'', manuscript in preparation.

M. Cameron and H.E. Williams, ``Comparing compressed sequences for faster nucleotide BLAST searches'', manuscript in preparation.

The original papers describing the BLAST algorithm are:

Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. & Lipman, D.J. (1997) "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs." Nucleic Acids Res. 25:3389-3402. Medline

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990) "Basic local alignment search tool." J. Mol. Biol. 215:403-410.  Medline

 Feedback & reporting bugs

If you have any feedback regarding FSA-BLAST including bug reports, questions or feature requests email them to mcam at cs.rmit.edu.au

When reporting a problem compiling or running FSA-BLAST please include as much information about the problem as possible, including:

 A detailed description of what caused the problem.
 The architecture and operating system of the system used.
 What collections and query sequences caused the problem.
 A copy of the output when the software crashed.












