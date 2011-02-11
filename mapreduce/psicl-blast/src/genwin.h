
/*****************************************************************************/
/***   (genwin.h)                                                          ***/
/*****************************************************************************/
/*--------------------------------------------------------------(includes)---*/

#include <stdio.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>

/*--------------------------------------------------------------(typedefs)---*/
//extern "C" {

typedef int Bool;

/*---------------------------------------------------------------(defines)---*/
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif
/*---------------------------------------------------------------(structs)---*/

struct Database
  {
   char *filename;
   FILE *fp;
   char *indexname;
   long int filepos;
  };

struct Sequence
  {
   struct Database *db;

   char *id;
   char *name;
   char *organism;
   char *header;

   char *seq;
   int start;                     /* for windows */
   int length;
   struct Alphabet *alphabet;

   struct Sequence *parent;       /* for windows */
   struct Sequence *root;         /* for subwindows */
   struct Sequence **children;    /* only the floaters? */

   Bool bogus;
   Bool punctuation;
   Bool rubberwin;             /* for windows */
   Bool floatwin;              /* for subwindows */
   Bool seedwin;

   int *state;
   double entropy;
   int *composition;
   Bool	*comptst;
   int	*compindex;

   char *classvec;               /* from ClaVec[aa] */
   struct Alphabet *clalphabet;
   double *scorevec;             /* from ScoVec[aa] or ScoFun(pos) */
   double score;                 /* from ScoFun(win) */
                          /* union, for integer scores? */
   struct PerScoreVec *foo;
   int *bar;
  };

struct Configuration
  {
   char *iseq;
   int ilength;

   int printper;
  };

struct Matrix
  {
   struct Sequence *parent;

   char **seq;
   int start;
   int period;
   int length;

   int **rowcomp;
   int **colcomp;
  };

/*------------------------------------------------------(alphabet structs)---*/

struct Alphabet
  {
   char *name;
   int size;

   int *index[128];
   char *chars;       /*  [size]  */
   char **charnames;  /*  [size]  */

   Bool caseinvariant;
  };

struct TransVector
  {
   struct Alphabet *from;
   struct Alphabet *to;

   int *index;     /*  [from->size]  */
  };

struct ScoreVector
  {
   struct Alphabet *from;

   double *score;   /*  [from->size]  */
   char *units;
  };

struct ScoreMatrix
  {
   struct Alphabet *from1;
   struct Alphabet *from2;

   double **score;  /*  [from1->size][from2->size]  */
   char *units;
  };

/*---------------------------------------------------------(bogus structs)---*/

struct PerScoreVec
  {
   int hits;
   int tot;
   double pct;
   double std;
   double prob;
  };

/*----------------------------------------------------------------(protos)---*/

extern struct Database *opendbase();
extern closedbase();

extern struct Sequence *openseq(), *firstseq(), *nextseq();
extern closeseq();

extern struct Sequence *openwin(), *nextwin();
extern shiftwin(), closewin();

extern compon(), stateon(), enton();
extern double entropy();

extern struct Matrix *openmat();
extern closemat();

extern upper(), lower(), findchar();

/*----------------------------------------------------------------(macros)---*/

/***   #define bogus(A) (aaindex[A]>=20)   ***/

/*---------------------------------------------------------------(globals)---*/

extern int aaindex[];
extern char aachar[];

//}
/*---------------------------------------------------------------------------*/
