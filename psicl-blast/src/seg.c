
/*****************************************************************************/
/***  (seg.c)                                                              ***/
/*** #include precomputed ln(fact(n)) array in "lnfac.h"                   ***/
/*****************************************************************************/

/*--------------------------------------------------------------(includes)---*/

#include "genwin.h"
#include "genwin.c"
#include "lnfac.h"

/*---------------------------------------------------------------(defines)---*/

#define LENCORRLIM 120
#define MIN(a,b)	((a) <= (b) ? (a) : (b))

/*---------------------------------------------------------------(structs)---*/

struct Segment
  {
   int begin;
   int end;
   struct Segment *next;
  };

/*---------------------------------------------------------------(globals)---*/

int window = 12;
int downset, upset;
double locut = 2.2;
double hicut = 2.5;

int hilenmin = 0;
int overlaps = FALSE;
int hionly = FALSE;
int loonly = FALSE;
int entinfo = TRUE;
int singleseq = FALSE;
int prettyseq = FALSE;
int prettytree = TRUE;
int charline = 60;
int maxtrim = 100;

double getprob(), lnperm(), lnass();

void seg_segSequence(char* sequence)
{
    struct Sequence *seq;
    struct Segment *segs;

   seq = (struct Sequence *) malloc(sizeof(struct Sequence));
    seq->seq = sequence;
    seq->length = strlen(sequence);
   seq->db = NULL;
   seq->parent = (struct Sequence *) NULL;
   seq->root = (struct Sequence *) NULL;
   seq->children = (struct Sequence **) NULL;
   seq->rubberwin = FALSE;
   seq->floatwin = FALSE;

   seq->punctuation = FALSE;
   seq->entropy = -2.;
   seq->state = (int *) NULL;
   seq->composition = (int *) NULL;
   seq->classvec = (char *) NULL;
   seq->scorevec = (double *) NULL;

   genwininit();

   downset = (window+1)/2 - 1;
   upset = window - downset;

    singleseq = TRUE;
    prettyseq = FALSE;
    prettytree = FALSE;
    hionly = TRUE;
    loonly = FALSE;

   entropy_init(window);

    segs = (struct Segment *) NULL;
    segseq(seq, &segs, 0);
    mergesegs(seq, segs);

    singreport(seq, segs);

    freesegs(segs);
}

getparams(argc, argv)
  int argc;
  char *argv[];

  {int i;
   int nargc;
   char **nargv;
   extern char *optarg;
   extern int optind;
   int c;

   if (argc<2)
     {
      usage();
      exit(1);
     }

   for (i=2; argc>i && argv[i][0]!='-'; i++)
     {
      if (i==2)
        {
         window = atoi(argv[2]);
        }
      else if (i==3)
        {
         locut = atof(argv[3]);
        }
      else if (i==4)
        {
         hicut = atof(argv[4]);
        }
     }

   if (locut>hicut)
     {
      fprintf(stderr, "Warning: trigger complexity greater than extension\n");
      hicut = locut;
     }

   downset = (window+1)/2 - 1;
   upset = window - downset;

    singleseq = TRUE;
    prettyseq = FALSE;
    prettytree = FALSE;
    hionly = TRUE;
    loonly = FALSE;
}
/*---------------------------------------------------------------(segment)---*/

segseq(seq, segs, offset)
  struct Sequence *seq;
  struct Segment **segs;
  int offset;

  {struct Segment *seg, *leftsegs;
   struct Sequence *leftseq;
   int first, last, lowlim;
   int loi, hii, i;
   int leftend, rightend, lend, rend;
   double *H, *seqent();

   H = seqent(seq);
   if (H==NULL) return;

   first = downset;
   last = seq->length - upset;
   lowlim = first;

   for (i=first; i<=last; i++)
     {
      if (H[i]<=locut && H[i]!=-1)
        {
         loi = findlo(i, lowlim, H);
         hii = findhi(i, last, H);

         leftend = loi - downset;
         rightend = hii + upset - 1;

         trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);

         if (i+upset-1<leftend)   /* check for trigger window in left trim */
           {
            lend = loi - downset;
            rend = leftend - 1;

            leftseq = openwin(seq, lend, rend-lend+1);
            leftsegs = (struct Segment *) NULL;
            segseq(leftseq, &leftsegs, offset+lend);
            if (leftsegs!=NULL)
              {
               if (*segs==NULL) *segs = leftsegs;
               else appendseg(*segs, leftsegs);
              }
            closewin(leftseq);

/*          trim(openwin(seq, lend, rend-lend+1), &lend, &rend);
            seg = (struct Segment *) malloc(sizeof(struct Segment));
            seg->begin = lend;
            seg->end = rend;
            seg->next = (struct Segment *) NULL;
            if (segs==NULL) segs = seg;
            else appendseg(segs, seg);  */
           }

         seg = (struct Segment *) malloc(sizeof(struct Segment));
         seg->begin = leftend + offset;
         seg->end = rightend + offset;
         seg->next = (struct Segment *) NULL;

         if (*segs==NULL) *segs = seg;
         else appendseg(*segs, seg);

         i = min(hii, rightend+downset);
         lowlim = i + 1;
/*       i = hii;     this ignores the trimmed residues... */
        }
     }

   free(H);
   return;
  }

/*----------------------------------------------------------------(seqent)---*/

double *seqent(seq)
  struct Sequence *seq;

  {struct Sequence *win;
   double *H;
   int i, first, last;

   if (window>seq->length)
     {
      return((double *) NULL);
     }

   H = (double *) malloc(seq->length*sizeof(double));

   for (i=0; i<seq->length; i++)
     {
      H[i] = -1.;
     }

   win = openwin(seq, 0, window);
   enton(win);

   first = downset;
   last = seq->length - upset;

   for (i=first; i<=last; i++)
     {
      if (seq->punctuation && hasdash(win))
        {H[i] = -1;
         shiftwin1(win);
         continue;}
      H[i] = win->entropy;
      shiftwin1(win);
     }

   closewin(win);
   return(H);
  }

/*---------------------------------------------------------------(hasdash)---*/

hasdash(win)
  struct Sequence *win;
{
	register char	*seq, *seqmax;

	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax) {
		if (*seq++ == '-')
			return TRUE;
	}
	return FALSE;
}

/*----------------------------------------------------------------(findlo)---*/

findlo(i, limit, H)
  int i, limit;
  double *H;

  {int j;

   for (j=i; j>=limit; j--)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j+1);
  }

/*----------------------------------------------------------------(findhi)---*/

findhi(i, limit, H)
  int i, limit;
  double *H;

  {int j;

   for (j=i; j<=limit; j++)
     {
      if (H[j]==-1) break;
      if (H[j]>hicut) break;
     }

   return(j-1);
  }

/*------------------------------------------------------------------(trim)---*/

trim(seq, leftend, rightend)
  struct Sequence *seq;
  int *leftend, *rightend;

  {struct Sequence *win;
   double prob, minprob;
   int shift, len, i;
   int lend, rend;
   int minlen;

/* fprintf(stderr, "%d %d\n", *leftend, *rightend);  */

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;

   minprob = 1.;
   for (len=seq->length; len>minlen; len--)
     {
      win = openwin(seq, 0, len);
      stateon(win);
      i = 0;

      shift = TRUE;
      while (shift)
        {
         prob = getprob(win->state, len);
         if (prob<minprob)
           {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
           }
         shift = shiftwin1(win);
         i++;
        }
      closewin(win);
     }

/* fprintf(stderr, "%d-%d ", *leftend, *rightend);  */

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

/* fprintf(stderr, "%d-%d\n", *leftend, *rightend);  */

   closewin(seq);
   return;
  }

/*---------------------------------------------------------------(getprob)---*/

double getprob(sv, total)
  int *sv;
  int total;

  {double ans, totseq;

#define LN20	2.9957322735539909
   totseq = ((double) total) * LN20;

   ans = lnass(sv) + lnperm(sv, total) - totseq;

   return(ans);
  }

/*----------------------------------------------------------------(lnperm)---*/

double lnperm(sv, tot)
  int *sv;
   int tot;

  {double ans;
   int i;

   ans = lnfac[tot];

   for (i=0; sv[i]!=0; i++) 
     {
      ans -= lnfac[sv[i]];
     }

   return(ans);
  }

/*-----------------------------------------------------------------(lnass)---*/

double lnass(sv)
	register int	*sv;
{
	double	ans;
	register int	svi, svim1;
	register int	class, total;
	register int    i;

	ans = lnfac[20];
	if (sv[0] == 0)
		return ans;

	total = 20;
	class = 1;
	svim1 = sv[0];
	for (i=0;; svim1 = svi) {
	        if (++i==20) {
		        ans -= lnfac[class];
                        break;
		      }
		else if ((svi = *++sv) == svim1) {
			class++;
			continue;
		}
		else {
			total -= class;
			ans -= lnfac[class];
			if (svi == 0) {
				ans -= lnfac[total];
				break;
			}
			else {
				class = 1;
				continue;
			}
		}
	}

	return ans;
}

/*-------------------------------------------------------------(mergesegs)---*/

mergesegs(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {struct Segment *seg, *nextseg;
   int len;

   if (overlaps) return;
   if (segs==NULL) return;

   if (segs->begin<hilenmin) segs->begin = 0;

   seg = segs;
   nextseg = seg->next;

   while (nextseg!=NULL)
     {
      if (seg->end>=nextseg->begin)               /* overlapping segments */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      len = nextseg->begin - seg->end - 1;
      if (len<hilenmin)                            /* short hient segment */
        {
         seg->end = nextseg->end;
         seg->next = nextseg->next;
         free(nextseg);
         nextseg = seg->next;
         continue;
        }
      seg = nextseg;
      nextseg = seg->next;
     }

   len = seq->length - seg->end - 1;
   if (len<hilenmin) seg->end = seq->length - 1;

   return;
  }

/*----------------------------------------------------------------(report)---*/

report(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

  {struct Sequence *subseq;
   struct Segment *seg, *nextseg;
   static int hi = 1;
   static int lo = 0;

   if (segs==NULL)
     {
      enton(seq);
      seqout(seq, hi, 1, seq->length);
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   if (segs->begin>0)
     {
      subseq = openwin(seq, 0, segs->begin);
      enton(subseq);
      seqout(subseq, hi, 1, segs->begin);
      closewin(subseq);
     }

   for (seg=segs; seg!=NULL; seg=seg->next)
     {
      subseq = openwin(seq, seg->begin, seg->end-seg->begin+1);
      enton(subseq);
      seqout(subseq, lo, seg->begin+1, seg->end+1);
      closewin(subseq);

      if (seg->next==NULL)
        {
         break;
        }

      nextseg = seg->next;
      
      if (nextseg->begin<=seg->end)
        {
         fprintf(stderr, "Overlapping segments: %s\n", seq->id);
         continue;
        }

      if (nextseg->begin==seg->end+1)
        {
         continue;
        }

      subseq = openwin(seq, seg->end+1, nextseg->begin-seg->end-1);
      enton(subseq);
      seqout(subseq, hi, seg->end+2, nextseg->begin);
      closewin(subseq);
     }

   if (seg->end+1==seq->length)
     {
/*    fputc('\n', stdout);   -- for spacing after each sequence */
      return;
     }

   subseq = openwin(seq, seg->end+1, seq->length-seg->end-1);
   enton(subseq);
   seqout(subseq, hi, seg->end+2, seq->length);
   closewin(subseq);

/* fputc('\n', stdout);   -- for spacing after each sequence */
   return;
  }

/*------------------------------------------------------------(singreport)---*/

singreport(seq, segs)
	struct Sequence	*seq;
	struct Segment	*segs;
{
	char	*proseq, *proseqmax;
	struct Segment	*seg;
	int	begin, end, i, ctr;

	proseq = seq->seq;
	proseqmax = proseq + seq->length;
	upper(proseq, seq->length);

	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		memset(proseq + begin, 'x', end - begin +1);
	}

/*	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==charline) {
			putc('\n', stdout);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}

	putc('\n', stdout);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}*/
}

/*----------------------------------------------------------(prettyreport)---*/

prettyreport(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

{
	char	*proseq, *proseqmax;
	char	format[10];
	struct Segment	*seg;
	int	begin, end, i, ctr;
	int	leftspace;

	leftspace = (int) ceil(log10((double)seq->length));
	sprintf(format, "%%%dd ", leftspace);

	upper(proseq = seq->seq, seq->length);

	for (seg=segs; seg!=NULL; seg=seg->next) {
		begin = seg->begin;
		end = seg->end;
		lower(proseq + begin, end - begin + 1);
	}

	fprintf(stdout, "%s\n", seq->header);

	space(leftspace+1);
	for (i=0, ctr=1; i<charline; ++i, ++ctr) {
		if (ctr==10) {
			putc('.', stdout);
			ctr = 0;
		}
		else
			putc(' ', stdout);
	}
	putc('\n', stdout);
	fprintf(stdout, format, 1);

	proseqmax = proseq + seq->length;
	for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq) {
		if (ctr==charline) {
			putc('\n', stdout);
			fprintf(stdout, format, i+1);
			ctr = 0;
		}
		putc(*proseq, stdout);
	}

	fprintf(stdout, " %d\n", seq->length);
	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*---------------------------------------------------------(pretreereport)---*/

pretreereport(seq, segs)
  struct Sequence *seq;
  struct Segment *segs;

{
	struct Sequence	*win;
	struct Segment	*seg;
	char	buffer[100], leftfmt[20], rightfmt[20];
	char	*curseq;
	int	i, left, right, len;
	int	current, nextloent;
	int	cline;

	cline = charline / 2;
   
	fprintf(stdout, "%s\n\n", seq->header);
	sprintf(leftfmt, "%%%ds", cline);
	sprintf(rightfmt, "%%-%ds", cline);

	current = 0;

	for (seg=segs; ; seg=seg->next) {
		if (seg==NULL)
			nextloent = seq->length;
		else
			nextloent = seg->begin;

		if (current < nextloent) {
			left = current;
			right = nextloent - 1;
			len = right - left + 1;
			win = openwin(seq, left, len);
			upper(curseq = win->seq, win->length);

			space(cline);
			fprintf(stdout, " %4d-%-4d ", left+1, right+1);

			while (len>0) {
				if (len<=cline) {
					fwrite(curseq, len, 1, stdout);
					putc('\n', stdout);
					break;
				}
				else {
					fwrite(curseq, cline, 1, stdout);
					putc('\n', stdout);
					space(cline+11);
					curseq += cline;
					len -= cline;
				}
			}
			closewin(win);
		}

		if (seg==NULL) break;

      left = seg->begin;
      right = seg->end;
      len = right - left + 1;
      win = openwin(seq, left, len);
      lower(curseq = win->seq, win->length);

		i = MIN(cline, len);
		if (i < cline)
			fprintf(stdout, "%*s", cline-i, "");
		fwrite(curseq, i, 1, stdout);
		fprintf(stdout, " %4d-%-4d ", left+1, right+1);
		putc('\n', stdout);

		len -= cline;
		if (len>0)
			curseq += cline;

		while (len>0) {
			i = MIN(cline, len);
			if (i < cline)
				space(cline - i);
			fwrite(curseq, i, 1, stdout);
			putc('\n', stdout);
			len -= i;
			if (len>0)
				curseq += i;
		}

		closewin(win);
		current = right + 1;
	}

	putc('\n', stdout);
}

/*-----------------------------------------------------------------(space)---*/

space(len)
	register int len;
{
	static char	spaces[] =
		{' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',
		' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '
		};
	register int i;

	while (len > 0) {
		i = MIN(len, sizeof(spaces)/sizeof(spaces[0]));
		fwrite(spaces, i, 1, stdout);
		len -= i;
	}
}

/*----------------------------------------------------------------(seqout)---*/

seqout(seq, hilo, begin, end)
  struct Sequence *seq;
  int hilo;
  int begin, end;

#define HDRLEN 60

{
	char	*proseq, *proseqmax, *id, *header;
   char outbuf[HDRLEN+1];
   static int hi = 1;
   static int lo = 0;
   int i, ctr, iend;

   if (hionly && hilo==lo) return;
   if (loonly && hilo==hi) return;

   proseq = seq->seq;
   proseqmax = proseq + seq->length;
   id = seq->id;
   if (id==NULL) id = seq->parent->id;
   header = seq->header;
   if (header==NULL) header = seq->parent->header;

   iend = findchar(header, ' ');
   if (iend!=-1) header = header+iend;

   if (entinfo)
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
/*    if (iend!=-1 && strlen(header)<=HDRLEN) fprintf(stdout, "%s", header);
      else if (iend!=-1) for (i=0; i<HDRLEN; i++) putc(header[i], stdout); */
      fprintf(stdout, " complexity=%4.2f (%d/%4.2f/%4.2f)\n",
           seq->entropy, window, locut, hicut);
     }
   else
     {
      fprintf(stdout, ">%s(%d-%d)", id, begin, end);
      if (iend!=-1)   /* fprintf(stdout, "%s\n", header); */
        {
		 i = MIN(HDRLEN, strlen(header));
		 fwrite(header, i, 1, stdout);
         putc('\n', stdout);
        }
      else putc('\n', stdout);
     }
   
   if (hilo==lo)
     {
      lower(proseq, seq->length);
     }
   else if (hilo==hi && seq->length>=hilenmin)
     {
      upper(proseq, seq->length);
     }
   else
     {
      lower(proseq, seq->length);
     }

   for (; proseq < proseqmax; proseq+=i) {
		i = MIN(charline, proseqmax - proseq);
		fwrite(proseq, i, 1, stdout);
		putc('\n', stdout);
	}

	if (putc('\n', stdout) == EOF) {
		fprintf(stderr, "premature EOF on write\n");
		exit(2);
	}
}

/*-------------------------------------------------------------(appendseg)---*/

appendseg(segs, seg)
  struct Segment *segs, *seg;

  {struct Segment *temp;

   temp = segs;
   while (1)
     {
      if (temp->next==NULL)
        {
         temp->next = seg;
         break;
        }
      else
        {
         temp = temp->next;
        }
     }

   return;
  }

/*--------------------------------------------------------------(freesegs)---*/

freesegs(segs)
  struct Segment *segs;

  {struct Segment *temp;

   while (segs!=NULL)
     {
      temp = segs->next;
      free(segs);
      segs = temp;
     }
  }

/*-----------------------------------------------------------------(usage)---*/

usage()

  {
   fprintf(stderr, "\
Usage: seg <file> <window> <locut> <hicut> <options>\n\
         <file>   - filename containing fasta-formatted sequence(s) \n\
         <window> - OPTIONAL window size (default 12) \n\
         <locut>  - OPTIONAL low (trigger) complexity (default 2.2) \n\
         <hicut>  - OPTIONAL high (extension) complexity (default 2.5) \n\
	 <options> \n\
            -x  each input sequence is represented by a single output \n\
                sequence with low-complexity regions replaced by \n\
                strings of 'x' characters \n\
            -c <chars> number of sequence characters/line (default 60)\n\
            -m <size> minimum length for a high-complexity segment \n\
                (default 0).  Shorter segments are merged with adjacent \n\
                low-complexity segments \n\
            -l  show only low-complexity segments (fasta format) \n\
            -h  show only high-complexity segments (fasta format) \n\
            -a  show all segments (fasta format) \n\
            -n  do not add complexity information to the header line \n\
            -o  show overlapping low-complexity segments (default merge) \n\
            -t <maxtrim> maximum trimming of raw segment (default 100) \n\
            -p  prettyprint each segmented sequence (tree format) \n\
            -q  prettyprint each segmented sequence (block format) \n");
   exit(1);
  }

