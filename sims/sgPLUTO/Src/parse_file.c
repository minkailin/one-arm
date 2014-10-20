#include "pluto.h"

static int nlines;
static char **fline;

/* ********************************************************* */
int ParOpen (char *fname)
/*
 * 
 *  PURPOSE:
 *
 *   Parse file *fname and store its content
 *   line by line in *fline
 *
 *********************************************************** */
{
  char  sline[512];
  FILE *fp;

  if (fline == NULL) fline = ARRAY_2D(128,128,char);

  fp = fopen(fname,"r");
  if (fp == NULL) {
    printf (" ! Error: file %s not found\n", fname);
    QUIT_PLUTO(1);
  }
  nlines = 0;

  while ( fgets(sline, 512, fp) != NULL ) {
    if (strlen(sline) > 0) {
      strcpy (fline[nlines],sline); 
      nlines++;
    }
  }
  fclose(fp);

  return(nlines);
}

/* ****************************************************** */
char *ParGet (const char *label, int pos)
/* 
 *
 *  PURPOSE
 *
 *   Search for *label in all the lines pointed to 
 *   by **fline. 
 *   If label exists, return a pointer to the
 *   string located pos words after label.
 *   Return null if the word does not exist.
 *
 ******************************************************** */
{
  int  k, success;
  char   sline[512], *str;
  static char par[128];
  const char  delimiters[] = " \t\r\f";

/* ---------------------------------------
     search if label exists
   --------------------------------------- */

  success = 0;
  for (k = 0; k < nlines; k++) {
    sprintf (sline,"%s",fline[k]);
    str = strtok(sline,delimiters);

    if (strcmp(str,label) == 0) {
      success = 1; 
      break;
    }
  }

  if (!success) {
    printf ("! Label '%s' was not found\n",label);
    QUIT_PLUTO(1);
  }

  for (k = 0; k < pos; k++) str = strtok (NULL, delimiters);
    
  if (str == NULL || (strlen(str) == 1 && str[0] == '\n')) {
    printf ("! A field (# %d) is missing after '%s'\n",pos,label);
    QUIT_PLUTO(1);
 /*   return (NULL); */
  }

  sprintf (par,"%s",strtok(str,"\n")); /* -- get rid of newline character -- */

  return (par);    
}

/* ****************************************************** */
int ParQuery (const char *label)
/* 
 *
 *  PURPOSE
 *
 *   Search for *label in all the lines pointed to 
 *   by **fline. Return (1) if it exists.
 *
 ******************************************************** */
{
  int    k;
  char   sline[512], *str;
  const char  delimiters[] = " \t\r\f";

/* ---------------------------------------
     search if label exists
   --------------------------------------- */

  for (k = 0; k < nlines; k++) {
    sprintf (sline,"%s",fline[k]);
    str = strtok(sline,delimiters);

    if (strcmp(str,label) == 0) return(1);
  }

  return (0);
}

