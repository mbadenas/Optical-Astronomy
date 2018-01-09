#include <stdio.h>
#include <stdlib.h>
#include "funcs.h"

/*
   Given a list of nd objects (position, magnitude)=(x,y,mag),
   return a similar list for nobj brightest objects.
*/

void bright_end(x,y,mag,nd,xs,ys,nobj,flag)
     double *x, *y, *mag, *xs, *ys;
     int nd, *nobj, *flag;
{
  void indexx();

  int i, *idx, gi;


  if((*nobj)>nd) printf("bright_end: error: not enough objects in the list \n");
  idx = (int *)malloc(nd*sizeof(int));

  /* use little trick to fool stupid 1,n indexing of numerical recipes */

  indexx(nd,mag-1,idx-1);

  gi = 0;
  for(i=0;i<(*nobj);i++)
  { 
    if (flag[idx[i]-1])
    {
      xs[gi] = x[idx[i]-1];
      ys[gi] = y[idx[i]-1];
      gi++;
    }
    else if (nd > (*nobj))
      (*nobj)++;
  }
  
  *nobj = gi;
  
  free(idx);

  return ;
}

