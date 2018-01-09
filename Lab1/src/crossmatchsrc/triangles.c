#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "funcs.h"

/*
   Given two coordinate lists (x1,y1) and (x2,y2) (each nobj long)
   return index xx of objects from list2 identified in list1
*/

void triangles(x1,y1,x2,y2,nobj1,nobj2,xx,param)
     double *x1, *y1, *x2, *y2, *param;
     int nobj1, nobj2, *xx;
{
  double LLIM, RLIM, RTOL, LTOL, CTOL, FVNO;
  
  long int ntri1, ntri2, i, j, k, *index, idx;
  double *perim1, *perim2, *cos1, *cos2, *rat1, *rat2;
  double dperim, drat, dcos, maxntri;
  int *ori1, *ori2, *ver11, *ver12, *ver13, *ver21, *ver22, *ver23;
  int vote[MAX_COUNT], vidx[MAX_COUNT], tver=0, count, found, best, n=0;


  LLIM = param[0];
  RLIM = param[1];
  LTOL = param[2];
  RTOL = param[3];
  CTOL = param[4];
  FVNO = param[5];

  ntri1 = nobj1*(nobj1-1)*(nobj1-2)/6;
  ntri2 = nobj2*(nobj2-1)*(nobj2-2)/6;

  perim1 = (double *)malloc(ntri1*sizeof(double));
  perim2 = (double *)malloc(ntri2*sizeof(double));
  rat1   = (double *)malloc(ntri1*sizeof(double));
  rat2   = (double *)malloc(ntri2*sizeof(double));
  cos1   = (double *)malloc(ntri1*sizeof(double));
  cos2   = (double *)malloc(ntri2*sizeof(double));

  ori1   = (int *)malloc(ntri1*sizeof(int));
  ori2   = (int *)malloc(ntri2*sizeof(int));

  ver11  = (int *)malloc(ntri1*sizeof(int));
  ver12  = (int *)malloc(ntri1*sizeof(int));
  ver13  = (int *)malloc(ntri1*sizeof(int));

  ver21  = (int *)malloc(ntri2*sizeof(int));
  ver22  = (int *)malloc(ntri2*sizeof(int));
  ver23  = (int *)malloc(ntri2*sizeof(int));

  index  = (long int *)malloc(ntri1*sizeof(long int));

  /* generate triangles from both lists of points */

  tri_gen(x1,y1,nobj1,perim1,rat1,cos1,ori1,ver11,ver12,ver13,&ntri1);
  tri_gen(x2,y2,nobj2,perim2,rat2,cos2,ori2,ver21,ver22,ver23,&ntri2);


  for(i=0L;i<ntri1;i++) {
    index[i] = -1L;

    for(j=0L;j<ntri2;j++) {

      /*
         for triangles big and "round" enough look for close neighbors
         in (perimeter, longest/shortest side ratio, cosine angle
	 between longest/shortest side) space
      */

      if(rat1[i]<RLIM && rat2[j]<RLIM && perim1[i]>LLIM && perim2[j]>LLIM) {

        dperim = fabs(perim1[i] - perim2[j]);
        drat   = fabs(1.0 - rat1[i]/rat2[j]);
        dcos   = fabs(1.0 - cos1[i]/cos2[j]);

        if( dperim<LTOL && drat<RTOL && dcos<CTOL && ori1[i]==ori2[j] )
          index[i]=j;

      }
    }
  }

  for (i=0;i<ntri1;i++)
    if ( index[i] != -1L )
      n++;

  free(perim1); free(rat1); free(cos1); free(ori1);
  free(perim2); free(rat2); free(cos2); free(ori2);


  /*
     voting procedure: for each point check how many triangles were
     matched and had this point as one of its vertices
  */

  for(i=0L;i<nobj1;i++) {
    count = 0;

    for(j=0L;j<ntri1;j++) {
      if( (idx = index[j]) != -1L) {

        if(i == ver11[j])      tver = ver21[idx];
        else if(i == ver12[j]) tver = ver22[idx];
        else if(i == ver13[j]) tver = ver23[idx];
        else idx = -1L;

        if( idx != -1L &&  count < MAX_COUNT) {
          found = 0;
          for(k=0L;k<count;k++) {
            if(vidx[k] == tver) {
              vote[k]++;
              found = 1;
            }
          }
          
          if(!found) {
            vidx[count] = tver;
            vote[count] = 1;
            count++;
          }
        }

      }
    }

    if(count >= MAX_COUNT) printf("triangles: warning: MAX_COUNT reached \n");

    best   = 0;
    for(k=0L;k<count;k++)
      if( vote[k] > vote[best] )
        best = k;

    maxntri = 3.0 * (double)n / ( cbrt(6.0*(double)n) + 1.0 );

    if( count > 0  && vote[best] > FVNO*maxntri ) 
      xx[i] = vidx[best];
    else
      xx[i] = -1;

  }

  free(ver11); free(ver12); free(ver13);
  free(ver21); free(ver22); free(ver23);

  return ;

}


