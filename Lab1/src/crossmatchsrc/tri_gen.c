#include <math.h>
#include "funcs.h"

#define PI 3.1415927

/*
   Triangle generation. From the list of nobj points (x,y) calculate
   perimeters, longest/shortest side ratios, cosines of the angle
   between longest/shortest sides and orientations for all possible
   triangles. Returns also indexes of vertices ver1-2-3 opposite
   to the shortest, mid and longest side and number of triangles.
*/

void tri_gen(x,y,nobj,perim,rat,cos,ori,ver1,ver2,ver3,ntri)
     double *x, *y, *perim, *rat, *cos;
     int nobj, *ori, *ver1, *ver2, *ver3;
     long int *ntri;
{
  int i, j, k, idx[3], vi[3];
  long int n;
  double xc, yc, dx, dy, r[3], s[3];


  n=0;

  for(i=0;i<nobj;i++) {
    for(j=i+1;j<nobj;j++) {
      for(k=j+1;k<nobj;k++) {

	vi[0] = k;
	vi[1] = j;
	vi[2] = i;

	dx    = x[j] - x[i];
	dy    = y[j] - y[i];
	r[0]  = sqrt(dx*dx + dy*dy);

	dx    = x[k] - x[i];
	dy    = y[k] - y[i];
	r[1]  = sqrt(dx*dx + dy*dy);

	dx    = x[k] - x[j];
	dy    = y[k] - y[j];
	r[2]  = sqrt(dx*dx + dy*dy);

	xc = (x[i] + x[j] + x[k])/3;
	yc = (y[i] + y[j] + y[k])/3;

	/* use little trick to fool stupid 1,n indexing of numerical recipes */

	indexx(3,r-1,idx-1);

	ver1[n] = vi[idx[0]-1];
	ver2[n] = vi[idx[1]-1]; 
	ver3[n] = vi[idx[2]-1]; 

	s[0] = r[idx[0]-1]; 
	s[1] = r[idx[1]-1]; 
	s[2] = r[idx[2]-1]; 

	ori[n]  = (x[ver2[n]] - xc)*(y[ver1[n]] - yc);
	ori[n] -= (x[ver1[n]] - xc)*(y[ver2[n]] - yc);
	ori[n] /= fabs(ori[n]);

	perim[n] = s[2];
	rat[n]   = s[2]/s[0];
	cos[n]   = (s[0]*s[0] + s[2]*s[2] - s[1]*s[1])/(2.0*s[0]*s[2]);

	n++;
      }
    }
  }

  *ntri = n;

  return;

}

