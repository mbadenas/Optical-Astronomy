#include <math.h>

/*
   Given first order transformation in coeffx and coeffy
   transform list2 of nobj2 objects (x2,y2) into approximate frame
   of reference of list1 of nobj1 objects (x1,y1) and look
   for close neighbors to make a refined matching of lists.
   TOL is size of matching box and upon return index contains
   subscripts of objects from list2 indentified in list1.
*/

void refine(coeffx,coeffy,x1,y1,x2,y2,nobj1,nobj2,index,TOL)
     double *coeffx, *coeffy, *x1, *y1, *x2, *y2, TOL;
     int nobj1, nobj2, *index;
{
  int i, j, goodi;
  double xt, yt, x, y, rr, rmin;


  for(i=0;i<nobj1;i++) index[i]=-1;

  for(j=0;j<nobj2;j++) {
    xt = coeffx[0]*x2[j] + coeffx[1]*y2[j] + coeffx[2];
    yt = coeffy[0]*x2[j] + coeffy[1]*y2[j] + coeffy[2];

    rmin = TOL * TOL;
    goodi = -1;
    for(i=0;i<nobj1;i++) 
    {
      x = xt - x1[i];
      y = yt - y1[i];
      rr = x * x + y * y;
      if(rr < rmin)
      {
        goodi = i;
	rmin = rr;
      }
    }

    if (goodi != -1 && index[goodi] == -1)
      index[goodi] = j;
    else if (goodi != -1)
    {
      xt = coeffx[0]*x2[index[goodi]] + coeffx[1]*y2[index[goodi]] + coeffx[2];
      yt = coeffy[0]*x2[index[goodi]] + coeffy[1]*y2[index[goodi]] + coeffy[2];

      x = xt - x1[goodi];
      y = yt - y1[goodi];
      rr = x * x + y * y;
      if(rmin < rr)
        index[goodi] = j;
    }   
  }
  
  return ;

}

