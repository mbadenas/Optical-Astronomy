/*
   Given (x1,y1) coordinates of nobj objects and (x2,y2) coordinates
   of the same objects on a different grid calculate first order
   coordinate transformation between two frames. coeffx and coeffy
   contain 3 coefficients each (in order x, y and shift).
*/

void xy_lin(x1,y1,x2,y2,nobj,coeffx,coeffy)
     double *x1, *x2, *y1, *y2, *coeffx, *coeffy;
     int nobj;
{
  double sum_x2x2=0.0, sum_x2y2=0.0, sum_y2y2=0.0;
  double sum_x2=0.0, sum_y2=0.0, sum_x1x2=0.0, sum_x1y2=0.0;
  double sum_x1=0.0, sum_y1=0.0, sum_y1x2=0.0, sum_y1y2=0.0;
  double a11, a12, a21, a22, b1, b2;
  int i;


  for(i=0;i<nobj;i++) {
    sum_x2x2 += x2[i]*x2[i];
    sum_x2y2 += x2[i]*y2[i];
    sum_y2y2 += y2[i]*y2[i];

    sum_x1x2 += x1[i]*x2[i];
    sum_y1y2 += y1[i]*y2[i];
    sum_x1y2 += x1[i]*y2[i];
    sum_y1x2 += y1[i]*x2[i];
    
    sum_x1 += x1[i];
    sum_x2 += x2[i];
    sum_y1 += y1[i];
    sum_y2 += y2[i];
  }
  

  a11 = sum_x2x2 - sum_x2*sum_x2/nobj;
  a12 = sum_x2y2 - sum_x2*sum_y2/nobj;
  a21 = sum_x2y2 - sum_x2*sum_y2/nobj;
  a22 = sum_y2y2 - sum_y2*sum_y2/nobj;


  b1  = sum_x1x2 - sum_x1*sum_x2/nobj;
  b2  = sum_x1y2 - sum_x1*sum_y2/nobj;

  coeffx[0] = ( b1/a12 - b2/a22 ) / ( a11/a12 - a21/a22 );
  coeffx[1] = ( b1/a12 - coeffx[0]*a11/a12 );
  coeffx[2] = ( sum_x1 - coeffx[0]*sum_x2 - coeffx[1]*sum_y2 )/nobj;

  b1  = sum_y1x2 - sum_y1*sum_x2/nobj;
  b2  = sum_y1y2 - sum_y1*sum_y2/nobj;

  coeffy[0] = ( b1/a12 - b2/a22 ) / ( a11/a12 - a21/a22 );
  coeffy[1] = ( b1/a12 - coeffy[0]*a11/a12 );
  coeffy[2] = ( sum_y1 - coeffy[0]*sum_x2 - coeffy[1]*sum_y2 )/nobj;

  return ;
}
