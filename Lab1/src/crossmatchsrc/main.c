#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "defs.h"
#include "funcs.h"

/* Program matches two (x,y) cordinate lists */

int main(int argc, char **argv)
{
  double param[6] = {5.0,5.0,0.3,0.01,0.01,0.9};

  char buf[256], *inp1name, *inp2name, *parfname, *outfname;

  double *x1, *y1, *x2, *y2, xc, yc;
  double *mag1, *mag2, *err1, *err2, *bg1, *bg2;
  double *xs1, *ys1, *xs2, *ys2, xmin1, xmax1, ymin1, ymax1;
  double *xm1, *ym1, *xm2, *ym2, xmin2, xmax2, ymin2, ymax2;
  double coeffx[3], coeffy[3], PTOL, LTOL, rr, RTOL;
  int nobj1, nobj2, nsub, nsub1, nsub2, nmatch, *index, max_nsub;
  int i, j, nend, oldend, *flag1, *flag2, niter, maxniter;
  int cutlist, cpy, error = 0, size, cutreg;

  FILE *inpf, *outf;


  /* IO stuff */

  if(argc != 5) {
    printf("usage: %s parameter_file ",argv[0]);
    printf("short_list long_list output_file\n");
    exit(-1);
  }

  parfname = argv[1];
  inp1name = argv[2];
  inp2name = argv[3];
  outfname = argv[4];

  inpf = fopen(parfname,"r");
  fscanf(inpf,"%*s = %d %*s = %lf",&cutreg,&RTOL);
  fscanf(inpf,"%*s = %d %*s = %lf",&cutlist,&LTOL);
  fscanf(inpf,"%*s = %lf %*s = %lf",&param[0],&param[1]);
  fscanf(inpf,"%*s = %lf %*s = %lf",&param[2],&param[3]);
  fscanf(inpf,"%*s = %lf %*s = %lf",&param[4],&param[5]);
  fscanf(inpf,"%*s = %d %*s = %d",&nsub,&max_nsub);
  fscanf(inpf,"%*s = %lf %*s = %d",&PTOL,&maxniter);
  fclose(inpf);

  if((x1   = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((x2   = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((y1   = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((y2   = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((mag1 = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((mag2 = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((err1 = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((err2 = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((bg1  = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;
  if((bg2  = (double *)malloc(INO*sizeof(double))) == NULL) error = 1;

  inpf = fopen(inp1name,"r");
  xmin1 = xmax1 = ymin1 = ymax1 = 0;
  xmin2 = xmax2 = ymin2 = ymax2 = 0;
  nobj1 = 0;
  size = 1;
  while(fgets(buf,256,inpf)!=NULL &&
	sscanf(buf,"%lf %lf %lf %lf %lf", 
        &x1[nobj1],&y1[nobj1],&mag1[nobj1],&err1[nobj1],&bg1[nobj1]) != 0) 
  {
    if(nobj1 == 0)
    {
      xmin1 = xmax1 = x1[nobj1];
      ymin1 = ymax1 = y1[nobj1];
    }
	 
    if (x1[nobj1] < xmin1) xmin1 = x1[nobj1];
    else if(x1[nobj1] > xmax1) xmax1 = x1[nobj1];
    
    if (y1[nobj1] < ymin1) ymin1 = y1[nobj1];
    else if(y1[nobj1] > ymax1) ymax1 = y1[nobj1];

   nobj1++;
    if (nobj1 >= size*INO) 
    {
      size++;
      if((x1  =(double *)realloc(x1,  size*INO*sizeof(double)))==NULL)error=1;
      if((y1  =(double *)realloc(y1,  size*INO*sizeof(double)))==NULL)error=1;
      if((mag1=(double *)realloc(mag1,size*INO*sizeof(double)))==NULL)error=1;
      if((err1=(double *)realloc(err1,size*INO*sizeof(double)))==NULL)error=1;
      if((bg1 =(double *)realloc(bg1, size*INO*sizeof(double)))==NULL)error = 1;
    }
  }
  fclose(inpf);

  inpf = fopen(inp2name,"r");
  nobj2 = 0;
  size = 1;
  while(fgets(buf,256,inpf)!=NULL &&
	sscanf(buf,"%lf %lf %lf %lf %lf",
        &x2[nobj2],&y2[nobj2],&mag2[nobj2],&err2[nobj2],&bg2[nobj2]) 
	!= 0) 
  {
    if(nobj2 == 0)
    {
      xmin2 = xmax2 = x2[nobj2];
      ymin2 = ymax2 = y2[nobj2];
    }
    
    if (x2[nobj2] < xmin2) xmin2 = x2[nobj2];
    else if(x2[nobj2] > xmax2) xmax2 = x2[nobj2];
    
    if (y2[nobj2] < ymin2) ymin2 = y2[nobj2];
    else if(y2[nobj2] > ymax2) ymax2 = y2[nobj2];

    nobj2++;
    if (nobj2 >= size*INO)
    {
      size++;
      if((x2  =(double *)realloc(x2,  size*INO*sizeof(double)))==NULL)error=1;
      if((y2  =(double *)realloc(y2,  size*INO*sizeof(double)))==NULL)error=1;
      if((mag2=(double *)realloc(mag2,size*INO*sizeof(double)))==NULL)error=1;
      if((err2=(double *)realloc(err2,size*INO*sizeof(double)))==NULL)error=1;
      if((bg2 =(double *)realloc(bg2, size*INO*sizeof(double)))==NULL)error=1;
    }
  }
  fclose(inpf);

  if((index = (int *)   malloc(nobj1*sizeof(int)))   == NULL) error = 1;
  if((flag1 = (int *)   malloc(nobj1*sizeof(int)))   == NULL) error = 1;
  if((flag2 = (int *)   malloc(nobj2*sizeof(int)))   == NULL) error = 1;

  if((xs1   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;
  if((ys1   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;
  if((xs2   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;
  if((ys2   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;
  if((xm1   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;
  if((ym1   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;
  if((xm2   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;
  if((ym2   = (double *)malloc(nsub*sizeof(double))) == NULL) error = 1;

  if (error == 1) 
  {
    fprintf(stderr,"Error: not enough memory available\n");
    exit (-1);
  }
  
  /* Select the stars to be matched*/

  /* First selection: all */
  for (i=0;i<nobj1;i++) flag1[i] = 1;
  for (i=0;i<nobj2;i++) flag2[i] = 1;

  /* Second selection: stars of one list around stars of the other list */
  if(cutlist)
    for(i=0;i<nobj2;i++)
    {
      cpy = 0;
      for(j=0;j<nobj1;j++)
      {
        rr = (x1[j]-x2[i])*(x1[j]-x2[i]) + (y1[j]-y2[i])*(y1[j]-y2[i]);
        if(rr < LTOL*LTOL) cpy = 1;
      }	
      if(!cpy) flag2[i] = 0;
    }

  /* Third selection: stars of both lists around the center of the field */

  if(cutreg)
  {
    xc = (xmin1 + xmax1) / 2.0;
    yc = (ymin1 + ymax1) / 2.0; 
    for(i=0;i<nobj1;i++)
    {
      rr = (x1[i]-xc)*(x1[i]-xc) + (y1[i]-yc)*(y1[i]-yc);
      if(rr > RTOL*RTOL) flag1[i] = 0;
    }
    
    xc = (xmin2 + xmax2) / 2.0;
    yc = (ymin2 + ymax2) / 2.0;
    for(i=0;i<nobj2;i++)
    {
      rr = (x2[i]-xc)*(x2[i]-xc) + (y2[i]-yc)*(y2[i]-yc);
      if(rr > RTOL*RTOL) flag2[i] = 0;
    }
  }

  /* Make iterations until all the reference stars appear in the matched list */
  
  nend = 1;
  oldend = 0;
  niter = 0;
  while (nend != oldend && niter < maxniter)
  {
    nmatch = 0;
    niter++;
    oldend = nend;

    /* Make iterations until more than 2 stars have been matched*/
	  
    while (nmatch < 3 && nsub <= max_nsub)
    {
      nmatch = 0;

      if (niter == maxniter) 
      {
        nsub += 10;
        niter = 1;
      }
     
      if((xs1 = (double *)realloc(xs1,nsub*sizeof(double))) == NULL) error = 2;
      if((ys1 = (double *)realloc(ys1,nsub*sizeof(double))) == NULL) error = 2;
      if((xs2 = (double *)realloc(xs2,nsub*sizeof(double))) == NULL) error = 2;
      if((ys2 = (double *)realloc(ys2,nsub*sizeof(double))) == NULL) error = 2;
      if((xm1 = (double *)realloc(xm1,nsub*sizeof(double))) == NULL) error = 2;
      if((ym1 = (double *)realloc(ym1,nsub*sizeof(double))) == NULL) error = 2;
      if((xm2 = (double *)realloc(xm2,nsub*sizeof(double))) == NULL) error = 2;
      if((ym2 = (double *)realloc(ym2,nsub*sizeof(double))) == NULL) error = 2;
     
      if (error == 2) 
      {
        fprintf(stderr,"Error: not enough memory available\n");
        exit (-1);
      }
      
      nsub1 = nsub;
      nsub2 = nsub;
      if(nsub > nobj1) nsub1 = nobj1;
      if(nsub > nobj2) nsub2 = nobj2;

      /* take nsub brightest objects from both lists */

      bright_end(x1,y1,mag1,nobj1,xs1,ys1,&nsub1,flag1);
      bright_end(x2,y2,mag2,nobj2,xs2,ys2,&nsub2,flag2);

      /* match nsub brightest stars for approximate transformation */

      triangles(xs1,ys1,xs2,ys2,nsub1,nsub2,index,param);

      for(i=0;i<nsub1;i++) {
        if(index[i] != -1) {      
          xm1[nmatch] = xs1[i];
          ym1[nmatch] = ys1[i];
          xm2[nmatch] = xs2[index[i]];
          ym2[nmatch] = ys2[index[i]];
          nmatch++;
        }
      } 

      nsub += 10;
    }
    nsub -= 10;
    
    if( nmatch < 3 ) {
      fprintf(stderr,"error: nmatch < 3 !\n");
      exit(-1);
    }

    /* linear fit to nmatch stars indentified by triangles */
  
    xy_lin(xm1,ym1,xm2,ym2,nmatch,coeffx,coeffy);
  
    /* using linear fit transform one list and look for close neighbors */
  
    refine(coeffx,coeffy,x1,y1,x2,y2,nobj1,nobj2,index,PTOL);
  
    nend = 0;
    for(i=0;i<nobj1;i++)
      if(index[i] != -1)
      {
        flag1[i] = 1;
	flag2[index[i]] = 2;
        nend++;
      }
      else
        flag1[i] = 0;

    for(i=0;i<nobj2;i++)
      if (flag2[i] == 2)
        flag2[i] = 1;
      else
        flag2[i] = 0;

    if(nend < 3 && nsub <= max_nsub) nend = oldend-1;
  }

    printf("%d %d\n",(int)(coeffx[2]+0.5),(int)(coeffy[2]+0.5));

  /* Write output file */

  outf = fopen(outfname,"w");
  for(i=0;i<nobj1;i++)
    if(index[i] != -1)
    {
      fprintf(outf,"%9.3lf %10.3lf "  ,x1[i],         y1[i]);
      fprintf(outf,"%10.3lf %10.3lf " ,x2[index[i]],  y2[index[i]]);
      fprintf(outf,"%10.3lf %10.4lf " ,mag1[i],       err1[i]);
      fprintf(outf,"%10.3lf %10.4lf " ,mag2[index[i]],err2[index[i]]);
      fprintf(outf,"%10.3lf %10.3lf\n",bg1[i],        bg2[index[i]]);
    }
  
  fclose(outf);

  free(x1);
  free(y1);
  free(mag1);
  free(err1);
  free(bg1);
  free(x2);
  free(y2);
  free(mag2);
  free(err2);
  free(bg2);
  free(xs1);
  free(xs2);
  free(xm1);
  free(xm2);
  free(ys1);
  free(ys2);
  free(ym1);
  free(ym2);
  free(index);
  free(flag1);
  free(flag2);

  exit(0);
}
