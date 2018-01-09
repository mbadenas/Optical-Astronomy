void indexx(int n, double *arrin, int *indx)
{
  int l,j,ir,indxt,i;
  double q;
  
  
  for (j=1;j<=n;j++) indx[j]=j;
  
  l=(n >> 1) + 1;
  ir=n;

  for (;;)
    {
      if (l > 1)
	q=arrin[(indxt=indx[--l])];
      else
	{
	  q=arrin[(indxt=indx[ir])];
	  indx[ir]=indx[1];
	  if (--ir == 1)
	    {
	      indx[1]=indxt;
	      return;
	    }
	}

      i=l;
      j=l << 1;

      while (j <= ir && j>0)
	{
	  if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
	  if (q < arrin[indx[j]])
	    {
	      indx[i]=indx[j];
	      j += (i=j);
	      if(i<1 || i>(n+1) || j<1 || j>(n+1)) break;
	    }
	  else j=ir+1;
	}
      indx[i]=indxt;
    }
}
