C -----------------------------------------------------------------------
C  Program to apply transformation coefficients to a table of 
C  instrumental photometry
C
C  ISYA 2004                               (c) Ignasi Ribas (IEEC, Spain)
C -----------------------------------------------------------------------
C   TODO: Possibility of including positional data to keep track
C
      PROGRAM TRANSSTD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      Open(1,file='coefficients.dat',status='old')
      Read(1,*)
      Read(1,*)
      Read(1,*)
      Read(1,100) coey1,coey2
      Read(1,100) ecoy1,ecoy2
      Read(1,*)
      Read(1,*)
      Read(1,*)
      Read(1,*)
      Read(1,100) coeb1,coeb2
      Read(1,100) ecob1,ecob2
      Close(1)
 100  Format(18x,f7.3,4x,f7.3)

      Open(1,file='stars.lst',status='old')
      Open(2,file='stars_std_pos.lst')
      Write(2,150)
 150  Format('# V_std ','  s_V  ',' (b-y)_s','  s_by ')
  1   Read(1,*,end=2) xceny,yceny,xcenb,ycenb,yinst,ery,binst,erb
      ystd=yinst+coey1+coey2*(binst-yinst)
      bystd=coeb1+coeb2*(binst-yinst)
      erby=dsqrt(erb*erb+ery*ery)
      erys=dsqrt(ery*ery+ecoy1*ecoy1+coey2*coey2*erby*erby+
     $ (binst-yinst)*(binst-yinst)*ecoy2*ecoy2)
      erbys=dsqrt(ecob1*ecob1+coeb2*coeb2*erby*erby+
     $ (binst-yinst)*(binst-yinst)*ecob2*ecob2)
      Write(2,200) ystd,erys,bystd,erbys,xceny,yceny
 200  Format(2(f8.3,f7.3),2f9.2)
      Goto 1
  2   Close(1)

      Stop
      End
