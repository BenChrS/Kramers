c----------------------------------------------------------
c                                                          
c             Die BROWNsche Bewegung                       
c            ========================                      
c                                                          
c  Es wird in diesem Programm folgende Gleichung geloest.  
c                                                          
c    m[D2]X(t)+[In(0-t)](gamma(t-t')*[D1]X(t'))=R(t)       
c                                                          
c     mit   X(t=0)=0 und V(t=0)=0                          
c                                                          
c  r(t) ist nichtweisses Rauschen.                         
c  <R(t)R(t')>=I(t-t')=A*Ie(t-t'),wobei Ie auf 1           
c  normiert ist. Wir haben zwei Ie gewaehlt.               
c  Ie(t)=1/(2*tau)*exp(-|t|/tau)                           
c  Ie(t)=1/(a*sqrt(pi))*exp(-(t/a)^2)                      
c  Durch Umformung loest man eine neue Gleiung:            
c                                                          
c    [D2]x(t)+Q[In(0--t)](Ie(t-t')*[D1]x(t'))=r(t)  (*)    
c                                                   ===    
c       X(t)=sqrt(T/m)*x(t)                                
c   <r(t)r(t')>=Q*Ie(t-t'), Q=A/(m*T)                      
c                                                          
c  Um die DGL zu loesen,benutzen wir das Adams-Bashforth   
c  Verfahren. Das Integration rechnen wir mit der Gross-   
c  trapezformel.                                           
c  Die entsprechende GREENsche Funktion und die Spektral-  
c  funktion koennen auch gerechnet werden, um die Eigen-   
c  schaften des Systems zu erhalten.                       
c  Bei der FFT ist das stepscale df. Wir koennen df so     
c  klein waehlen,dass dw=2*pi*df die gleiche Groesseord-   
c  nung wie d hat.                                         
c                                                          
c  file1='corrl.cnt'   : die Korrelationsfunktion          
c  file2='puls.cnt'    : der einzelne Puls                 
c  file3='spect.cnt'   : die Spektralfunktion              
c  file4='x2.cnt'      : <X*X>                             
c  file5='kinetic.cnt' : 1/2*m*<V*V>                       
c  file8='x.cnt'       : X(t)                              
c  file9='v.cnt'       : V(t)                              
c file10='noise.cnt'   : nichtweisses Rauschen             
c  green=1/0           : GREENsche Funktion rechnen ja/nein
c file11='greenx.cnt'  : GREENsche Funktion fuer X         
c file12='greenv.cnt'  : GREENsche Funktion fuer V         
c                                                          
c----------------------------------------------------------

      PROGRAM brown

CU    USES corf,puls,spekt,noise

      INTEGER nt,nv,nn,it,nwd,nw
      REAL T,temp,A,masse,pi,eps,x0,v0
      PARAMETER (T=30., temp=1., A=4., masse=1.)
      PARAMETER (pi=3.1415926,eps=0., x0=0., v0=0.)
      PARAMETER (nt=30*2**7, nv=20*2**7, nn=2**13)
      PARAMETER (it=500, nwd=2**4, nw=2*2**6*nwd)
      INTEGER nb,idg
      REAL d,df,dw,Q,tm,c,ck,kin,tt,ss
      REAL r(0:nt),f(nn+1),fp(2*nv+1),gamma(nn+1)
      REAL x(0:nt),y(0:nt),v(0:nt),x2(0:nt),v2(0:nt),xg(0:nt),vg(0:nt)
      REAL relr(0:nt),relv(nt-nt/2),s(0:nw)

      open(1,status='unknown',file='corrl7.cnt')
      open(2,status='unknown',file='puls7.cnt')
      open(3,status='unknown',file='spect7.cnt')
      open(4,status='unknown',file='x27.cnt')
      open(5,status='unknown',file='kinetic7.cnt')
      open(8,status='unknown',file='x7.cnt')
      open(9,status='unknown',file='v7.cnt')
      open(10,status='unknown',file='noise7.cnt')
      open(11,status='unknown',file='greenx7.cnt')
      open(12,status='unknown',file='greenv7.cnt')
      open(13,status='unknown',file='corrlf7.cnt')
      open(14,status='unknown',file='corrlvf7.cnt')

c=====Konstante=====================
 
      d=T/nt
      idg=10
      Q=A/(temp*masse)
      tm=temp/masse
 
c=====Null Anweisung==========

      do 10 l=0,nt
        v2(l)=0.
        x2(l)=0.
        relr(l)=0.
 10   continue

      do 20 j=1,nt-nt/2
         relv(j)=0.
 20   continue  
     
c=====Spektralfunktion==============

      df=2*d/nwd
c      call spekt(nwd*nn,tm,Q,d,df,nw,s)
      dw=2.*pi*df
      tt=0.
      ss=0.
      do 30 i=0,nw
        ss=ss+s(i)
        write(3,*) tt, s(i)
        tt=tt+dw
 30   continue
      ss=(ss-0.5*s(0))*dw/(2.*pi)/tm
      write(*,*)ss

c=====Rechenzeitsparen==============
c
c      call corf(d,nn,f)
c      nb=nn+1
c      do 40 i=1,nn+1
c        if(f(i).le.-0.01) then
c           nb=i
c           i=nn+1
c        end if
c 40   continue
c      do 45 j=nb,nn+1
c         nb=nn+1
c         if(abs(f(j)).le.1.E-4) then
c            nb=j
c            j=nn+1
c         end if
c 45   continue
c      write(*,*)'nb=',nb
     
c=====Die GREENsche Funktion=====

      x0g=0.
      v0g=1.
      do 50 i=1,2*nv+1
         fp(i)=0.
 50   continue

c      call solve(Q,f,nb,x0g,v0g,d,fp,nv,idg,r,xg,vg,nt)

c=====einzelner Puls================

      call corf(d,nn,f)
      tt=0.
      do 60 i=1,nn+1
         gamma(i)=f(i)
         write(1,*)tt,gamma(i)
         tt=tt+d
 60   continue 
      call puls(f,nn,nv,d,eps,fp)
      tt=0.
      do 70 l=1,2*nv+1
         write(2,*)tt,fp(l)
         tt=tt+d
 70   continue
      
c=====Hauptschleife===============
      
      idg=10
      c=sqrt(tm)
      do 80 i=1,it         
         call solve(Q,gamma,nb,x0/c,v0/c,d,fp,nv,idg,r,x,v,nt)
         do 90 j=0,nt
            x2(j)=x2(j)+x(j)*x(j)
            v2(j)=v2(j)+v(j)*v(j)
            relr(j)=relr(j)+r(0)*r(j)
 90      continue
         do 100 l=nt,nt/2+1,-1
            relv(nt+1-l)=relv(nt+1-l)+v(nt)*v(l)
 100     continue

         if(i.eq.it/2) then
            do 110 j=0,nt
               y(j)=x(j)
 110        continue
         end if
         write(*,*)i
 80   continue

c=========Output==============

      c=tm/it
      ck=0.5*masse

      tt=0.
      do 200 l=0,nt
         x(l)=sqrt(tm)*x(l)
         y(l)=sqrt(tm)*y(l)
         v(l)=sqrt(tm)*v(l)
         x2(l)=c*x2(l)
         v2(l)=c*v2(l)
         kin=ck*v2(l)
c         kin=v2(l)
         relr(l)=relr(l)/it
         write(4,*) tt,x2(l)
         write(5,*) tt,kin
         write(8,*) tt,x(l),y(l)
         write(9,*) tt,v(l)
         write(10,*) tt,r(l)
         write(11,*) tt,xg(l)
         write(12,*) tt,vg(l)
         write(13,*) tt,relr(l)
         tt=tt+d
200    continue
      
       tt=0.
       do 210 i=1,nt-nt/2
          relv(i)=tm*relv(i)/it
          write(14,*) tt,relv(i)
          tt=tt+d
210    continue

       tt=3.*d
       do 220 j=3,nt
          write(15,*)alog(tt),alog(x2(j)),
     &        (alog(x2(j))-alog(x2(j-1)))/alog(tt/(tt-d))
          tt=tt+d
220    continue

      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      close(8)
      close(9)
      close(10)
      close(11)
      close(12)
      close(13)
      close(14)

      STOP
      END


c===================subroutine=====================
      
c-----------------------------
c    correlationsfuntion Ie   
c-----------------------------
      SUBROUTINE corf(d,nn,f)
      INTEGER nn
      REAL x,d,a,tau
      PARAMETER(pi=3.1415926, a=1., tau=2)
      REAL f(nn+1)
      x=0.
      do 10 i=1,nn+1
c        f(i)=1./(a*sqrt(pi))*exp(-(x/a)**2)
c        f(i)=0.5/tau*exp(-x/tau)
c        f(i)=(1.0-a*a/2.*x*x)*exp(-a*x)
        f(i)=0.25*a*a*(1.0-a*x)*exp(-a*x)
        x=x+d
10    continue
      return
      END


c----------------------------------------------------------
c                 solve the equation                       
c                                                          
c [D2]x(t)+Q[In(0->t)](gamma(t-t')*[D1]x(t'))=r(t)  (*)    
c                                                   ===    
c   <r(t)>=0,    <r(t)r(t')>=Q*gamma(t-t').                
c----------------------------------------------------------
      SUBROUTINE solve(Q,gamma,nb,x0,v0,d,fp,nv,idg,r,x,v,nt)
CU    USERS noise
      INTEGER nb,nv,nt,idg
      REAL Q,x0,v0,d
      REAL gamma(nb),fp(2*nv+1)
      INTEGER jb
      REAL c1,c2,g1,g2,g3,z1,z2,z3,b1,b2,b3,b4,b5
      REAL s1,s2,s3
      REAL coef(3),r(0:nt),x(0:nt),v(0:nt)

c=====Anfangsbedingungen==========

      x(0)=x0 
      v(0)=v0

c=====color noise=================

      call noise(fp,d,nt,nv,idg,r)
     
c=====Koeffizienten==============

      coef(1)=23./12.
      coef(2)=-16./12.
      coef(3)=5./12.

      c1=-Q*d**2
      c2=d*sqrt(Q)

      g1=coef(1)*c1
      g2=coef(2)*c1
      g3=coef(3)*c1

      z1=coef(1)*c2
      z2=coef(2)*c2
      z3=coef(3)*c2

      b1=1.5*d
      b2=-0.5*d
      b3=1.5*c1
      b4=1.5*c2
      b5=-0.5*c2

c=====Adams-Bashforth-Verfahren===
      
      s3=0.
      x(1)=x(0)+d*v(0)
      v(1)=v(0)+c2*r(0)
      s2=0.5*(gamma(1)*v(1)+gamma(2)*v(0))
      x(2)=x(1)+b1*v(1)+b2*v(0)
      v(2)=v(1)+b3*s2+b4*r(1)+b5*r(0)

      do 10 j=3,nt
         jb=min0(j,nb)
         s1=0.5*(gamma(1)*v(j-1)+gamma(jb)*v(j-jb))
         do 20 l=2,jb-1
            s1=s1+gamma(l)*v(j-l)
20       continue
         x(j)=x(j-1)+d*(coef(1)*v(j-1)+coef(2)*v(j-2)+coef(3)*v(j-3))
         v(j)=v(j-1)+g1*s1+g2*s2+g3*s3+z1*r(j-1)+z2*r(j-2)+z3*r(j-3)
         s3=s2
         s2=s1
10    continue

      return
      END

      
c------------------------------
c      noise r(E)              
c------------------------------
      SUBROUTINE noise(fs,d,np,nv,idg,r)
CU    USERS wn
      INTEGER np,nv,nd,idg
      REAL d,c
      REAL fs(2*nv+1),p(0:np+2*nv),r(0:np)
      nd=np+2*nv
      c=sqrt(d)
      call wn(p,nd,idg)
      do 10 i=0,np
        r(i)=p(i)*fs(1)
        do 20 j=2,2*nv+1
          r(i)=r(i)+p(i+j-1)*fs(j)
20      continue
        r(i)=c*r(i)
10    continue
      return
      END
      

c------------------------------
c     white noise              
c------------------------------
      SUBROUTINE wn(p,np,idg) 
CU    USES gasdev      
      INTEGER np,idg
      REAL p(0:np)
      do 10 l=0,np
          p(l)=gasdev(idg)
10    continue
      return
      END


c------------------------------------------------------
c               Spektralfunktion                       
c                                                      
c   n=2*nn,da man von Null bis Unendlichen integriert. 
c   Der erste Punkt ist deswegen f(1)/2,da hier Trapez-
c   formel angewandet ist. ff ist Frequenz.            
c   Umfangfrequenz w=2*pi*ff                           
c   Deswegen ist dw=2*pi*df                            
c------------------------------------------------------
      SUBROUTINE spekt(n,tm,Q,d,df,nw,s)
CU    USERS corf,realft1
      INTEGER n,nw
      REAL pi
      PARAMETER(pi=3.1415926)
      REAL tm,Q,d,ff,real,real2,imag
      REAL f(n),s(0:nw)
      call corf(d,n-1,f)
      f(1)=f(1)/2.
      call realft(f,n)
      f(2)=0.
      ff=df*(2.*pi)
c      ff=0.
      do 20 i=2,nw+1
        real=d*Q*f(2*i-1)
        imag=d*Q*f(2*i)
        write(21,*)ff,real/Q,imag/Q
        real2=real*real
        s(i-1)=2.*tm*(real/(real2+(imag-ff)**2))
        ff=ff+df*(2.*pi)
20    continue
c      s(0)=2.*tm*5.*Q/(2.*Q+1)**2
      return
      END


c-----------------------------
c         one puls            
c-----------------------------
      SUBROUTINE puls(b,n,nv,d,eps,bs)
CU    USERS cosft1
      INTEGER n,nv
      REAL b(n+1),bs(2*nv+1)
      REAL d,eps
      call cosft1(b,n)
      do 10 i=1,n+1
        b(i)=abs(b(i))
        if (b(i).lt.eps) then 
          b(i)=0.
         else 
          b(i)=sqrt((b(i)))
         end if
10    continue
      call cosft1(b,n)
      do 20 i=nv+1,1,-1
        bs(nv+i)=(2*d)**1.5*b(i)
20    continue 
      do 30 l=nv,1,-1
        bs(l)=bs(2*nv-l+2)
30    continue
      return
      END

       
c-------------------------------------
c        Cosinetransformation         
c-------------------------------------
      SUBROUTINE cosft1(y,n)
      INTEGER n
      REAL y(n+1)
CU    USES realft
      INTEGER j
      REAL sum,y1,y2
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/n
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      sum=0.5*(y(1)-y(n+1))
      y(1)=0.5*(y(1)+y(n+1))
      do 11 j=1,n/2-1
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=0.5*(y(j+1)+y(n-j+1))
        y2=(y(j+1)-y(n-j+1))
        y(j+1)=y1-wi*y2
        y(n-j+1)=y1+wi*y2
        sum=sum+wr*y2
11    continue
      call realft(y,n)
      y(n+1)=y(2)
      y(2)=sum
      do 12 j=4,n,2
        sum=sum+y(j)
        y(j)=sum
12    continue
      return
      END
      
c---------------------------------------------
c Fourietransformation einer reellen Funktion 
c---------------------------------------------
      SUBROUTINE realft(data,n)
      INTEGER n
      REAL data(n)
CU    USES four1
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      c2=-0.5
      call four1(data,n/2)
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      h1r=data(1)
      data(1)=h1r+data(2)
      data(2)=h1r-data(2)
      return
      END
 
c------------------------------------
c  Fast Fourietransformation (FFT)   
c------------------------------------
      SUBROUTINE four1(data,nn)
      INTEGER nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/mmax
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END 

c---------------------------------------------
c   eine zufaellige Zahl zwischen 0 und 1     
c---------------------------------------------
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
     
c------------------------------------------
c        Gauss-Verteilung                  
c------------------------------------------
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
      

















