$debug
C    Cross Spectrum ANALYSIS PROGRAM BY Lyu SangJin 1998.6.6
C
C     THIS PROGRAM IS BLACKMAN & TURKEY METHOD 
C     REF. JENKINS & WATTS
C     LAG WINDOW: BATTLET WINDOW(TURKEY WINDOW)
C
C    phase range : 0 ~ 360 degree
      parameter(ix1=50000,ix2=2,ix3=2,mlg=10000,lg=5000)
      dimension mdim(ix2),xm(ix2),ddt(ix2),freq(ix2),
     *nrpp(8),xvar(ix2),peri(ix2)
      common/a/x(ix1,ix2)
      common/b/cov(lg,ix2)
      common/c/pspec(mlg,ix2)
      common/d/w(lg)
      common/com1/nrf
      character*30 fmt
      character*30 fcom
      character*20 infile(2),outfile
c
      data nrpp/13,14,15,16,17,18,19,20/
c
c
      pi=3.141592653
      nrf=3
      ll1=0
      iw=0
c
c     ns;no.data sets(<=ix2)
c
      open(11,file='crosp.con',status='old')
      read(11,'(40x,a20)') infile(1)
      read(11,'(40x,a20)') infile(2)
      read(11,'(40x,a20)') outfile
      open(12,file=outfile)
      
      read(11,'(40x,a30)') fcom
      read(11,'(40x,i10)') ns
      read(11,'(40x,i10)') ic
c
      do 7 j=1,ns
         do 9 i=1,lg
9        cov(i,j)=0.0
         do 8 i=1,mlg
    8    pspec(i,j)=0.0
    7 continue
      do 601 i=1,ns
         nr1=nrpp(i)
         open(nr1,file=infile(i),status='old')
  601 continue
c
c     data read in
c
c     na; no.raw data in each series,
c     lag; no. of lags in correl. cal.( less than lg(250) and
c          <0.1*na)
c     ic=0 ;no cross spect, others;cross spect
c     ndt; data interval used in filtering, dt; raw data time interval
c          if raw spec. is taken, ndt=1
c     infile; input file name( less than 15character)
c
  100 format(/,40x,2i10)
  101 format(40x,3i10)
  102 format(/,40x,a30)
  103 format(40x,3f10.6)
  104 format(40x,i10)
 1044 format(40x,f10.5)
c
c ******main looping*****************
c
      mmax=0
      do 10 j=1,ns
      read(11,102) fmt
      read(11,101) na,ndt,lag
      read(11,103) dt,xm(j),xvar(j)
      read(11,104) infl
      read(11,104) noco
      read(11,104) iw
      read(11,1044) scl
      nr1=nrpp(j)
      do i = 1, na
      read(nr1,fmt) x(i,j)
      enddo
  110 format(//,2x,'*** spectral analysis of the raw data ',/,3x,a30)
c
c      scaling the x(i) values
c

      if(scl.eq.1.0) go to 567
      do 345 i=1,na
  345 x(i,j)=scl*x(i,j)
  567 continue
c
c     do sum,mean removal,calculate varience of each series
c
      x0=0.0
      j1=j
      do 11 i=1,na
   11 x0=x0+x(i,j)
      xm(j)=x0/na
c
      var=0.0
      do 12 i=1,na
      x(i,j)=x(i,j)-xm(j)
   12 var=var+x(i,j)*x(i,j)
      xvar(j)=var/na
c
       ma=na
c
      if(infl.eq.0) then
      write(12,110) fcom
      else
      write(12,109) fcom
          endif
 109  format(//,2x,'** spectral analysis of a filtered data ',/,3x,a30)
      write(12,200) j,infile(j),xm(j),xvar(j),scl,ma
  200 format(//,1x,'station number = ',i5,3x,'station = ',a15,//,3x,
     *' data mean = ',f12.5,3x,'varience = ',f12.5,
     */,3x,'mag. of x scaled by',2x,f10.5,'no.data used = ',i10)
c
      call detrnd(j1,na)
c
c     calculate the autocovariance
c
      do 13 k=1,lag
      coz=0.0
      nk=ma-k
      do 14 i=1,nk
   14 coz=coz+x(i,j)*x(i+k,j)
   13 cov(k,j)=coz/na
      cov(lg,j)=xvar(j)
      if(infl.eq.0) ndt=1
      dnt=dt*ndt
      ddt(j)=dnt
c
      write(12,201) dt,ddt(j)
  201 format(/,3x,'sample interval(dt) =',f10.5,1x,'hrs',/,
     *3x,'dt after filtering = ',f10.5)
c
c     power spectrum
c
      l1=lag
c
      if(l1.ne.ll1) go to 334
      if(iw.eq.iww) go to 335
  334 if(iw.ne.1) go to 333
c
c       hanning window
      do 332 i=1,l1-1
      pim=(pi*i)/float(l1)
  332 w(i)=0.5*(1+cos(pim))
      iww=1
      dgf=(8.*na)/(3.*l1)
      bw=4./(3.*l1*dnt)
      go to 335
  333 continue
c
c       paren window
      l2=l1/2
      do 336 i=1,l2
      ti=float(i)/float(l1)
      td=ti*ti
  336 w(i)=1-6.*td+6.*ti*td
      do 337 i=l2+1,l1-1
      ti=float(i)/float(l1)
 337  w(i)=2.*(1.-ti)**3
      iww=2
      dgf=(3.71*na)/l1
      bw=1.86/(l1*dnt)
  335 continue
      ll1=l1
c
      call pspect(j1,l1,dnt)
c
      xm(j)=float(na)
      mdim(j)=l1
      if(mmax.lt.l1) mmax=l1
c
      write(12,202) dgf,bw,iw,l1
  202 format(/,3x,'degree of freedom = ',f10.5,3x,'bandwidth = ',f10.5
     *,/,3x,'kinds of window= ',i5,2x,'no. of the window-weight = ',i10)
c
c
   10 continue

c
c     writting the cov & power spectrum
c
      i0=0
      p0=0.0
c
      if(noco.eq.0) go to 578
      write(12,203)
 203  format(///,10x,'autocorrelations ',//,2x,'lags',10x,'autocor11',
     *10x,'autocor22')
c
      write(12,204) i0,(cov(lg,j),j=1,ns)
      do 30 i=1,mmax
      write(12,204) i, (cov(i,j),j=1,ns)
   30 continue
  204 format(i6,5(5x,f14.8))
c
  578 continue

c
      write(12,205)
  205 format(///,10x,'power spectrums ',//,2x,'no.',6x,'freq1',
     *2x,'period',2x,'p-spect11',2x,'freq2',2x,'period',2x,'p-spect22')

      write(12,1000) i0,p0,(pspec(mlg,j),j=1,ns)
 1000 format(i5,2x,f7.4,9x,f15.5,16x,f15.5)
      do 40 i=1,mmax
      do 45 j=1,ns
      freq(j)=i/(mdim(j)*2.0)/ddt(j)
      peri(j)=1./freq(j)
   45 continue
   40 write(12,206) i,(freq(j),peri(j),pspec(i,j),j=1,ns)
  206 format(i5,2x,2(f7.4,f9.3,f15.5))
c
c
c     bivariate cross-cov.
c
c
      if(ic.eq.0) go to 20
      do 21 j=1,ic
      read(11,100) j1,j2
c
      ml1=mdim(j1)
      n=ifix(xm(j1))
      nf=mdim(j1)
c
c     na(j1) should be eq. na(j2)
c
      call bcrcov(n,j1,j2,ml1,ks,noco)
      dt=ddt(j1)
      call crspec(ml1,j1,j2,nf,dt)
c
   21 continue
c
   20 continue
      stop
      end
c
c
c
      subroutine detrnd(j,n)
      parameter(ix1=50000,ix2=2)
      common/a/x(ix1,ix2)
c
c     detrending by the l.s.m
c
      tbar=float(n+1)/2.
      zn=float(n)
      smtt=zn*(zn*zn-1.0)/12.0
      sumx=0.0
      do 40 i=1,n
   40 sumx=sumx+x(i,j)*(float(i)-tbar)
      beta=sumx/smtt
      do 50 i=1,n
   50 x(i,j)=x(i,j)-beta*(float(i)-tbar)
      return
      end
c
c
c
      subroutine pspect(j,m,dt)
      parameter(mlg=10000,lg=5000,ix2=2,ix3=2)
      common/d/w(lg)
      common/b/cov(lg,ix2)
      common/c/pspec(mlg,ix2)

      m1=m-1
      pi=3.141592653
      av2=0.0
      do 30 k=1,m1
 30   av2=av2+w(k)*cov(k,j)
c
      clg=cov(lg,j)
      do 20 i=1,m
      f2=(pi*i)/float(m)
      c=cos(f2)
      v1=0.0
      v0=0.0
      do 21 k=m1,1,-1
      v2=2.*c*v1-v0+w(k)*cov(k,j)
      v0=v1
      v1=v2
   21 continue
   20 pspec(i,j)=2.*dt*(clg+2.*(v1*c-v0))
      pspec(mlg,j)=2.0*dt*(clg+2.*av2)
      return
      end
c
c
c
      subroutine bcrcov(n,j1,j2,ml1,ks,noco)
      parameter(ix1=50000,ix2=2,ix3=2,lg=5000,mlg=10000)
      common/a/x(ix1,ix2)
      common/b/cov(lg,ix2)
      common/com1/nrf
c
      ks=0
      cxm=0.0
      cym=0.0
      do 10 k=1,ml1
      cx=0.0
      cy=0.0
      nk=n-k
      do 11 i=1,nk
      cx=cx+x(i,j1)*x(i+k,j2)
  11  cy=cy+x(i+k,j1)*x(i,j2)
      cov(k,1)=cx/n
      cov(k,2)=cy/n
c      if(cov(k,1).gt.cxm) go to 9
c      go to 10
c   9  ks=k
c      cxm=cov(k,1)
   10 continue
c
        do 30 i=1,n
   30   cxm=cxm+x(i,j1)*x(i,j2)
        clg=cxm/n
        cov(lg,2)=clg
c      if(ks.gt.(ml1/10)) ks=0
      if(noco.eq.0) go to 29
      write(12,99) j1,j2,ks
      write(12,100) ks,cov(lg,2),cov(lg,2)
      do 15 k=1,ml1
   15 write(12,100) k,cov(k,1),cov(k,2)
   29 continue
c
      do 12 k=1,ml1
      cov(k,2)=0.5*(cov(k,1)+cov(k,2))
      cov(k,1)=0.5*(cov(k,1)-cov(k,2))
  12  continue
      cov(lg,1)=0.0
c
c      ml1=ml1-ks
c
  100 format(i5,4(5x,f13.5))
   99 format(///,5x,'cross-cov. btwn ',i2 ' and',i2,10x,'droped lag;'
     *,i5,//,2x,'no.',8x,'cr-cov. 12',8x,'cr-cov.21',8x,'---even---',
     *8x,'----odd---')
      return
      end
c
c
c
      subroutine crspec(ml1,j1,j2,nf,dt)
      parameter(mlg=10000,lg=5000,ix2=2,ix3=2)
      real*8 qspec,cospec
      real gain2(mlg)
      common/d/w(lg)
      common/b/cov(lg,ix2)
      common/c/pspec(mlg,ix2)
      dimension cohsq(mlg),phas(mlg)
      common/com1/nrf
c
      epl=1.0e-8
      epu=-1.0e-8
      pi=3.141592653
c
      tco=0.0
      do 10 i=1,ml1-1
   10 tco=tco+cov(i,2)*w(i)
c
      qlg=cov(lg,2)
      deg=360./(2*pi)
      do 20 i=1,nf
      pis=(i*pi)/float(nf)
      c=cos(pis)
      sn=sin(pis)
      v0=0.0
       v1=0.0
      z0=0.0
      z1=0.0
      do 21 k=ml1-1,1,-1
      v2=2.*c*v1-v0+w(k)*cov(k,2)
      z2=2.*c*z1-z0+w(k)*cov(k,1)
      v0=v1
      v1=v2
      z0=z1
      z1=z2
   21 continue
      cospec=2.0*dt*(qlg+2.0*(v1*c-v0))
      qspec=4.0*dt*z1*sn
      sq=cospec*cospec+qspec*qspec
      if(cospec.le.epl.and.cospec.ge.epu) cospec=0.0
      if(cospec.eq.0.0) then
       if(qspec.gt.0.0) then
       ph=90.0
       goto 199
       else
       ph=270.
       goto 199
       endif
      endif
      ph=deg*atan(qspec/cospec)
199   continue       
       if(cospec.gt.0.0) goto 19
          ph=180.+ph
  19  continue
      phas(i)=ph
      if(pspec(i,j1).le.1.0e-8) go to 333
      if(pspec(i,j2).le.1.0e-8) go to 333
      cohsq(i)=sq/(pspec(i,j1)*pspec(i,j2))
      gain2(i)=sq/(pspec(i,j2)**2)
      go to 20
  333 cohsq(i)=33.3
      gain2(i)=33.3
   20 continue
c
c
      sqm=2.*dt*(qlg+2.*tco)
      cohsq(mlg)=sqm/(pspec(mlg,j1)*pspec(mlg,j2))
      gain2(mlg)=sqm/(pspec(mlg,j2)**2)
c
c
      write(12,100) j1,j2
  100 format(///,5x,'cross-spectrum btwn ',i5,2x,'and',i5,2x,
     *//,2x,'no.',6x,'freq.',6x,'period',5x,'coherency',5x,'phase',5x,
     .'Gain relative to 2')
c
      write(12,101) sqrt(abs(cohsq(mlg))),sqrt(abs(gain2(mlg)))
  101 format(4x,'0',11x,'-infinite-',f14.7,'---zero---',f14.7)
      ff2=nf*2.
      do 200 i=1,nf
      fi=float(i)
      peri=dt*ff2/fi
      freq=fi/ff2/dt
      write(12,210) i,freq,peri,sqrt(cohsq(i)),phas(i),sqrt(gain2(i))
  200 continue
 210  format(i5,f11.5,f12.5,f14.7,f10.3,f14.7)
      return
      end