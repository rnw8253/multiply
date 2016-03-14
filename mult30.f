CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Cheat Sheet:
C     a1 - "ATOM"
C     a2 - Atom number
C     a3 - Atom name
C     a4 - Residue name
C     a5 - Residue number
C     a6 - x
C     a7 - y
C     a8 - z
C     a9 - empty
C     a0 - empty
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM multiplier
        integer max
        parameter(max=100000)
        character*66 a
        character*4  a1(max)
        character*7  a2(max)
        character*6  a3(max)
        character*3  a4(max)
        character*6  a5(max)
        character*12 a6(max)
        character*8  a7(max)
        character*8  a8(max)
        character*6  a9(max)
        character*6  a0(max)
        integer i,j,k,ip,npar,alen(10),nnew,inew
        integer ix,iy,iz,imax
        double precision x(3,max),xcom(3),xmax,xnow
        double precision u(3),th,mat(3,3)
        integer nres(max),restot
        double precision x1,x2,r2
        integer inext,inextp,ma(55)
        double precision ran3
C
        open(20,FILE='ran3.dat',STATUS='old')
        read(20,*) inext,inextp
        do i=1,55
          read(20,*) ma(i)
        enddo
        close(20)
C
        nnew=90
        imax=int((dble(nnew)-0.001d0)**(1.d0/3.d0)) - 1
C
        alen(1)=4
        alen(2)=alen(1)+7
        alen(3)=alen(2)+6 
        alen(4)=alen(3)+3 
        alen(5)=alen(4)+6 
        alen(6)=alen(5)+12
        alen(7)=alen(6)+8 
        alen(8)=alen(7)+8 
        alen(9)=alen(8)+6 
        alen(10)=alen(9)+6 
C
        open(20,FILE="mFHF_linear.pdb",STATUS="old")
        i=1
        read(20,'(a66)') a
        do while (a(1:3).ne.'TER')
          a1(i)=a(1:alen(1))
          a2(i)=a(alen(1)+1:alen(2))
          a3(i)=a(alen(2)+1:alen(3))
          a4(i)=a(alen(3)+1:alen(4))
          a5(i)=a(alen(4)+1:alen(5))
          a6(i)=a(alen(5)+1:alen(6))
          a7(i)=a(alen(6)+1:alen(7))
          a8(i)=a(alen(7)+1:alen(8))
          a9(i)=a(alen(8)+1:alen(9))
          a0(i)=a(alen(9)+1:alen(10))
          read(20,'(a66)') a
          i=i+1
        enddo
        npar=i-1
C
        do i=1,npar
          read(a5(i),'(i7)') nres(i)
          read(a6(i),'(f12.3)') x(1,i)
          read(a7(i),'(f8.3)') x(2,i)
          read(a8(i),'(f8.3)') x(3,i)
        enddo
        close(20)
        restot=nres(npar)
C
        do j=1,nnew-1
          do i=1,npar
            ip=i+j*npar
            a1(ip)=a1(i)
            a3(ip)=a3(i)
            a4(ip)=a4(i)
            a9(ip)=a9(i)
            a0(ip)=a0(i)
            nres(ip)=nres(i)+restot*j
          enddo
        enddo
C
        do i=1,3
          xcom(i)=0.d0
        enddo
        do i=1,npar
          do j=1,3
            xcom(j)=xcom(j)+x(j,i)
          enddo
        enddo
        do i=1,3
          xcom(i)=xcom(i)/dble(npar)
        enddo
C
        xmax=0.d0
        do i=1,npar
          do j=i+1,npar
            do k=1,3
              xnow=dabs(x(k,i)-x(k,j))
              if(xnow.gt.xmax) xmax=xnow
            enddo
          enddo
        enddo
        xmax=xmax+0.5d0
        write(0,*) xmax
C        xmax=26.0d0
C
        do i=1,npar
          x(1,i)=x(1,i)-xcom(1)
          x(2,i)=x(2,i)-xcom(2)
          x(3,i)=x(3,i)-xcom(3)
        enddo
C
        ix=1
        iy=0
        iz=0
        do inew=1,nnew-1
          th=2.d0*3.1415926536d0*ran3(ma,inext,inextp)
          x1=1.d0-2.d0*ran3(ma,inext,inextp)
          x2=1.d0-2.d0*ran3(ma,inext,inextp)
          r2=x1*x1+x2*x2
          do while (r2.gt.1.d0)
            x1=1.d0-2.d0*ran3(ma,inext,inextp)
            x2=1.d0-2.d0*ran3(ma,inext,inextp)
            r2=x1*x1+x2*x2
          enddo
          u(1)=2.d0*x1*sqrt(1.d0-r2)
          u(2)=2.d0*x2*sqrt(1.d0-r2)
          u(3)=1.d0-2.d0*r2
          do i=1,3
            do j=1,3
              mat(j,i)=u(i)*u(j)*(1.d0-dcos(th))
              if(i.eq.j) then
                mat(j,i)=mat(j,i)+dcos(th)
              else
                k=i+1
                if(k.eq.4) k=1
                if(j.eq.k) then
                  k=k+1
                  if(k.eq.4) k=1
                  mat(j,i)=mat(j,i)-u(k)*dsin(th)
                else
                  mat(j,i)=mat(j,i)+u(k)*dsin(th)
                endif
              endif
            enddo
          enddo
C
          do i=1,npar
            ip=i+inew*npar
            do j=1,3
              x(j,ip)=mat(j,1)*x(1,i)+mat(j,2)*x(2,i)+mat(j,3)*x(3,i)
            enddo
            x(1,ip)=x(1,ip)+xmax*dble(ix)
            x(2,ip)=x(2,ip)+xmax*dble(iy)
            x(3,ip)=x(3,ip)+xmax*dble(iz)
            do j=1,3
              x(j,ip)=x(j,ip)+xcom(j)
            enddo
          enddo
C
          if(ix.eq.imax) then
            ix=0
            if(iy.eq.imax) then
              iy=0
              iz=iz+1
            else
              iy=iy+1
            endif
          else
            ix=ix+1
          endif
        enddo
C
        do i=1,npar
          do j=1,3
              x(j,i)=x(j,i)+xcom(j)
          enddo
        enddo
C
        do j=1,nnew
          do i=1,npar
             ip=i+(j-1)*npar
             write(6,998)  a1(i),     ip ,  a3(i),a4(i),nres(ip),
     &                  x(1,ip),x(2,ip),x(3,ip),a9(i),a0(i)
          enddo
          write(6,'(a6)') 'TER   '
        enddo
998     format(a4,i7,a6,a3,i6,f12.3,f8.3,f8.3,a6,a6)
        write(6,'(a6)') 'END   '
C
        open(20,FILE='ran3.dat',STATUS='unknown')
        write(20,*) inext,inextp
        do i=1,55
          write(20,*) ma(i)
        enddo
        close(20)
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function ran3(ma,inext,inextp)
        integer mbig,mseed,mz,mbig1,inext,inextp,mj
        double precision fac
        parameter(mbig=1000000000,mseed=161803398,mz=0)
        parameter(mbig1=mbig-1,fac=1.d0/mbig)
        integer ma(55)
C
        inext=inext+1
        if(inext.eq.56) inext=1
        inextp=inextp+1
        if(inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)
        if(mj.le.mz)mj=mj+mbig1
        ma(inext)=mj
        ran3=mj*fac
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

