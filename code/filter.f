c-----------------------------------------------------------------------
      subroutine pod_proj(uu,r1,nn,msg)

c     uu: output/input, vector to be filtered
c     r1: input, number of modes to be filtered
c     nn: input, length of vector uu
c     msg: input, used as an indicator for
c     different filter function. Currently, we
c     support msg = step, linear, parabo, and cubic.

      real uu(nn)
      real a1,a2,a3,a4
      integer r1,nn,ncut
      character*6  msg

      if (msg.eq.'step  ') then
         call rzero(uu(nn-r1+1),r1)
      elseif (msg.eq.'linear') then
         ncut = nn-r1
         do i=ncut+1,nn
            uu(i) = (-uu(ncut)/r1)*(i-nn)
         enddo
      elseif (msg.eq.'parabo') then
         ncut = nn-r1
         do i=ncut+1,nn
            uu(i) = (-uu(ncut)/r1**2)*(i-ncut)**2 + uu(ncut)
         enddo
      elseif (msg.eq.'cubic ') then
         ncut = nn-r1
         a1 = 2
         a2 = -3*(ncut+nn)
         a3 = 6*ncut*nn
         a4 = nn**3-3*nn**2*ncut
         do i=ncut+1,nn
            uu(i) = (a1*i**3+a2*i**2+a3*i+a4)*uu(ncut)/(ncut+nn)**3
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pod_df(uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      integer icalld
      save    icalld
      data    icalld /0/

      real tmp(nb)
      real uu(nb)

      if (icalld.eq.0) then
         call setdf(dfops,au,bu,rdft*rdft)
         call dgetrf(nb,nb,dfops,nb,ipiv,info)
         icalld=1
      endif

      if (rdft.le.1e-12) then
      else
         call mxm(bu,nb,uu,nb,tmp,1)
         call dgetrs('N',nb,1,dfops,nb,ipiv,tmp,nb,info)
         call copy(uu,tmp,nb)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setdf(flu,a,b,ad_diff)

      include 'SIZE'
      include 'MOR'


      real flu(nb*nb),a(nb*nb),b(nb*nb)
      real ad_diff
      
      call cmult2(flu,a,ad_diff,nb*nb)
      call add2(flu,b,nb*nb)
         
      return
      end
