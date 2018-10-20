c-----------------------------------------------------------------------
      subroutine readevec(evec)

      include 'SIZE'
      include 'MOR'

      common /scrk5/ w(lx1*ly1*lz1*lelt)

      real evec(ls,nb)

      call rzero(evec,ls*nb)

      if (nid.eq.0) then
         open (unit=12,file='evectors.dat')
         read (12,*) (evec(k,1),k=1,ls*nb)
         close (unit=12)
      endif

      call gop(evec,w,'+  ',ls*nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine readc0(c0,n)

      include 'SIZE'
      include 'PARALLEL'

      real c0(n)

      if (nid.eq.(np-1)) then
         open (unit=12,file='cten')
         read (12,*) (c0(k),k=1,n)
         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine reada0(a0,n)

      include 'SIZE'

      real a0(n)

      common /scrk5/ w(lx1*ly1*lz1*lelt)

      call rzero(a0,n)

      if (nid.eq.0) then
         open (unit=12,file='amat')
         read (12,*) (a0(k),k=1,n)
         close (unit=12)
      endif

      call gop(a0,w,'+  ',n)

      return
      end
c-----------------------------------------------------------------------
      subroutine readb0(b0,n)

      include 'SIZE'

      real b0(n)

      common /scrk5/ w(lx1*ly1*lz1*lelt)

      call rzero(b0,n)

      if (nid.eq.0) then
         open (unit=12,file='amat')
         read (12,*) (b0(k),k=1,n)
         close (unit=12)
      endif

      call gop(b0,w,'+  ',n)

      return
      end
c-----------------------------------------------------------------------
      subroutine readic(ic,n)

      include 'SIZE'

      real ic(n)

      common /scrk5/ w(lx1*ly1*lz1*lelt)

      call rzero(ic,n)

      if (nid.eq.0) then
         open (unit=12,file='ic')
         read (12,*) (ic(k),k=1,n)
         close (unit=12)
      endif

      call gop(ic,w,'+  ',n)

      return
      end
c-----------------------------------------------------------------------
      subroutine readbasis(ub,vb,wb,nb)

      ! only works for np.eq.1

      include 'SIZE'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ub(lt,0:nb), vb(lt,0:nb), wb(lt,0:nb)

      n=lx1*ly1*lz1*nelt

      if (nid.eq.0) then
         open (unit=12,file='basis')

         do j=0,nb
            read (12,*) (ub(i,j),i=1,n)
         enddo

         do j=0,nb
            read (12,*) (vb(i,j),i=1,n)
         enddo

         if (ldim.eq.3) then
         do j=0,nb
            read (12,*) (wb(i,j),i=1,n)
         enddo
         endif

         close (unit=12)
      endif

      return
      end
c-----------------------------------------------------------------------