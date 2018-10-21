c-----------------------------------------------------------------------
      ! Unit Tests
c-----------------------------------------------------------------------
      include 'rom.f'
c-----------------------------------------------------------------------
      subroutine grammian_test

      include 'SIZE'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real vv(lt,ls)

      call gengram
      call readgram(vv,ls)

      s1=0.
      s2=0.

      do j=1,ls
      do i=1,ls
         s1=s1+0.25*(uu(i,j)-uu(j,i))**2
         s2=s2+uu(i,j)**2
      enddo
      enddo

      esym=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'esym',esym,s1,s2

      s1=0.
      s2=0.

      do i=1,ls*ls
         s1=s1+(uu(i,1)-vv(i,1))**2
         s2=s2+uu(i,1)**2
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.1e-16.and.esym.lt.2-e15) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine eigenvector_test

      include 'SIZE'
      include 'MOR'

      real evec2(ls,nb)

      call readgram(uu,ls)
      call genevec(evec)
      call readevec(evec2,ls,nb)

      s1=0.
      s2=0.

      do j=1,nb
      do i=1,ls
         s1=s1+(evec(i,j)-evec2(i,j))**2
         s2=s2+evec2(i,j)**2
      enddo
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.1e-16) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine bases_test

      include 'SIZE'
      include 'MOR'
      include 'MASS'

      parameter (lt=lx1*ly1*lz1*lelt)

      real ubb(lt,0:nb), vbb(lt,0:nb), wbb(lt,0:nb)
      real du(lt,0:nb), dv(lt,0:nb), dw(lt,0:nb)

      n=lx1*ly1*lz1*nelt

      call readevec(evec,ls,nb)
      call genbases
      call readbases(ubb,vbb,wbb,nb)

      s1=0.
      s2=0.

      do i=0,nb
         call opsub3(du(1,i),dv(1,i),dw(1,i),ub(1,i),vb(1,i),wb(1,i),
     $                                       ubb(1,i),vbb(1,i),wbb(1,i))
         s1=s1+op_glsc2_wt(du(1,i),dv(1,i),dw(1,i),
     $                     du(1,i),dv(1,i),dw(1,i),bm1)
         s2=s2+op_glsc2_wt(ubb(1,i),vbb(1,i),wbb(1,i),
     $                     ubb(1,i),vbb(1,i),wbb(1,i),bm1)
      enddo

      edif=sqrt(s1/s2)
      if (nio.eq.0) write (6,*) 'edif',edif,s1,s2

      iexit=1
      if (edif.lt.1e-16) iexit=0

      call exit(iexit)

      return
      end
c-----------------------------------------------------------------------
      subroutine initial_condition_test

      include 'SIZE'
      include 'MOR'

      real u0(0:nb)

      call readbases(ub,vb,wb,nb)
      call makeic
      call readic(u0,nb+1)

      return
      end
c-----------------------------------------------------------------------
      subroutine a_operator_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
      subroutine b_operator_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
      subroutine c_operator_test(ifl2)

      logical ifl2

      return
      end
c-----------------------------------------------------------------------
