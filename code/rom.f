c-----------------------------------------------------------------------
      include 'pod.f'
      include 'read.f'
      include 'aux.f'
      include 'dump.f'
      include 'time.f'
      include 'conv.f'
c-----------------------------------------------------------------------
      subroutine rom_update_v

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      ad_step = istep
      jfield=ifield
      ifield=1

      if (ifheat) then
         if (ifflow) call exitti(
     $   'error: running rom_update_v with ifflow = .true.$',nelv)
         if (istep.eq.0) then
            call rom_setup
         else
            call rom_step
            call recon(vx,vy,vz,u) ! reconstruct velocity to be used in h-t
            time=time+dt
         endif
      else
         call rom_setup

         if (nio.eq.0) write (6,*) 'starting rom_step loop',ad_nsteps

         ad_step = 1
         do i=1,ad_nsteps
            call rom_step
            time=time+dt
            ad_step=ad_step+1
         enddo
      endif

      ifield=jfield

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_setup

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside rom_setup'

      setup_start=dnekclock()

      call rom_init_params
      call rom_init_fields

      if (.not.ifread) then
         call gengram
         call genevec
      endif

      call setbases
      call setops

      if (ifdumpops) call dump_all

      call qoisetup

      setup_end=dnekclock()

      if (nio.eq.0) write (6,*) 'exiting rom_setup'
      if (nio.eq.0) write (6,*) 'setup_time:', setup_end-setup_start

      return
      end
c-----------------------------------------------------------------------
      subroutine setops

      include 'SIZE'

      if (nio.eq.0) write (6,*) 'inside setops'

      call seta
      call setb
      call setc
      call setu

      if (nio.eq.0) write (6,*) 'exiting setops'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_params

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside rom_init_params'

      ad_nsteps=nsteps
      ad_iostep=iostep

      ad_dt = dt
      ad_re = 1/param(2)

      ifl2=.false.
      if (param(33).eq.0) ifl2=.true.

      ifavgic=.false.
      if (param(34).ne.0) ifavgic=.true.

      ifdumpops=.false.
      ifread=.false.
      np35=nint(param(35))
      if (np35.eq.1) then
         ifdumpops=.true.
      else if (np35.eq.2) then
         ifread=.true.
      endif

      ifdrago=.false.
      if (param(36).ne.0) ifdrago=.true.

      ifvort=.false. ! default to false for now
      ifdump=.true.
      ifrecon=.true.
      ifpart=.false.
      ifravg=.false.

      call compute_BDF_coef(ad_alpha,ad_beta)

      if (nio.eq.0) write (6,*) 'exiting rom_init_params'

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_init_fields

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk0/ t1(lt),t2(lt),t3(lt),t4(lt),t5(lt),t6(lt),
     $               u0(lx1*ly1*lz1*lelt,3)

      logical alist

      character*127 fname1

      if (nio.eq.0) write (6,*) 'inside rom_init_fields'

      n=lx1*ly1*lz1*nelv
      call rone(wm1,n)

      ns = ls

      if (ifrecon) then
         call opzero(u0,u0(1,2),u0(1,3))
         call get_saved_fields(us,vs,ws,ns,u0,ifvort)
      endif

      fname1='avg.list'

      if (ifavgic.and..not.ifread) then
         call opcopy(t1,t2,t3,vx,vy,vz)

         inquire (file=fname1,exist=alist)
         if (alist) then
            call auto_averager(fname1)
            call opcopy(ub,vb,wb,vx,vy,vz)
         else
            call opzero(t4,t5,t6)
            s=1./real(ns)
            do i=1,ns
               call opadds(t4,t5,t6,us(1,i),vs(1,i),ws(1,i),s,n,2)
            enddo
            call opcopy(ub,vb,wb,t4,t5,t6)
         endif
         call opcopy(vx,vy,vz,t1,t2,t3)
      else
         call opcopy(ub,vb,wb,vx,vy,vz)
      endif

      if (ifvort) then
         call comp_vort3(u0,t1,t2,ub,vb,wb)
         call comp_vort3(t1,t2,t3,vx,vy,vz)

         call copy(ub,u0,n)
         call rzero(vb,n)
         call rzero(wb,n)

         call copy(vx,t1,n)
         call rzero(vy,n)
         call rzero(vz,n)
      else
         ! copy zero mode to u0(1,1:3)
         call opcopy(u0,u0(1,2),u0(1,3),ub,vb,wb)
      endif

      if (ifrecon) then
         do i=1,ns
            call opsub3(ust(1,i),vst(1,i),wst(1,i),
     $                  us(1,i),vs(1,i),ws(1,i),ub,vb,wb)
         enddo
      endif

      if (nio.eq.0) write (6,*) 'exiting rom_init_fields'

      return
      end
c-----------------------------------------------------------------------
      subroutine setc

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      real cux(lt), cuy(lt), cuz(lt)

      common /scrk1/ t1(lt), binv(lt)
      common /scrcwk/ wk(lcloc)

      call invers2(binv,bm1,lx1*ly1*lz1*nelv)
      call rone(binv,lx1*ly1*lz1*nelv)

      if (nio.eq.0) write (6,*) 'inside setc'

      l=0
      mid = 0
      nlocmin = lcglo/np
      npmin = np-lcglo+(lcglo/np)*np

      if (ifread.and.nid.eq.0) open (unit=12,file='ops/c')

      do k=0,nb
         if (.not.ifread) call setcnv_c(ub(1,k),vb(1,k),wb(1,k))
         do j=0,nb
            if (.not.ifread) call setcnv_u(ub(1,j),vb(1,j),wb(1,j))
            if (.not.ifread) call ccu(cux,cuy,cuz)
            do i=1,nb
               l=l+1
               if (.not.ifread) cltmp(l) = 
     $            op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),cux,cuy,cuz,binv)
               icltmp(1,l) = i
               icltmp(2,l) = j
               icltmp(3,l) = k
               mcloc = nlocmin + mid / npmin
               if (l.eq.mcloc) then
                  if (ifread) then
                     if (nid.eq.0) then
                        read (12,*) (cltmp(kk),kk=1,mcloc)
                     else
                        call rzero(cltmp,mcloc)
                     endif
                  endif
                  if (ifread) call gop(cltmp,wk,'+  ',mcloc)
                  if (nid.eq.mid) then
                     ncloc = mcloc
                     call copy(clocal,cltmp,ncloc)
                     call icopy(icloc,icltmp,ncloc*3)
                  endif
                  mid=mid+1
                  l = 0
               endif
            enddo
         enddo
      enddo

      if (ifread.and.nid.eq.0) close (unit=12)

      if (nio.eq.0) write (6,*) 'exiting setc'

      return
      end
c-----------------------------------------------------------------------
      subroutine seta

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrvh/ h1(lt),h2(lt)
      common /scrns/ usave(lt),vsave(lt),wsave(lt),wk1(lt)
      common /scrread/ tab((nb+1)**2)

      if (nio.eq.0) write (6,*) 'inside seta'

      n=lx1*ly1*lz1*nelt
      call rone (h1,n)
      call rzero(h2,n)

      if (ifread) then
         call read_serial(a0,(nb+1)**2,'ops/a ',wk1,nid)
      else
         do j=0,nb ! Form the A matrix for basis function
            call axhelm(usave,ub(1,j),h1,h2,1,1)
            call axhelm(vsave,vb(1,j),h1,h2,1,1)
            if (ldim.eq.3) call axhelm(wsave,wb(1,j),h1,h2,1,1)
            do i=0,nb
               a0(i,j) = glsc2(ub(1,i),usave,n)+glsc2(vb(1,i),vsave,n)
               if (ldim.eq.3) a0(i,j) = a0(i,j)+glsc2(wb(1,i),wsave,n)
            enddo
         enddo
      endif

      do j=1,nb
      do i=1,nb
         a(i,j)=a0(i,j)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting seta'

      return
      end
c-----------------------------------------------------------------------
      subroutine setb

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrread/ tab((nb+1)**2)

      if (nio.eq.0) write (6,*) 'inside setb'

      if (ifread) then
         call read_serial(b0,(nb+1)**2,'ops/b ',tab,nid)
      else
         do j=0,nb
         do i=0,nb
            b0(i,j) = op_glsc2_wt(ub(1,i),vb(1,i),wb(1,i),
     $                            ub(1,j),vb(1,j),wb(1,j),bm1)
         enddo
         enddo
      endif

      do j=1,nb
      do i=1,nb
         b(i,j)=b0(i,j)
      enddo
      enddo

      if (nio.eq.0) write (6,*) 'exiting setb'

      return
      end
c-----------------------------------------------------------------------
      subroutine setu

      include 'SIZE'
      include 'MOR'

      if (nio.eq.0) write (6,*) 'inside setu'

      call proj2bases(u,uic,vic,wic)

      if (nio.eq.0) write (6,*) 'exiting setu'

      return
      end
c-----------------------------------------------------------------------
      subroutine qoisetup

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrk1/ ux(lt),uy(lt),uz(lt)

      if (ifdrago) then
         call opcopy(ux,uy,uz,vx,vy,vz)

         do i=0,nb
            nio = -1
            call opcopy(vx,vy,vz,ub(1,i),vb(1,i),wb(1,i))
            call torque_calc(1.,x0,.true.,.false.)
            rdgx(i)=dragvx(1)
            rdgy(i)=dragvy(1)
            rdgz(i)=dragvz(1)
            nio = nid
            if (nio.eq.0) then
               write (6,*) i,rdgx(i),rdgy(i),rdgz(i),'dragi'
            endif
         enddo

         call opcopy(vx,vy,vz,ux,uy,uz)
      endif

      return
      end
c-----------------------------------------------------------------------
