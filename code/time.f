c-----------------------------------------------------------------------
      subroutine rom_step_v

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),rhs(nb)

      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)

c     Variable for vorticity
      real vort(lt,3)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (ad_step.eq.1) then
         step_time = 0.
      endif

      last_time = dnekclock()

      n  = lx1*ly1*lz1*nelt

      icount = min0(ad_step,3)

      call setr_v(rhs,icount)

      if (ad_step.le.3) then
         call cmult2(fluv,bv,ad_beta(1,icount)/ad_dt,nb*nb)
         call add2s2(fluv,av,1/ad_re,nb*nb)
         call lu(fluv,nb,nb,irv,icv)
      endif

      if (isolve.eq.0) then ! standard matrix inversion
         call solve(rhs,fluv,1,nb,nb,irv,icv)
      else if (isolve.eq.1) then ! constrained solve
         !call csolve(rhs,flu,...
      else
         call exitti('incorrect isolve specified...')
      endif

      call shiftu(rhs)
      call setavg
      call setj

      call comp_drag
c     call comp_rms ! old

      step_time=step_time+dnekclock()-last_time

      if (ifheat) call recon(vx,vy,vz,u)

      if (mod(ad_step,ad_iostep).eq.0) then
!        This output is to make sure the ceof matches with matlab code

         if (nio.eq.0) then
            write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid,step_time
            if (ad_step.eq.ad_nsteps) then
               do j=1,nb
                  write(6,*) j,u(j,1),'final'
               enddo
               write (6,*) 'step_time: ',step_time
            else
               do j=1,nb
                  write(6,*) j,u(j,1)
               enddo
            endif
         endif

         if (.true.) then
            idump=ad_step/ad_iostep
            call recon(vx,vy,vz,u)

            ! compute the vorticity of the ROM reconstructed field
            call opcopy(t1,t2,t3,vx,vy,vz)
            call comp_vort3(vort,work1,work2,t1,t2,t3)
            ifto = .true. ! turn on temp in fld file
            call outpost(vx,vy,vz,pavg,vort,'rom')
            if (nio.eq.0) write (6,*) 'inside ifdump'
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rom_step_t

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
      real tmp(0:nb),rhs(nb)

      common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      if (ad_step.eq.1) then
         step_time = 0.
      endif

      last_time = dnekclock()

      n=lx1*ly1*lz1*nelt

      icount = min0(ad_step,3)

      call setr_t(rhs,icount)

      if (ad_step.le.3) then
         call cmult2(flut,bt,ad_beta(1,icount)/ad_dt,nb*nb)
         call add2s2(flut,at,1/ad_re,nb*nb)
         call lu(flut,nb,nb,irt,ict)
      endif

      if (isolve.eq.0) then ! standard matrix inversion
         call solve(rhs,flut,1,nb,nb,irt,ict)
      else if (isolve.eq.1) then ! constrained solve
         !call csolve(rhs,flu,...
      else
         call exitti('incorrect isolve specified...')
      endif

      call shiftt(rhs)

      step_time=step_time+dnekclock()-last_time

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_BDF_coef(ad_alpha,ad_beta)

      real ad_alpha(3,3), ad_beta(4,3)

      call rzero(ad_alpha,3*3)
      call rzero(ad_beta,3*4)

      ad_beta(1,1) = 1.
      ad_beta(2,1) = -1.

      ad_beta(1,2) = 1.5
      ad_beta(2,2) = -2
      ad_beta(3,2) = 0.5

      ad_beta(1,3) = 11./6
      ad_beta(2,3) = -3
      ad_beta(3,3) = 1.5
      ad_beta(4,3) = -1./3.

      ad_alpha(1,1)=1

      ad_alpha(1,2)=2
      ad_alpha(2,2)=-1

      ad_alpha(1,3)=3
      ad_alpha(2,3)=-3
      ad_alpha(3,3)=1

      return
      end
c-----------------------------------------------------------------------
      subroutine evalc(cu,cl,icl,uu)

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

      real cu(nb)

      integer icalld
      save    icalld
      data    icalld /0/

      common /scrk4/ work(lx1*ly1*lz1*lelt)

      real cl(lcloc),uu(0:nb)
      integer icl(3,lcloc)

      if (icalld.eq.0) then
         evalc_time=0.
         icalld=1
      endif

      stime=dnekclock()

      l=1

      call rzero(cu,nb)

      do l=1,ncloc
         i=icl(1,l)
         j=icl(2,l)
         k=icl(3,l)
         cu(i)=cu(i)+cl(l)*uu(j)*u(k,1)
      enddo

      call gop(cu,work,'+  ',nb)

      evalc_time=evalc_time+dnekclock()-stime

      return
      end
c-----------------------------------------------------------------------
c     subroutine rom_const
c DEPRECATED
c This subroutine is solving rom with constrains
c The subroutine is based on BFGS method with barrier function
c-----------------------------------------------------------------------
c     include 'SIZE'
c     include 'TOTAL'
c     include 'MOR'

c     parameter (lt=lx1*ly1*lz1*lelt)

c     Matrices and vectors for advance
c     real tmp(0:nb),tmat(nb,nb+1)
c     real coef(1:nb)

c     common /scrrstep/ t1(lt),t2(lt),t3(lt),work(lt)

c     Variable for vorticity
c     real vort(lt,3)

c     Working arrays for LU

c     common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

c     n  = lx1*ly1*lz1*nelt

c     time=time+ad_dt

c     count = min0(ad_step,3)

c     if (ad_step.le.3) then
c        call cmult2(helm,b,ad_beta(1,count)/ad_dt,nb*nb)
c        call add2s2(helm,a,1/ad_re,nb*nb)
c     endif

c     ONE = 1.
c     ZERO= 0.

c     call mxm(u,nb+1,ad_beta(2,count),3,tmp,1)
c     call mxm(b,nb,tmp(1),nb,opt_rhs(1),1)

c     call cmult(opt_rhs,-1.0/ad_dt,nb+1)

c     s=-1.0/ad_re

c     call add2s2(rhs,av0,s,nb+1) ! not working...
c     do i=0,nb
c        opt_rhs(i)=opt_rhs(i)+s*av0(i,0)
c     enddo
c      call add2s2(rhs,av0(1,0),-1/ad_re,nb)

c     call copy(conv(1,3),conv(1,2),nb)
c     call copy(conv(1,2),conv(1,1),nb)

c     if (np.eq.1) then
c        call mxm(c,nb*(nb+1),u,nb+1,tmat,1)
c        call mxm(tmat,nb,u,nb+1,conv,1)
c     else
c        call evalc(conv)
c     endif

c     call mxm(conv,nb,ad_alpha(1,count),3,tmp(1),1)

c     call sub2(opt_rhs,tmp,nb+1)

c     call copy(u(1,3),u(1,2),nb)
c     call copy(u(1,2),u(1,1),nb)

c     call opt_const

c     if (mod(ad_step,ad_iostep).eq.0) then

c        This output is to make sure the ceof matches with matlab code

c        if (nio.eq.0) then
c           write (6,*)'ad_step:',ad_step,ad_iostep,npp,nid
c           if (ad_step.eq.ad_nsteps) then
c              do j=1,nb
c                 write(6,*) j,u(j,1),'final'
c              enddo
c           else
c              do j=1,nb
c                 write(6,*) j,u(j,1)
c              enddo
c           endif
c        endif

c        if (ifdump) then
c           call opzero(vx,vy,vz)
c           do j=1,nb
c              call opadds(vx,vy,vz,ub(1,j),vb(1,j),wb(1,j),coef(j),n,2)
c           enddo
c           call opadd2  (vx,vy,vz,ub,vb,wb)

c           ! compute the vorticity of the ROM reconstructed field
c           call opcopy(t1,t2,t3,vx,vy,vz)
c           call comp_vort3(vort,work1,work2,t1,t2,t3)
c           ifto = .true. ! turn on temp in fld file
c           call copy(t,vort,n)

c           call outpost (vx,vy,vz,pr,t,'rom')
c        endif
c     endif

c     return
c     end
c-----------------------------------------------------------------------
      subroutine opt_const

      include 'SIZE'
      include 'TOTAL'
      include 'MOR'

c     parameter for barrier function
c     it should start from value greater than one and decrease
      real B_qn(nb,nb)
      real go(nb),fo,qndf
      real tmp(nb,nb),tmp1(nb,nb),tmp2(nb,nb),tmp3(nb,nb)
      real tmp4(nb),tmp5(nb),tmp6(nb,nb),tmp7(nb,nb)
      real yy(nb,nb),ys,sBs
      integer par_step

      if (nio.eq.0) write (6,*) 'inside opt_const'

      par_step = 1
c      par = 1.0 
      par = 0.01 


c     invhelm for computing qnf
c     not changing during BFGS
c     only changes during the time stepper
      if (ad_step.le.3) then 
         call copy(invhelm(1,1),helm(1,1),nb*nb)
         call lu(invhelm,nb,nb,ir,ic)
      endif

c     BFGS method with barrier function starts
      do k=1,par_step

c     use helm from BDF3/EXT3 as intial approximation
         do i=1,nb
            call copy(B_qn(1,i),helm(1,i),nb)
         enddo

         call comp_qnf
         call comp_qngradf

c     compute quasi-Newton step
         do j=1,500

            call copy(tmp3(1,1),B_qn(1,1),nb*nb)
            call lu(tmp3,nb,nb,ir,ic)
            call copy(qns,qngradf,nb)
            call chsign(qns,nb)
            call solve(qns,tmp3,1,nb,nb,ir,ic)
            call add2(u(1,1),qns,nb)

            ! check whether solution exceed the boundary
            do ii=1,nb
               if ((u(ii,1)-sample_max(ii)).ge.1e-10) then
                  u(ii,1) = 0.9 * sample_max(ii)
               elseif ((sample_min(ii)-u(ii,1)).ge.1e-10) then
                  u(ii,1) = 0.9 * sample_min(ii)
               endif
            enddo

            call copy(go,qngradf,nb) ! store old qn-gradf
            call comp_qngradf        ! update qn-gradf
            call sub3(qny,qngradf,go,nb)

            ! update approximate Hessian by two rank-one update
            ! first rank-one update
            ! outer product: s_k * s_k^T               
            call mxm(qns,nb,qns,1,tmp,nb)
             
            ! s_k * s_k^T * B_k
            call mxm(tmp,nb,B_qn,nb,tmp1,nb)

            ! B_k * s_k * s_k^T * B_k 
            call mxm(B_qn,nb,tmp1,nb,tmp2,nb)

            ! s_k^T * B_k * s_k 
            call mxm(B_qn,nb,qns,nb,tmp5,1)
            sBs = glsc2(qns,tmp5,nb)

            ! second rank-one update
            ! outer product: y_k * y_k^T               
            call mxm(qny,nb,qny,1,yy,nb)

            ys = glsc2(qny,qns,nb)

            do ii=1,nb
               call cmult(tmp2(1,ii),-1.0/sBs,nb)
               call cmult(yy(1,ii),1.0/ys,nb)
            enddo

            do ii=1,nb
               call add4(B_qn(1,ii),B_qn(1,ii),tmp2(1,ii),yy(1,ii),nb)
            enddo

            fo = qnf      ! store old qn-f
            call comp_qnf ! update qn-f
            qndf = abs(qnf-fo) 
            write(6,*)'f and old f',qnf,fo,qndf

            if (qndf .lt. 1e-10) goto 900

c     update solution
         enddo
  900    write(6,*)'criterion reached, number of iteration:',j,par 
         par = par*0.1

      enddo
      if (nio.eq.0) write (6,*) 'exitting opt_const'

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_qngradf
      
      include 'SIZE'
      include 'MOR'
      real tmp1(nb),tmp2(nb),tmp3(nb),tmp4(nb)
      real mpar

!      if (nio.eq.0) write (6,*) 'inside com_qngradf'

      call sub3(tmp1,u(1,1),sample_max,nb)  
      call sub3(tmp2,u(1,1),sample_min,nb)  
      call invcol1(tmp1,nb)
      call invcol1(tmp2,nb)

      call add3(tmp3,tmp1,tmp2,nb)

      mpar = -1.0*par

      call add3s12(qngradf,opt_rhs(1),tmp3,-1.0,mpar,nb)

      ONE = 1.
      ZERO= 0.
      call dgemv( 'N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp4,1)
      call add2(qngradf,tmp4,nb)


!      if (nio.eq.0) write (6,*) 'exiting com_qngradf'

      return 
      end
c-----------------------------------------------------------------------
      subroutine comp_qnf
      
      include 'SIZE'
      include 'MOR'
      real tmp1(nb),tmp2(nb),tmp3(nb)
      real tmp4(nb),tmp5(nb),tmp6(nb)
      real term1,term2,term3,term4

c     bar1 and bar2 are the barrier function for two constrains
      real bar1,bar2

!      if (nio.eq.0) write (6,*) 'inside com_qnf'

c     evaluate quasi-newton f

c     term1 represents 0.5*x'*H*x
      ONE = 1.
      ZERO= 0.
      call dgemv( 'N',nb,nb,ONE,helm,nb,u(1,1),1,ZERO,tmp6,1)
      term1 = 0.5 * glsc2(tmp6,u(1,1),nb)
      write(6,*)'term1',term1

c     term2 represents x'*g
      term2 = glsc2(u(1,1),opt_rhs(1),nb)
      write(6,*)'term2',term2

c     term3 represents 0.5*g'*inv(H)*g
      call copy(tmp5,opt_rhs(1),nb)
      call solve(tmp5,invhelm,1,nb,nb,ir,ic)
      term3 = 0.5 * glsc2(opt_rhs(1),tmp5,nb)
      write(6,*)'term3',term3

c     barrier term
      call sub3(tmp1,sample_max,u(1,1),nb)  
      call sub3(tmp2,u(1,1),sample_min,nb)  

c     currently can only come up with this way to compute log for an array
      do i=1,nb
         tmp3(i) = log(tmp1(i))
         tmp4(i) = log(tmp2(i))
      enddo

      bar1 = vlsum(tmp3,nb)
      bar2 = vlsum(tmp4,nb)
      term4 = par*(bar1+bar2)
      write(6,*)'term4',term4

      qnf = term1 - term2 + term3 - term4

!      if (nio.eq.0) write (,*) 'exitting com_qnf'

      return
      end
c-----------------------------------------------------------------------
      subroutine setr_t(rhs,icount)

      include 'SIZE'
      include 'MOR'

      common /scrrhs/ tmp(0:nb)

      real rhs(nb)

      call mxm(ut,nb+1,ad_beta(2,icount),3,tmp,1)
c     call mxm(bv0,nb+1,tmp,nb+1,rhs,1)
      call mxm(bt,nb,tmp(1),nb,rhs,1)

      call cmult(rhs,-1.0/ad_dt,nb)

      s=-1.0/ad_re

c     call add2s2(rhs,av0,s,nb+1) ! not working...
      do i=1,nb
         rhs(i)=rhs(i)+s*at0(i,0)
      enddo

      call copy(ctr(1,3),ctr(1,2),nb)
      call copy(ctr(1,2),ctr(1,1),nb)

      call evalc(ctr,ctl,ictl,ut)

      call mxm(ctr,nb,ad_alpha(1,icount),3,tmp(1),1)

      call sub2(rhs,tmp(1),nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine setr_v(rhs,icount)

      include 'SIZE'
      include 'MOR'

      common /scrrhs/ tmp(0:nb)

      real rhs(nb)

      call mxm(u,nb+1,ad_beta(2,icount),3,tmp,1)
c     call mxm(bv0,nb+1,tmp,nb+1,rhs,1)
      call mxm(bv,nb,tmp(1),nb,rhs,1)

      call cmult(rhs,-1.0/ad_dt,nb)

      s=-1.0/ad_re

c     call add2s2(rhs,av0,s,nb+1) ! not working...
      do i=1,nb
         rhs(i)=rhs(i)+s*av0(i,0)
      enddo

      call copy(cvr(1,3),cvr(1,2),nb)
      call copy(cvr(1,2),cvr(1,1),nb)

      call evalc(cvr,cvl,icvl,u)

      call mxm(cvr,nb,ad_alpha(1,icount),3,tmp(1),1)

      call sub2(rhs,tmp(1),nb)
      if (ifforce) call add2(rhs,bg(1),nb)

      return
      end
c-----------------------------------------------------------------------
      subroutine setavg

      include 'SIZE'
      include 'MOR'
      include 'AVG'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scravg/ ux(lt),uy(lt),uz(lt)

      if (ad_step.eq.1) then
         call rzero(ua,nb+1)
         call rzero(u2a,(nb+1)**2)
      endif

      call add2(ua,u,nb+1)

      do j=0,nb
      do i=0,nb
         u2a(i,j)=u2a(i,j)+u(i,1)*u(j,1)
      enddo
      enddo

      if (ad_step.eq.ad_nsteps) then
         s=1./real(ad_nsteps)
         call cmult(ua,s,nb+1)
         call cmult(u2a,s,(nb+1)**2)

         call recon(ux,uy,uz,ua)
         call outpost(ux,uy,uz,pavg,tavg,'avg')

         call opzero(ux,uy,uz)
         n=lx1*ly1*lz1*nelv
         do j=0,nb
         do i=0,nb
            call col3(ubt,ub(1,i),ub(1,j),n)
            call col3(vbt,vb(1,i),vb(1,j),n)
            if (ldim.eq.3) call col3(wbt,wb(1,i),wb(1,j),n)
            call opadds(ux,uy,uz,ubt,vbt,wbt,u2a(i,j),n,1)
         enddo
         enddo
         call outpost(ux,uy,uz,pavg,tavg,'rms')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setj

      include 'SIZE'
      include 'MOR'

      if (ad_step.eq.3) call copy(uj,u,3*(nb+1))
      if (ad_step.eq.ad_nsteps) then
         call copy(uj(0,4),u,3*(nb+1))
         do k=1,6
         do j=0,nb
         do i=0,nb
            u2j(i,j,k)=uj(i,k)*uj(j,k)
         enddo
         enddo
         enddo
         s=1./real(ad_nsteps)
      endif

      return
      end
c-----------------------------------------------------------------------
