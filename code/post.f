c-----------------------------------------------------------------------
      subroutine drivep

      include 'POST'
      include 'MASS'
      include 'GEOM'
      include 'PARALLEL'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      character*127 flist

c     nelp=512
      nelp=47
c     nelp=32
c     nelp=1
      nsnap=ns

      ns=ls
      call rflist(fnames,ns)

c     iftherm=.false.
      iftherm=.true.

      call ilgls_setup(ilgls,msr,ms,ns,np,nid)
      call iglls_setup(iglls,itmp,ms,ns,np,nid)

      call ilgls_setup(ilglsp,msrp,msp,nb+1,min(np,nb+1),nid)

      nsmax=ivlmax(ms,np)
      call shift_setup(igsh,nekcomm,itmp,nsmax*nelp*lxyz*ldim,np)

      do id=0,np-1
         if (id.eq.nid) then
            do ip=1,np
               write (6,*) id,ip,ms(ip),'mip'
            enddo
         endif
         call nekgsync()
      enddo

      call sleep(1)

      inel=1
      ieg1=0

      ns=ls
      nsg=lsg

      call rzero_ops

      ng=(ndim-1)*3

      do while (ieg1+1.le.nelgv)
         ieg0=ieg1+1
         ieg1=min(ieg1+inel+nelp-1,nelgv)
         write (6,*) nid,'working on elements ',ieg0,ieg1
         nel=ieg1-ieg0+1
         n=lxyz*(ieg1-ieg0+1)
         call rsnapsm(uu,tt,ieg0,ieg1)

         call setgeom(gfac,w9,ieg0,ieg1,lxyz,ng,nid)
         call setvisc(visc,w,ieg0,ieg1,lxyz,nid)
         call setmass(mass,wv1,ieg0,ieg1,lxyz)
         call setrxp(rxp,rxpt,ieg0,ieg1)

         call setbb(gub,uu,mass,wvf1,wvf2,wvf12,ilgls(1),ms,n,ndim,igsh)
         call setaa(gua,uu,visc,gfac,wvf1,wvf2,wvf12,ilgls(1),
     $      ms,n,nel,ndim,ng,igsh)
c        call setcc(guc,uu,uu,rxp,wvf1,wvf2,wvf3,wvf4,ilgls(1),
c    $      ms,msr,n,ndim,ndim,nel,igsh)

         if (iftherm) then
            call setbb(gtb,tt,mass,wvf1,wvf2,wvf12,ilgls(1),ms,n,1,igsh)
            call setaa(gta,tt,visc,gfac,wvf1,wvf2,wvf12,ilgls(1),
     $         ms,n,nel,1,ng,igsh)
c           call setcc(gtc,uu,tt,rxp,wvf1,wvf2,wvf3,wvf4,ilgls(1),
c    $         ms,msr,n,ndim,1,nel,igsh)
         endif
      enddo

      call setcc_snap(gc2)

      call dump_parallel(gub,ms(nid+1)*ns,'ops/gub ',nid)
      call dump_parallel(gua,ms(nid+1)*ns,'ops/gua ',nid)
      call dump_parallel(guc,nl,'ops/guc ',nid)

      if (iftherm) then
         call dump_parallel(gtb,ms(nid+1)*ns,'ops/gtb ',nid)
         call dump_parallel(gta,ms(nid+1)*ns,'ops/gta ',nid)
         call dump_parallel(gtc,nl,'ops/gtc ',nid)
      endif

      if (np.eq.1)
     $   call dump_parallel(guc2,ms(nid+1)*ns*ns,'ops/gc2 ',nid)

      call setgg(gram,gub,eevec,ilgls,nsg,ms(mid+1))
      call setqq(qu(nsg+1),eevec,gram,nsg)

      call dump_serial(gram,nsg*nsg,'ops/gu_ ',nid)

      call rzero(qu,nsg)

      call dump_serial(qu,nsg*nsg,'ops/qu ',nid)

      mmm=(ns+1)*ns

      call gsub0(gub,bbt,aat,ns,nsg)
      call dump_parallel(gub,ms(nid+1)*ns,'ops/gub0 ',nid)

      call gsub0(gua,bbt,aat,ns,nsg)
      call dump_parallel(gua,ms(nid+1)*ns,'ops/gua0 ',nid)

      call setgg(gram,gub,eevec,ilgls,nsg,ms(mid+1))
      call setqq(qu0(nsg+1),eevec,gram,nsg)

      call cfill(qu0,1./nsg,nsg)
      do i=1,ns
         s=-(1./nsg)*vlsum(qu0(1+i*nsg),nsg)
         call cadd(qu0(1+i*nsg),s,nsg)
      enddo

      call dump_serial(gram,nsg*nsg,'ops/gu0_ ',nid)

      if (iftherm) then
         call setgg(gram,gtb,eevec,ilgls,nsg,ms(mid+1))
         call setqq(qt(nsg+1),eevec,gram,nsg)

         call dump_serial(gram,nsg*nsg,'ops/gt_ ',nid)

         call rzero(qt,nsg)

c        call dump_serial(qt,nsg*nsg,'ops/qt ',nid)

         mmm=(ns+1)*ns

         call gsub0(gtb,bbt,aat,ns,nsg)
         call dump_parallel(gtb,ms(nid+1)*ns,'ops/gtb0 ',nid)

         call gsub0(gta,bbt,aat,ns,nsg)
         call dump_parallel(gta,ms(nid+1)*ns,'ops/gta0 ',nid)

         call setgg(gram,gtb,eevec,ilgls,nsg,ms(mid+1))
         call setqq(qt0(nsg+1),eevec,gram,nsg)

         call cfill(qt0,1./nsg,nsg)
         do i=1,ns
            s=-(1./nsg)*vlsum(qt0(1+i*nsg),nsg)
            call cadd(qt0(1+i*nsg),s,nsg)
         enddo

         call dump_serial(gram,nsg*nsg,'ops/gt0_ ',nid)
      endif

      call rzero_ops

      inel=1
      ieg1=0

      do while (ieg1+1.le.nelgv)
         ieg0=ieg1+1
         ieg1=min(ieg1+inel+nelp-1,nelgv)
         nel=ieg1-ieg0+1
         n=lxyz*(ieg1-ieg0+1)
         call rsnapsm(uu,tt,ieg0,ieg1)

         call setgeom(gfac,w9,ieg0,ieg1,lxyz,ng,nid)
         call setvisc(visc,w,ieg0,ieg1,lxyz,nid)
         call setmass(mass,wv1,ieg0,ieg1,lxyz)
         call setrxp(rxp,rxpt,ieg0,ieg1)

         m=n*ndim
         call setzz(zz,uu,qu,wvf1,wvf2,ilgls,iglls,m,ms(nid+1),nsg)
         call setbb(bb,zz,mass,wvf1,wvf2,wvf12,ilgls(1),ms,n,ndim,igsh)
         call setaa(aa,zz,visc,gfac,wvf1,wvf2,wvf12,ilgls(1),
     $      ms,n,nel,ndim,ng,igsh)
         call setcc(cc,zz,zz,rxp,wvf1,wvf2,wvf3,wvf4,ilgls(1),
     $      ms,msr,n,ndim,ndim,nel,igsh)
         
         if (iftherm) then
         call dump_serial(qt,nsg*nsg,'ops/qt ',nid)
         call setzz(zt,tt,qt,wvf1,wvf2,ilgls,iglls,n,ms(nid+1),nsg)
         call setbb(bbt,zt,mass,wvf1,wvf2,wvf12,ilgls(1),ms,n,1,igsh)

         call setaa(aat,zt,visc,gfac,wvf1,wvf2,wvf12,ilgls(1),
     $      ms,n,nel,1,ng,igsh)
         call setcc(cct,zz,zt,rxp,wvf1,wvf2,wvf3,wvf4,ilgls(1),
     $      ms,msr,n,ndim,1,nel,igsh)
         endif

         call setzz(zz,uu,qu0,wvf1,wvf2,ilgls,iglls,m,ms(nid+1),nsg)
         call setbb(bb0,zz,mass,wvf1,wvf2,wvf12,ilgls(1),ms,n,ndim,igsh)
         call setaa(aa0,zz,visc,gfac,wvf1,wvf2,wvf12,ilgls(1),
     $      ms,n,nel,ndim,ng,igsh)
         call setcc(cc0,zz,zz,rxp,wvf1,wvf2,wvf3,wvf4,ilgls(1),
     $      ms,msr,n,ndim,ndim,nel,igsh)

         if (iftherm) then
            call setzz(zt,tt,qt0,wvf1,wvf2,ilgls,iglls,m,ms(nid+1),nsg)
            call setbb(bbt0,zt,mass,
     $         wvf1,wvf2,wvf12,ilgls(1),ms,n,1,igsh)
            call setaa(aat0,zt,visc,gfac,wvf1,wvf2,wvf12,ilgls(1),
     $         ms,n,nel,1,ng,igsh)
            call setcc(cct0,zz,zt,rxp,wvf1,wvf2,wvf3,wvf4,ilgls(1),
     $         ms,msr,n,ndim,1,nel,igsh)
         endif
      enddo

      nl=ms(nid+1)*ns*ns
      call setcc_transfer(cc,nl)

      call dump_parallel(bb,ms(nid+1)*ns,'ops/bbu_z ',nid)
      call dump_parallel(aa,ms(nid+1)*ns,'ops/aau_z ',nid)
      call dump_parallel(cc,nl,'ops/ccu_z ',nid)

      call setcc_transfer(cc0,nl)

      call dump_parallel(bb0,ms(nid+1)*ns,'ops/bbu0_z ',nid)
      call dump_parallel(aa0,ms(nid+1)*ns,'ops/aau0_z ',nid)
      call dump_parallel(cc0,nl,'ops/ccu0_z ',nid)

      if (iftherm) then
         call setcc_transfer(cct,nl)

         call dump_parallel(bbt,ms(nid+1)*ns,'ops/bbt_z ',nid)
         call dump_parallel(aat,ms(nid+1)*ns,'ops/aat_z ',nid)
         call dump_parallel(cct,nl,'ops/cct_z ',nid)

         call setcc_transfer(cct0,nl)

         call dump_parallel(bbt0,ms(nid+1)*ns,'ops/bbt0_z ',nid)
         call dump_parallel(aat0,ms(nid+1)*ns,'ops/aat0_z ',nid)
         call dump_parallel(cct0,nl,'ops/cct0_z ',nid)
      endif

      call exitt0
    1 format(i8,i8,1p2e15.6,' zcomp')

      return
      end
c-----------------------------------------------------------------------
      subroutine gsub0(gram,gram0,wk,ns,nsg)

      real gram(1),gram0(1),wk(1)

      call rzero(gram0,nsg)

      do i=1,ns
         call add2(gram0,gram(1+(i-1)*nsg),nsg)
      enddo

      call gop(gram0,wk,'+  ',nsg)
      s=1./real(nsg)

      call cmult(gram0,s,nsg)
      gram00=s*vlsum(gram0,nsg)

      do i=1,ns
         call sub2(gram(1+(i-1)*nsg),gram0,nsg)
      enddo

      call cadd(gram0,-gram00,nsg)

      do i=1,ns
         call cadd(gram(1+(i-1)*nsg),-gram0(i),nsg)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setzz(z,u,evec,w1,w2,ilgls,iglls,n,ns,nsg)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer ilgls(1),iglls(1)

      real z(n,1),u(n,1)
      real evec(nsg,nsg)
      real w1(n),w2(n)

      do isg=1,nsg
         is=iglls(isg)
         i=ilgls(1)
         call mxm(u,n,evec(i,isg),ns,w1,1)
         call gop(w1,w2,'+  ',n)
         if (is.gt.0.and.is.le.ns) call copy(z(1,is),w1,n)
      enddo

c     call mxm(u,n,evec,nsg,z,nsg)

      return
      end
c-----------------------------------------------------------------------
      subroutine ilgls_setup(ilgls,msr,ms,msg,np,nid)

      integer ilgls(1),ms(1),msr(2,1)

      n=0
      msmin=msg/np

      l=1
      do id=0,np-1
         ms(id+1)=msmin+max(min(msg-msmin*np-id,1),0)
         if (id.le.nid-1) n=n+ms(id+1)
         msr(1,id+1)=l
         msr(2,id+1)=l+ms(id+1)-1
         l=l+ms(id+1)
      enddo

      do is=1,msmin+max(min(msg-msmin*np-nid,1),0)
         ilgls(is)=is+n
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine iglls_setup(iglls,itmp,ms,nsg,np,nid)

      integer iglls(nsg)
      integer ms(1)

      call izero(iglls,nsg)

      ic=1
      do id=0,np-1
         jc=1
         do is=1,ms(id+1)
            if (id.eq.nid) then
               iglls(ic)=jc
            endif
            ic=ic+1
            jc=jc+1
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setaa(a,z,visc,gfac,w1,w2,w3,igs,ns,n,nel,ndim,ng,igsh)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer igsh(2)
      real a(1),z(n,ndim,1),w1(n,ndim,1),w2(n,ndim,1)
      real w3(n,ndim,1,2)
      real visc(n),gfac(n,ng)

      integer ns(1)

      nsg=ivlsum(ns,mp)

      call copy(w1,z,n*ndim*ns(mid+1))
      call copy(w2,z,n*ndim*ns(mid+1))

      imesh=1
      isd=1

      do i=1,ns(mid+1)*ndim
         call aop(w1(1,i,1),z(1,i,1),visc,gfac,imesh,isd,nel)
      enddo

      j=igs

      nsmax=ivlmax(ns,mp)

      do id=0,mp-1
         do k=1,ns(mod(mid+id,mp)+1)
            do i=1,ns(mid+1)
               a(j+(i-1)*nsg)=a(j+(i-1)*nsg)+
     $            vlsc2(w1(1,1,i),w2(1,1,k),n*ndim)
            enddo
            j=mod(j,nsg)+1
         enddo
         call shift(igsh,w2,w3,n*ndim*nsmax)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setbb(b,u,mass,w1,w2,w3,igs,ns,n,ndim,igsh)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer ns(1),igsh(2)

      real b(1),u(n,ndim,1),mass(n),
     $     w1(n,ndim,1),w2(n,ndim,1)

      real w3(n,ndim,1,3)

      nsg=ivlsum(ns,mp)

      call copy(w1,u,n*ndim*ns(mid+1))
      call copy(w2,u,n*ndim*ns(mid+1))

      do i=1,ns(mid+1)*ndim
         call col2(w1(1,i,1),mass,n)
      enddo

      j=igs

      nsmax=ivlmax(ns,mp)

      do id=0,mp-1
         do k=1,ns(mod(mid+id,mp)+1)
            do i=1,ns(mid+1)
               b(j+(i-1)*nsg)=b(j+(i-1)*nsg)+
     $            vlsc2(w1(1,1,i),w2(1,1,k),n*ndim)
            enddo
            j=mod(j,nsg)+1
         enddo
         call shift(igsh,w2,w3,n*ndim*nsmax)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc_snap(cc)

      include 'SIZE'
      include 'MOR'
      include 'TSTEP'

      parameter (lt=lx1*ly1*lz1*lelt)

      common /scrcc/ t1(lt),t2(lt),t3(lt)

      real cc(ns,ns,ns)

      n=lx1*ly1*lz1*nelv

      call rone(ones,lx1*ly1*lz1*nelv)
      ifield=1

      do ks=1,ns
         do js=1,ns
            call convect_new(t1,us0(1,1,js),.false.,
     $         us0(1,1,ks),us0(1,2,ks),us0(1,mdim,ks),.false.)
            call convect_new(t2,us0(1,2,js),.false.,
     $         us0(1,1,ks),us0(1,2,ks),us0(1,mdim,ks),.false.)
            do is=1,ns
               cc(is,js,ks)=cc(is,js,ks)
     $            +vlsc2(t1,us0(1,1,is),n)+vlsc2(t2,us0(1,2,is),n)
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc_lgc(c,z,t,rxp,w1,w2,w3,w4,ns,nsg,n,mdim,ndim,nel)

      real rxp(1)
      real c(ns,nsg,nsg),z(n,ndim,ns),t(n,mdim,ns),
     $     w1(n,mdim,ns),w2(n,ndim,ns),w3(n,mdim,ns),w4(n,mdim,ns)

      do ks=1,ns
         do js=1,ns
c           call convect_new(w3,t(1,1,js),.false.,
c    $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.)
c           call convect_new(w3(1,2,1),t(1,2,js),.false.,
c    $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.)
            call conv(w3,t(1,1,js),.false.,
     $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.,rxp,nel)
            call conv(w3(1,2,1),t(1,2,js),.false.,
     $         z(1,1,ks),z(1,2,ks),z(1,mdim,ks),.false.,rxp,nel)
            do is=1,ns
               c(is,js,ks)=c(is,js,ks)
     $            +vlsc2(w3,t(1,1,is),n*mdim)
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc(
     $   c,z,t,rxp,w1,w2,w3,w4,igs,ns,nsr,n,ndim,mdim,nel,igsh)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer ns(1),nsr(2,1),igsh(2)

      real rxp(1)
      real c(1),z(n,ndim,1),t(n,mdim,1),
     $     w1(n,mdim,1),w2(n,mdim,1),w3(n,ndim,1),w4(n,mdim,1)

      nsg=ivlsum(ns,mp)
      nsmax=ivlmax(ns,mp)

      call copy(w1,t,n*mdim*ns(mid+1))
      call copy(w2,t,n*mdim*ns(mid+1))
      call copy(w3,z,n*ndim*ns(mid+1))

      ms=ns(mid+1)

      j=igs
      k=igs

      do kid=0,mp-1
         k1=nsr(1,mod(kid+mid,mp)+1)
         k2=nsr(2,mod(kid+mid,mp)+1)
         do jid=0,mp-1
            j1=nsr(1,mod(jid+mid,mp)+1)
            j2=nsr(2,mod(jid+mid,mp)+1)
            k=1
            do ks=k1,k2
               j=1
               do js=j1,j2
                  call conv(w4,w2(1,1,j),.false.,
     $               w3(1,1,k),w3(1,2,k),w3(1,ndim,k),.false.,rxp,nel)
                  call conv(w4(1,2,1),w2(1,2,j),.false.,
     $               w3(1,1,k),w3(1,2,k),w3(1,ndim,k),.false.,rxp,nel)
                  do i=1,ms
                    c(i+(js-1)*ms+(ks-1)*ms*nsg)=
     $              c(i+(js-1)*ms+(ks-1)*ms*nsg)
     $                  +vlsc2(w4(1,1,1),w1(1,1,i),n*mdim)
                  enddo
                  j=j+1
               enddo
               k=k+1
            enddo
            call shift(igsh,w2,w4,n*mdim*nsmax)
         enddo
         call shift(igsh,w3,w4,n*ndim*nsmax)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc_transfer(c,nl)

      include 'POST'
      include 'PARALLEL'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real c(1)

      nsg=ns

      nl=0

      n=ms(mid+1)
      mmax=2*lsg*lsg
      ntot=nsg*nsg*nsg

      nmin=ntot/mp

      id=0
      msum=nmin+max(min(ntot-nmin*mp-id,1),0)
      mm=msum

      id=0
      l=1

      isg0=ilgls(1)
      isg1=ilgls(ms(mid+1))

      write (6,*) isg0,isg1,'isg'

      m=1
      do k=1,nsg
      do j=1,nsg
      do i=1,nsg
         if (mm.lt.l) then
            id=id+1
            write (6,*) 'mm=',mm,l
            mm=mm+nmin+max(min(ntot-nmin*mp-id,1),0)
         endif
         if (i.ge.isg0.and.i.le.isg1) then
            iw1(m)=id
            iw2(m)=l
            m=m+1
         endif
         l=l+1
      enddo
      enddo
      enddo

      nl=n*nsg*nsg

      do id=0,mp-1
         if (mid.eq.id) then
            do i=1,nl
               write (6,1) id,i,nl,iw1(i),iw2(i),c(i)
            enddo
         endif
         call nekgsync
      enddo

      call fgslib_crystal_tuple_transfer(cr_h,nl,mmax,iw1,1,iw2,1,c,1,1)
      call fgslib_crystal_tuple_sort(cr_h,nl,iw1,1,iw2,1,c,1,2,1)

      do id=0,mp-1
         if (mid.eq.id) then
            do i=1,nl
               write (6,2) id,i,nl,iw1(i),iw2(i),c(i)
            enddo
         endif
         call nekgsync
      enddo

    1 format(i5,i8,i8,i5,i8,1pe13.4,' crystal1')
    2 format(i5,i8,i8,i5,i8,1pe13.4,' crystal2')

      return
      end
c-----------------------------------------------------------------------
      subroutine setgeom(gfac,w,ieg0,ieg1,lxyz,ng)!,nid)

      include 'SIZE' ! nid
      include 'GEOM'
      include 'PARALLEL'

      real gfac(lxyz,ieg1-ieg0+1,ng),w(lxyz*(ieg1-ieg0+1)*ng)

      call rzero(gfac,lxyz*(ieg1-ieg0+1)*max(ng,4))

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
            call copy(gfac(1,ieg-ieg0+1,1),g1m1(1,1,1,ie),lxyz)
            call copy(gfac(1,ieg-ieg0+1,2),g2m1(1,1,1,ie),lxyz)
            call copy(gfac(1,ieg-ieg0+1,4),g4m1(1,1,1,ie),lxyz)
            if (ng.eq.6) then
               call copy(gfac(1,ieg-ieg0+1,3),g3m1(1,1,1,ie),lxyz)
               call copy(gfac(1,ieg-ieg0+1,5),g5m1(1,1,1,ie),lxyz)
               call copy(gfac(1,ieg-ieg0+1,6),g6m1(1,1,1,ie),lxyz)
            endif
         endif
      enddo

      call gop(gfac,w,'+  ',lxyz*(ieg1-ieg0+1)*max(ng,4))

      return
      end
c-----------------------------------------------------------------------
      subroutine setvisc(visc,w,ieg0,ieg1,lxyz)!,nid)

      include 'SIZE'
      include 'PARALLEL'

      real visc(lxyz,ieg1-ieg0+1), w(lxyz*(ieg1-ieg0+1))

      call rone(visc,lxyz*(ieg1-ieg0+1))
      return

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
c           call copy(visc(1,i-ieg0+1),h1(1,1,1,ie),lxyz)
         endif
      enddo

      call gop(visc,w,'+  ',lxyz*(ieg1-ieg0+1))

      return
      end
c-----------------------------------------------------------------------
      subroutine setmass(mass,w,ieg0,ieg1,lxyz)!,nid)

      include 'SIZE' ! nid
      include 'MASS'
      include 'PARALLEL'

      real mass(lxyz,ieg1-ieg0+1), w(lxyz*(ieg1-ieg0+1))

      ! TODO: set bm1

      call rzero(mass,lxyz*(ieg1-ieg0+1))

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
            call copy(mass(1,ieg-ieg0+1),bm1(1,1,1,ie),lxyz)
         endif
      enddo

      call gop(mass,w,'+  ',lxyz*(ieg1-ieg0+1))

      return
      end
c-----------------------------------------------------------------------
      subroutine setconv(rxd,w,ieg0,ieg1,lxyz)!,ndim,nid)

      include 'SIZE' ! ndim & nid
      include 'GEOM'
      include 'PARALLEL'

      real rxd(lxyz,ndim*ndim,ieg1-ieg0+1), w1(lxyz,9,(ieg1-ieg0+1)),
     $     w(lxyz*(ieg1-ieg0+1))

c     logical if3d

      call rzero(w,lxyz*(ieg1-ieg0+1))
c     if3d=ndim.eq.3

      do ieg=ieg0,ieg1
         if (gllnid(ieg).eq.nid) then
            ie=gllel(ieg)
            if (ndim.eq.2) then
               call copy(w1(1,1,ieg-ieg0+1),rxm1(1,1,1,ie),lxyz)
               call copy(w1(1,2,ieg-ieg0+1),rym1(1,1,1,ie),lxyz)
               call copy(w1(1,3,ieg-ieg0+1),sxm1(1,1,1,ie),lxyz)
               call copy(w1(1,4,ieg-ieg0+1),sym1(1,1,1,ie),lxyz)
            else
               call copy(w1(1,1,ieg-ieg0+1),rxm1(1,1,1,ie),lxyz)
               call copy(w1(1,2,ieg-ieg0+1),rym1(1,1,1,ie),lxyz)
               call copy(w1(1,3,ieg-ieg0+1),rzm1(1,1,1,ie),lxyz)
               call copy(w1(1,4,ieg-ieg0+1),sxm1(1,1,1,ie),lxyz)
               call copy(w1(1,5,ieg-ieg0+1),sym1(1,1,1,ie),lxyz)
               call copy(w1(1,6,ieg-ieg0+1),rzm1(1,1,1,ie),lxyz)
               call copy(w1(1,7,ieg-ieg0+1),txm1(1,1,1,ie),lxyz)
               call copy(w1(1,8,ieg-ieg0+1),tym1(1,1,1,ie),lxyz)
               call copy(w1(1,9,ieg-ieg0+1),tzm1(1,1,1,ie),lxyz)
            endif
         endif
      enddo

      call gop(w1,w2,'+  ',lxyz*(ieg1-ieg0+1)*ndim*ndim)

      do ie=1,ieg1-ieg0+1
         call intp_rstd(rxd(1,1,ie),w1(1,1,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         call intp_rstd(rxd(1,2,ie),w1(1,2,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         call intp_rstd(rxd(1,3,ie),w1(1,3,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         call intp_rstd(rxd(1,4,ie),w1(1,4,ie),lx1,lxd,if3d,0) ! 0 --> fwd

         if (ndim.eq.3) then
            call intp_rstd(rxd(1,5,ie),w1(1,5,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,6,ie),w1(1,6,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,7,ie),w1(1,7,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,8,ie),w1(1,8,ie),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rxd(1,9,ie),w1(1,9,ie),lx1,lxd,if3d,0) ! 0 --> fwd
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
c     subroutine ax(au,u,helm1,helm2,imesh,isd)
      subroutine aop(au,u,visc,gfac,imesh,isd,mel)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
c     include 'PARALLEL'
c     include 'CTIMER'

      common /fastax/ wddx(lx1,lx1),wddyt(ly1,ly1),wddzt(lz1,lz1)
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv

      real           au    (lx1,ly1,lz1,1)
     $ ,             u     (lx1,ly1,lz1,1)
     $ ,             visc (lx1,ly1,lz1,1)
     $ ,             helm2 (lx1,ly1,lz1,1)
      common /ctmp1/ dudr  (lx1,ly1,lz1)
     $ ,             duds  (lx1,ly1,lz1)
     $ ,             dudt  (lx1,ly1,lz1)
     $ ,             tmp1  (lx1,ly1,lz1)
     $ ,             tmp2  (lx1,ly1,lz1)
     $ ,             tmp3  (lx1,ly1,lz1)

      real           tm1   (lx1,ly1,lz1)
      real           tm2   (lx1,ly1,lz1)
      real           tm3   (lx1,ly1,lz1)
      real           duax  (lx1)
      real           ysm1  (lx1)
      equivalence    (dudr,tm1),(duds,tm2),(dudt,tm3)

      real gfac(lx1,ly1,lz1,mel,ng)

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nxy=lx1*ly1
      nyz=ly1*lz1
      nxz=lx1*lz1
      nxyz=lx1*ly1*lz1
      ntot=nxyz*mel

c     if (.not.ifsolv) call setfast(visc,helm2,imesh)
      call rzero (au,ntot)

      do ie=1,mel
 
        if (ifaxis) call setaxdy ( ifrzer(ie) )

        if (ldim.eq.2) then

c       2-d case ...............

c          if (iffast(ie)) then
           if (.false.) then

c          Fast 2-d mode: constant properties and undeformed element

           h1 = visc(1,1,1,ie)
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           call mxm   (u(1,1,1,ie),lx1,wddyt,ly1,tm2,ly1)
           call col2  (tm1,gfac(1,1,1,ie,4),nxyz)
           call col2  (tm2,gfac(1,1,1,ie,5),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)

           else


           call mxm  (dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           call mxm  (u(1,1,1,ie),lx1,dytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,gfac(1,1,1,ie,1),nxyz)
           call col3 (tmp2,duds,gfac(1,1,1,ie,2),nxyz)
           call addcol3 (tmp1,duds,gfac(1,1,1,ie,4),nxyz)
           call addcol3 (tmp2,dudr,gfac(1,1,1,ie,4),nxyz)
           call col2 (tmp1,visc(1,1,1,ie),nxyz)
           call col2 (tmp2,visc(1,1,1,ie),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,ie),tm1,nxyz)
           call add2 (au(1,1,1,ie),tm2,nxyz)
        endif

        else

c       3-d case ...............

           if (iffast(ie)) then

c          Fast 3-d mode: constant properties and undeformed element

           h1 = visc(1,1,1,ie)
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           do 5 iz=1,lz1
           call mxm   (u(1,1,iz,ie),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
 5         continue
           call mxm   (u(1,1,1,ie),nxy,wddzt,lz1,tm3,lz1)
           call col2  (tm1,gfac(1,1,1,ie,4),nxyz)
           call col2  (tm2,gfac(1,1,1,ie,5),nxyz)
           call col2  (tm3,gfac(1,1,1,ie,6),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call add2  (au(1,1,1,ie),tm3,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)

           else

           call mxm(dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,ie),lx1,dytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,ie),nxy,dztm1,lz1,dudt,lz1)
           call col3    (tmp1,dudr,gfac(1,1,1,ie,1),nxyz)
           call col3    (tmp2,duds,gfac(1,1,1,ie,2),nxyz)
           call col3    (tmp3,dudt,gfac(1,1,1,ie,3),nxyz)
           if (ifdfrm(ie)) then
              call addcol3 (tmp1,duds,gfac(1,1,1,ie,4),nxyz)
              call addcol3 (tmp1,dudt,gfac(1,1,1,ie,5),nxyz)
              call addcol3 (tmp2,dudr,gfac(1,1,1,ie,4),nxyz)
              call addcol3 (tmp2,dudt,gfac(1,1,1,ie,6),nxyz)
              call addcol3 (tmp3,dudr,gfac(1,1,1,ie,5),nxyz)
              call addcol3 (tmp3,duds,gfac(1,1,1,ie,6),nxyz)
           endif
           call col2 (tmp1,visc(1,1,1,ie),nxyz)
           call col2 (tmp2,visc(1,1,1,ie),nxyz)
           call col2 (tmp3,visc(1,1,1,ie),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
           call add2 (au(1,1,1,ie),tm1,nxyz)
           call add2 (au(1,1,1,ie),tm2,nxyz)
           call add2 (au(1,1,1,ie),tm3,nxyz)

           endif

        endif

      enddo

c     If axisymmetric, add a diagonal term in the radial direction (ISD=2)

      if (ifaxis.and.(isd.eq.2)) then
         do ie=1,mel

            if (ifrzer(ie)) then
               call mxm(u  (1,1,1,ie),lx1,datm1,ly1,duax,1)
               call mxm(ym1(1,1,1,ie),lx1,datm1,ly1,ysm1,1)
            endif

            do 190 j=1,ly1
            do 190 i=1,lx1
                if (ym1(i,j,1,ie).ne.0.) then
                  if (ifrzer(ie)) then
                     term1 = 0.0
                     if(j.ne.1) 
     $             term1 = bm1(i,j,1,ie)*u(i,j,1,ie)/ym1(i,j,1,ie)**2
                     term2 =  wxm1(i)*wam1(1)*dam1(1,j)*duax(i)
     $                       *jacm1(i,1,1,ie)/ysm1(i)
                  else
                   term1 = bm1(i,j,1,ie)*u(i,j,1,ie)/ym1(i,j,1,ie)**2
                     term2 = 0.
                  endif
                  au(i,j,1,ie) = au(i,j,1,ie)
     $                          + visc(i,j,1,ie)*(term1+term2)
                endif
  190       continue
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
c     subroutine setk(uk,u,z,mass,w1,w2,ns,nsg,n,ndim)
c
c     real uk(ns,nsg),u(n,ndim,ns),z(n,ndim,ns),
c    $     mass(n),w1(n,ndim,ns),w2(n,ndim,ns)
c
c     call copy(w1,z,n*ndim*ns)
c     call copy(w2,u,n*ndim*ns)
c
c     do i=1,ns*ndim
c        call col2(w2(1,i,1),mass,n)
c     enddo
c
c     do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
c        do js=1,ns
c        do is=1,ns
c           b(is,js+ns*ioff)=b(is,js+ns*ioff)
c    $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
c        enddo
c        enddo
c        
c        call shift(w1,n*ndim*ns)
c     enddo
c
c     return
c     end
c-----------------------------------------------------------------------
      subroutine conv(bdu,u,ifuf,cx,cy,cz,ifcf,rxp,mel)

c     Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians
c
      include 'SIZE'
      include 'TOTAL'

      real bdu(1),u(1),cx(1),cy(1),cz(1)
      real rxp(lxd*lyd*lzd,ldim*ldim,mel)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)

      integer e

      nxyz1=lx1*ly1*lz1
      nxyzd=lxd*lyd*lzd

      nxyzu=nxyz1
      if (ifuf) nxyzu=nxyzd

      nxyzc = nxyz1
      if (ifcf) nxyzc=nxyzd

      iu=1 ! pointer to scalar field u
      ic=1 ! pointer to vector field C
      ib=1 ! pointer to scalar field Bdu

      do e=1,mel
         if (ifcf) then
            call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
            call copy(tr(1,2),cy(ic),nxyzd)
            if (if3d) call copy(tr(1,3),cz(ic),nxyzd)
         else  ! map coarse velocity to fine mesh (C-->F)
            call intp_rstd(fx,cx(ic),lx1,lxd,if3d,0)
            call intp_rstd(fy,cy(ic),lx1,lxd,if3d,0)
            if (if3d) call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0)

            if (if3d) then  ! Convert convector F to r-s-t coordinates
               do i=1,nxyzd
                  tr(i,1)=rxp(i,1,e)*fx(i)+
     $                    rxp(i,2,e)*fy(i)+
     $                    rxp(i,3,e)*fz(i)
                  tr(i,2)=rxp(i,4,e)*fx(i)+
     $                    rxp(i,5,e)*fy(i)+
     $                    rxp(i,6,e)*fz(i)
                  tr(i,3)=rxp(i,7,e)*fx(i)+
     $                    rxp(i,8,e)*fy(i)+
     $                    rxp(i,9,e)*fz(i)
               enddo
            else
               do i=1,nxyzd
                  tr(i,1)=rxp(i,1,e)*fx(i)+rxp(i,2,e)*fy(i)
                  tr(i,2)=rxp(i,3,e)*fx(i)+rxp(i,4,e)*fy(i)
               enddo
           endif
         endif

         if (ifuf) then
            call grad_rst(ur,us,ut,u(iu),lxd,if3d)
         else
            call intp_rstd(uf,u(iu),lx1,lxd,if3d,0)
            call grad_rst(ur,us,ut,uf,lxd,if3d)
         endif

         if (if3d) then
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i)=tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
            enddo
         else
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i)=tr(i,1)*ur(i)+tr(i,2)*us(i)
            enddo
         endif

         call intp_rstd(bdu(ib),uf,lx1,lxd,if3d,1)

         ic=ic+nxyzc
         iu=iu+nxyzu
         ib=ib+nxyz1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setrxp(rxp,tmp,ieg0,ieg1)
c
c     Eulerian scheme, add convection term to forcing function
c     at current time step.
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP' ! for istep
      include 'PARALLEL'

      common /dealias1/ zd(lxd),wd(lxd)
      real rxp(lxd*lyd*lzd,ldim*ldim,ieg1-ieg0+1)
      real tmp(lxd*lyd*lzd,ldim*ldim,ieg1-ieg0+1)
      integer e

      integer ilstep
      save    ilstep
      data    ilstep /-1/

c     if (.not.ifgeom.and.ilstep.gt.1) return  ! already computed
c     if (ifgeom.and.ilstep.eq.istep)  return  ! already computed
      ilstep=istep

      nxyz1=lx1*ly1*lz1
      nxyzd=lxd*lyd*lzd

      call zwgl(zd,wd,lxd)  ! zwgl -- NOT zwgll!
      call rzero(rxp,lxd*lyd*lzd*ldim*ldim*(ieg1-ieg0+1))

      if (if3d) then
         do ieg=ieg0,ieg1
            ie=ieg-ieg0+1
            if (gllnid(ieg).eq.nid) then
               e=gllel(ieg)

c              Interpolate z+ and z- into fine mesh, translate to r-s-t coords

               call intp_rstd(rxp(1,1,ie),rxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,2,ie),rym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,3,ie),rzm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,4,ie),sxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,5,ie),sym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,6,ie),szm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,7,ie),txm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,8,ie),tym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,9,ie),tzm1(1,1,1,e),lx1,lxd,if3d,0)

               l=0
               do k=1,lzd
               do j=1,lyd
               do i=1,lxd
                  l=l+1
                  w=wd(i)*wd(j)*wd(k)
                  do ii=1,9
                     rxp(l,ii,ie)=w*rxp(l,ii,ie)
                  enddo
               enddo
               enddo
               enddo
            endif
         enddo

      else ! 2D
         do ieg=ieg0,ieg1
            ie=ieg-ieg0+1
            if (gllnid(ieg).eq.nid) then
               e=gllel(ieg)

c              Interpolate z+ and z- into fine mesh, translate to r-s-t coords

               call intp_rstd(rxp(1,1,ie),rxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,2,ie),rym1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,3,ie),sxm1(1,1,1,e),lx1,lxd,if3d,0)
               call intp_rstd(rxp(1,4,ie),sym1(1,1,1,e),lx1,lxd,if3d,0)

               l=0
               do j=1,lyd
               do i=1,lxd
                  l=l+1
                  w=wd(i)*wd(j)
                  do ii=1,4
                     rxp(l,ii,ie)=w*rxp(l,ii,ie)
                  enddo
               enddo
               enddo
            endif
         enddo
      endif

      call gop(rxp,tmp,'+  ',lxd*lyd*lzd*ldim*ldim*(ieg1-ieg0+1))

      return
      end
c-----------------------------------------------------------------------
      subroutine rzero_ops

      include 'POST'

      call rzero(gua,ms(nid+1)*ns)
      call rzero(gub,ms(nid+1)*ns)
      call rzero(guc,ms(nid+1)*ns*ns)
      call rzero(guc2,ms(nid+1)*ns*ns)

      if (iftherm) then
         call rzero(gta,ms(nid+1)*ns)
         call rzero(gtb,ms(nid+1)*ns)
         call rzero(gtc,ms(nid+1)*ns*ns)
         call rzero(gtc2,ms(nid+1)*ns*ns)
      endif

      call rzero(aa,ms(nid+1)*ns)
      call rzero(bb,ms(nid+1)*ns)
      call rzero(cc,ms(nid+1)*ns*ns)

      call rzero(aa0,ms(nid+1)*ns)
      call rzero(bb0,ms(nid+1)*ns)
      call rzero(cc0,ms(nid+1)*ns*ns)

      call rzero(aat,ms(nid+1)*ns)
      call rzero(bbt,ms(nid+1)*ns)
      call rzero(cct,ms(nid+1)*ns*ns)

      call rzero(aat0,ms(nid+1)*ns)
      call rzero(bbt0,ms(nid+1)*ns)
      call rzero(cct0,ms(nid+1)*ns*ns)

      return
      end
c-----------------------------------------------------------------------
      subroutine setgg(gram,gub,eevec,ilgls,nsg,ns)

      real gram(nsg,nsg),gub(nsg,nsg),eevec(nsg*nsg)
      integer ilgls(nsg)

      call rzero(gram,nsg*nsg)

      do is=1,ns
         isg=ilgls(is)
         call copy(gram(1,isg),gub(1,is),nsg)
      enddo

      call gop(gram,eevec,'+  ',nsg*nsg)

      return
      end
c-----------------------------------------------------------------------
      subroutine setqq(eevec,evecp,gram,nsg)

      parameter (l=3)

      common /test/ t1(l,l),t2(l,l)

      real eevec(nsg,nsg),evecp(nsg),gram(nsg,nsg)

      call genevec(eevec,evecp,gram,1)

      call dump_serial(eevec,nsg*nsg,'ops/eevec1 ',0)

      call mxm(gram,nsg,eevec,nsg,t1,nsg)

      do i=1,nsg
         s=1./sqrt(evecp(i))
         call cmult(eevec(1,i),s,nsg)
      enddo

      return
      end
c-----------------------------------------------------------------------
