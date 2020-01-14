c-----------------------------------------------------------------------
      subroutine drivep

      include 'POST'
      include 'MASS'
      include 'PARALLEL'

      character*127 flist

      nelp=512
      nelp=256
      nelp=4
      nsnap=ns

      ! assigned snapshot range

      isg0=1
      isg1=nsg

      call reade_init(flist,ns)

      inel=1
      ieg1=0

      ns=ls
      nsg=lsg

      call rzero(gram,ns*nsg)
      call rzero(aa,ns*nsg)
      call rzero(bb,ns*nsg)
      call rzero(cc,ns*nsg*nsg)

      do while (ieg1+1.le.nelgv)
         ieg0=ieg1+1
         ieg1=min(ieg1+inel+nelp-1,nelgv)
         nel=ieg1-ieg0+1
         n=lxyz*(ieg1-ieg0+1)
         call reade_dummy(uu,ieg0,ieg1)
         call setmass(mass,wv1,ieg0,ieg1,lxyz)
         call setgg(gram,uu,mass,wvf1,wvf2,ns,nsg,n,ndim)
      enddo

      call dump_serial(gram,ns*ns,'ops/gramp ',nid)

      ! eigendecomposition here or external process

      write (6,*) 'ns=',ns
      mmm=ns*ns
      call read_serial(evecp,mmm,'ops/evecp ',ug,nid)

      inel=1
      ieg1=0

      call rone(visc,n)

      do while (ieg1+1.le.nelgv)
         ieg0=ie1+1
         ieg1=min(ieg1+inel+nelp-1,nelgv)
         nel=ieg1-ieg0+1
         n=lxyz*(ieg1-ieg0+1)

c        call reade(uu,ieg0,ieg1,'U',flist)
         call reade_dummy(uu,ieg0,ieg1)

c        call setvisc(visc,w,ieg0,ieg1,lxyz,nid) ! not required
         call setgeom(gfac,w9,ieg0,ieg1,lxyz,3*(ndim-1),nid)
c        call setconv(rxd) ! TODO: implement

         call setmass(mass,wv1,ieg0,ieg1,lxyz,nid)

         call setzz(zz,uu,evecp,wvf1,wvf2,n*ndim,ns,nsg)

c        do k=1,ns
c        do i=1,lxyz*nel
c           write (6,1) i,k,ub(i,k),zz(i+(k-1)*lxyz*nel*ndim),
c    $                      ub(i,k)/zz(i+(k-1)*lxyz*nel*ndim)
c        enddo
c        enddo

         call setbb(bb,zz,mass,wvf1,wvf2,ns,nsg,n,ndim)

c        call setaa(aa,zz,visc,gfac,wvf1,wvf2,ns,nsg,n,ndim)
c        call setcc(c,z,wvf1,wvf2,ns,nsg,n,ndim,ndim) ! TODO: implement

c        call setkk(uk,zz,mass,wvf1,wvf2,ns,nsg,n,ndim)
      enddo

      call rzero(zsc,nsg)

      do i=1,ns
         zsc(i)=sqrt(bb(i+(i-1)*ns))
      enddo

      call invcol1(zsc,ns)

      call gop(zsc,ug,'+  ',nsg)

      do j=1,nsg
      do i=1,ns
         bb(i+(j-1)*ns)=bb(i+(j-1)*ns)*zsc(i)*zsc(j)
      enddo
      enddo

      do j=1,10
      do i=1,10
         write (6,*) bb(i+(j-1)*10),'bb'
      enddo
      enddo

      call dump_serial(bb,mmm,'ops/bup ',nid)
      call exitt0

    1 format (i8,i8,1p3e16.5,' zz')

      return
      end
c-----------------------------------------------------------------------
      subroutine reade_init(flist,ns)

      character*127 flist

      return
      end
c-----------------------------------------------------------------------
      subroutine reade_dummy(v,ieg0,ieg1)

      include 'POST'

      parameter (ll=lx1*ly1*lz1)

      real v(lxyz,ieg1-ieg0+1,ldim,ns)

      do is=1,ns
      do idim=1,ndim
      do ie=ieg0,ieg1
         call copy(v(1,ie-ieg0+1,idim,is),us0(1+(ie-1)*ll,idim,is),ll)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine reade(u,ieg0,ieg1,clist,flist)

      real u(1)
      logical ifread(0:4)
      character*127 clist
      character*127 fname

      call parse_clist(ifread,clist)
c     call reade_helper(u,ieg0,ieg1,ifread,fname)

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_clist(ifread,clist)

      logical ifread(0:4)
      character*128 clist

      ! no support for thermal for now

      do i=0,4
         ifread(i)=.false.
      enddo

c     if (clist matches 'U') then
         ifread(1)=.true.
c     else
c        call exitti('unsupported clist$',1)
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setzz(z,u,evec,w1,w2,n,ns,nsg)

c     include 'SIZE'
c     include 'PARALLEL'

      real z(n,ns),u(n,ns)
      real evec(ns,nsg)
      real w1(n),w2(n)

c     do isg=1,nsg
c        is=iglls(isg)
c        is=isg
c        call mxm(u,n,evec(1,isg),ns,w1,1)
c        call gop(w1,w2,'+  ',n)
c        if (is.gt.0.and.is.le.ns) call copy(z(1,is),w1,n)

c     enddo

      do i=1,ns
         call mxm(u,n,evec(1,i),ns,z(1,i),1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      function ilgls(is)

      ilgls=0

      return
      end
c-----------------------------------------------------------------------
      function iglls(isg)

      iglls=0

      return
      end
c-----------------------------------------------------------------------
      subroutine setaa(a,z,w1,w2,ns,nsg,n,ndim)

      real a(ns,nsg),z(n,ndim,ns),w1(n,ndim,ns),w2(n,ndim,ns)

      call copy(w1,z,n*ndim*ns)
      call copy(w2,z,n*ndim*ns)

      do i=1,ns*ndim
         call ax(w2(1,i,1))
      enddo

      do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
         do js=1,ns
         do is=1,ns
            a(is,js+ns*ioff)=a(is,js+ns*ioff)
     $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
         enddo
         enddo
         
         call shift(w1,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setbb(b,u,mass,w1,w2,ns,nsg,n,ndim)

      real b(ns,nsg),u(n,ndim,ns),mass(n),w1(n,ndim,ns),w2(n,ndim,ns)

      write (6,*) 'ns,nsg,n,ndim',ns,nsg,n,ndim

      call copy(w1,u,n*ndim*ns)
      call copy(w2,u,n*ndim*ns)

      do i=1,ns*ndim
         call col2(w2(1,i,1),mass,n)
      enddo

      do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
         do js=1,ns
         do is=1,ns
            b(is,js+ns*ioff)=b(is,js+ns*ioff)
     $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
         enddo
         enddo
         
         call shift(w1,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setcc(c,z,t,w1,w2,w3,w4,ns,nsg,n,mdim,ndim)

      real c(ns,nsg,nsg),z(n,mdim,ns),t(n,ndim,ns),
     $     w1(n,mdim,ns),w2(n,ndim,ns),w3(n,mdim,ns),w4(n,mdim,ns)

      call copy(w1,t,m*mdim*ns)
      call copy(w2,z,n*ndim*ns)
      call copy(w3,t,n*mdim*ns)

      do i=1,ns*ndim
c        call col2(w2(1,i,1),mass,n)
         ! apply grad on w1
      enddo

      do jof=0,nsg/ns-1 ! assume ns is the same across all processors
         do iof=0,nsg/ns-1
            do ks=1,ns
            do js=1,ns
c              call conv(w4,w1,w2) ! w4= mass * (w1 * w4)
               do is=1,ns
                  c(is,js+ns*iof,ks+ns*jof)=c(is,js+ns*iof,ks+ns*jof)
     $               +vlsc2(w3(1,1,is),w4(1,1,js),n*ndim)
               enddo
            enddo
            enddo
            
            call shift(w3,n*ndim*ns)
         enddo
         call shift(w2,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setgeom(gfac,w,ieg0,ieg1,lxyz,ng)!,nid)

      include 'SIZE' ! nid
      include 'GEOM'
      include 'PARALLEL'

      real gfac(lxyz,ieg1-ieg0+1,ng),w(lxyz*(ieg1-ieg0+1)*ng)

      call rzero(gfac,lxyz*(ieg1-ieg0+1)*ng)

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

      call gop(gfac,w,'+  ',lxyz*(ieg1-ieg0+1)*ng)

      return
      end
c-----------------------------------------------------------------------
      subroutine setvisc(visc,w,ieg0,ieg1,lxyz)!,nid)

      include 'SIZE'
      include 'PARALLEL'

      real visc(lxyz,ieg1-ieg0+1), w(lxyz*(ieg1-ieg0+1))

      call rzero(visc,lxyz*(ieg1-ieg0+1))

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
      subroutine setgg(gram,u,mass,w1,w2,ns,nsg,n,ndim)

      real gram(ns,nsg),u(n,ndim,ns),mass(n),w1(n,ndim,ns),w2(n,ndim,ns)

      write (6,*) 'ns,nsg,n,ndim',ns,nsg,n,ndim

      call copy(w1,u,n*ndim*ns)
      call copy(w2,u,n*ndim*ns)

      do i=1,ns*ndim
         call col2(w2(1,i,1),mass,n)
      enddo

      do ioff=0,nsg/ns-1 ! assume ns is the same across all processors
         do js=1,ns
         do is=1,ns
            gram(is,js+ns*ioff)=gram(is,js+ns*ioff)
     $         +vlsc2(w1(1,1,is),w2(1,1,js),n*ndim)
         enddo
         enddo
         
         call shift(w1,n*ndim*ns)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine shift_setup(iwk1,iwk2,n)

      include 'SIZE'

      common /pcomm/ igsh1,igsh2
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer*8 iwk1(n),iwk2(n)

      igsh1=0
      igsh2=0

      do i=1,n
         i1=nid/2
         i2=mod(nid-1+mp,mp)/2
         irem=mod(nid,2).eq.0
         iwk1(i)=i+i1*n
         iwk2(i)=i+i2*n
      enddo

      do i=1,mp
         if (nid.eq.(i-1)) then
            do j=1,n
               write (6,*) nid,i,j,iwk1(j),iwk2(j)
            enddo
         endif
         call nekgsync
      enddo

      ntot=(mp*n)/2

      call fgslib_gs_setup(igsh1,iwk1,ntot,nekcomm,mp)
      call fgslib_gs_setup(igsh2,iwk2,ntot,nekcomm,mp)

      return
      end
c-----------------------------------------------------------------------
      subroutine shift(u,w1,w2,n)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /pcomm/ igsh1,igsh2

      real u(n),w1(n),w2(n)

      call copy(w1,u,n)
      call copy(w2,u,n)

      call fgslib_gs_op(igsh1,w1,1,1,0)
      call fgslib_gs_op(igsh2,w2,1,1,0)

      if (mod(mid,2).eq.0) then
         call sub2(u,w2,n)
         call chsign(u,n)
      else
         call sub2(u,w1,n)
         call chsign(u,n)
      endif

      return
      end
c-----------------------------------------------------------------------
c     subroutine ax(au,u,helm1,helm2,imesh,isd)
      subroutine aop(au,u,visc,gfac,imesh,isd,nel)
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
c
      common /fastax/ wddx(lx1,lx1),wddyt(ly1,ly1),wddzt(lz1,lz1)
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv
c
      real           au    (lx1,ly1,lz1,1)
     $ ,             u     (lx1,ly1,lz1,1)
     $ ,             helm1 (lx1,ly1,lz1,1)
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

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nxy=lx1*ly1
      nyz=ly1*lz1
      nxz=lx1*lz1
      nxyz=lx1*ly1*lz1
      ntot=nxyz*nel

      if (.not.ifsolv) call setfast(helm1,helm2,imesh)
      call rzero (au,ntot)

      do ie=1,nel
 
        if (ifaxis) call setaxdy ( ifrzer(ie) )

        if (ldim.eq.2) then

c       2-d case ...............

           if (iffast(ie)) then

c          Fast 2-d mode: constant properties and undeformed element
 
           h1 = visc(1,1,1,ie)
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           call mxm   (u(1,1,1,ie),lx1,wddyt,ly1,tm2,ly1)
           call col2  (tm1,g4m1(1,1,1,ie),nxyz)
           call col2  (tm2,g5m1(1,1,1,ie),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)
 
           else
 
           call mxm  (dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           call mxm  (u(1,1,1,ie),lx1,dytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,g1m1(1,1,1,ie),nxyz)
           call col3 (tmp2,duds,g2m1(1,1,1,ie),nxyz)
           if (ifdfrm(ie)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,ie),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,ie),nxyz)
           endif
           call col2 (tmp1,helm1(1,1,1,ie),nxyz)
           call col2 (tmp2,helm1(1,1,1,ie),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,ie),tm1,nxyz)
           call add2 (au(1,1,1,ie),tm2,nxyz)

        endif

        else

c       3-d case ...............
 
           if (iffast(ie)) then
 
c          Fast 3-d mode: constant properties and undeformed element
 
           h1 = helm1(1,1,1,ie)
           call mxm   (wddx,lx1,u(1,1,1,ie),lx1,tm1,nyz)
           do 5 iz=1,lz1
           call mxm   (u(1,1,iz,ie),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
 5         continue
           call mxm   (u(1,1,1,ie),nxy,wddzt,lz1,tm3,lz1)
           call col2  (tm1,g4m1(1,1,1,ie),nxyz)
           call col2  (tm2,g5m1(1,1,1,ie),nxyz)
           call col2  (tm3,g6m1(1,1,1,ie),nxyz)
           call add3  (au(1,1,1,ie),tm1,tm2,nxyz)
           call add2  (au(1,1,1,ie),tm3,nxyz)
           call cmult (au(1,1,1,ie),h1,nxyz)
 
           else
 
 
           call mxm(dxm1,lx1,u(1,1,1,ie),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,ie),lx1,dytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,ie),nxy,dztm1,lz1,dudt,lz1)
           call col3    (tmp1,dudr,g1m1(1,1,1,ie),nxyz)
           call col3    (tmp2,duds,g2m1(1,1,1,ie),nxyz)
           call col3    (tmp3,dudt,g3m1(1,1,1,ie),nxyz)
           if (ifdfrm(ie)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,ie),nxyz)
              call addcol3 (tmp1,dudt,g5m1(1,1,1,ie),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,ie),nxyz)
              call addcol3 (tmp2,dudt,g6m1(1,1,1,ie),nxyz)
              call addcol3 (tmp3,dudr,g5m1(1,1,1,ie),nxyz)
              call addcol3 (tmp3,duds,g6m1(1,1,1,ie),nxyz)
           endif
           call col2 (tmp1,helm1(1,1,1,ie),nxyz)
           call col2 (tmp2,helm1(1,1,1,ie),nxyz)
           call col2 (tmp3,helm1(1,1,1,ie),nxyz)
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
         do ie=1,nel
 
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
     $                          + helm1(i,j,1,ie)*(term1+term2)
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
