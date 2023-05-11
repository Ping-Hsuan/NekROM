       program EOFF

       include 'SIZE'
       include 'OFFLINE'

       common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
       integer nb,mel
       logical iftherm,ifdebug,ifgf,ifbuoy,ifm1,ifavg0
       integer comm_out

       call offline_init(comm_out)
        
       write(6,*)'Checking Parameters'
       write(6,*)'lxyz:',lxyz
       write(6,*)'lel',lel
       write(6,*)'lsg',lsg
       write(6,*)'leg',leg
       write(6,*)'lgld',lgld
       write(6,*)'lx1,ly1,lz1',lx1,ly1,lz1
       write(6,*)'mp',mp

       nb = 20
       mel = 10
       nsg = 20 ! number of ls 
       neg = 512
       ifavg0=.true.
       ifm1=.true.
       
       ! set gvec to be identitiy matrix
       call gengrams1(ga,gb,gc,gm,gf,gt,buf,tmpf,gvec,ieg,indxr,
     $     nsg,nb,mp,neg,mel,ldim,lx1,lxyz,iftherm,ifbuoy,ifgf,ifm1)

       ! remove 0th mode contribution, also support arbitrary lifting functions
       ! gb still has 0th mode contribution
       ! gg is the one that does not have 0th mode contribution
       call setg(gg,gb,gt,nsg,ifavg0) 
       
       call write_ops(ga,gg,gc,nsg**2,nsg**3,mid,' gu',.false.)

       ! compute eigenvector
       call setq(gvec,gvect,gvecc,gval,gg,nsg,nsc,mp,mid,ifavg0,nsg1)
       write(6,*)'nsg1',nsg1
       write(6,*)'nsc',nsc

       do j=1,nsg
       do i=1,nsg
          write (6,*) i,j,gvec(i+(j-1)*nsg,1),'q'
       enddo
       enddo

       ! Form the rom operators
       call qop2(ga,gt,gvec,gvect,nsg,nsg1)
       call qop2(gb,gt,gvec,gvect,nsg,nsg1)
       call write_ops(ga,gb,gc,(nsg1)**2,nsc*(nsg1)**2,
     $                mid,'  u',.false.)
       call qop3(gc,gt,gvec,gvect,gvecc,nsg,nsg1,nsc,mid,'u')

c      if(mid.eq.0) then
c      do i=1,nb*nb
c      write(6,*)'gram checking',i,gb(i)
c      enddo
c      endif

c      call genevec(gvec,gval,gb,nsg,1)
c      if(mid.eq.0) then
c      do i=1,nb*nb
c      write(6,*)'gram checking again',i,gb(i)
c      enddo
c      endif

c      if(mid.eq.0) then
c      do i=1,nb
c      write(6,*)'gval checking',i,gval(i)
c      enddo
c      do i=1,nb*nb
c      write(6,*)'gvec checking',i,gvec(i,1)
c      enddo
c      endif

c      call transpose(gvect,nb,gvec,nsg)
c      do i=1,nb*nb
c      write(6,*)'gvec checking',i,gvect(i,1)
c      enddo
c      call qop2(gb,gt,gvec,gvect,nsg,nb)
c      if (mid.eq.0) then
c      do i=1,nb*nb
c      write(6,*)'B_rom checking',i,gb(i)
c      enddo
c      endif

       call exitt0

       end
