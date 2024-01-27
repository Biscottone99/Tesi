program hamiltonian
  implicit none
  integer:: a, b, c, i, j, k, nsiti,nso, count, ne, n, lda, lwork, nf,ldvl,ldvr,info
  integer,allocatable:: config(:)
  real*8::ax, ay, az, temp, uc,energy, t, temp2, hbar, sx(2,2), sz(2,2), cl, temp66,me
  real*8, allocatable:: cord(:,:),rad(:,:,:), pot(:,:), charge(:,:), opera(:,:),eig(:,:),opera2(:,:), rwork(:)
  logical::bool, bool1
  character*1:: jobvl, jobvr
  complex*16::cplx,pf
  complex*16::sy(2,2), spin(2,2)
  complex*16,allocatable::opera3(:,:), ham(:,:),  mom(:,:), w(:),vl(:,:), vr(:,:),work(:)
  
  nsiti=4
  hbar=1!6.58d-16
  me=1!9.11d-31
  cl=1!299792458
  ne=2
  nso=nsiti*2
  nf=28
  cplx=cmplx(0.d0, 1.d0)

 temp66=me*hbar*cl*cl*2d0
!=========================PARAMETERS=========================
  uc=11.26d0
  t=2.50d0
  pf=cplx!(cplx*t)/(temp66)
 
!============================================================  
!  write(*,*) pf, temp66
  allocate (cord(nsiti,3),config(nf), charge(nf,nsiti), pot(nsiti,nsiti), rad(3,nsiti,nsiti),opera(nso,nso), ham(nso,nso), eig(nso,nso),opera2(nso,nso), opera3(nso,nso), mom(nso,nso))


  open(1,file='geom.dat')
  open(2,file='configurations.dat')
  open(3, file='results.dat')
  open(4,file='operatore.dat')
  open(5,file='ham.dat')
  open(6,file='check.dat')


  do i=1,4
     read(1,*) cord(i,1), cord(i,2), cord(i,3)
  enddo
  close(1)

!!$  do i=1, nf
!!$     read(2,*) config(i)
!!$  enddo
!!$  do k=1,nf
!!$     count=0
!!$     do i=0,nso-1,2
!!$        bool=btest(config(k),i)
!!$        bool1=btest(config(k),i+1)
!!$        if(bool)count=count+1
!!$        if(bool1)count=count+1
!!$        charge(k,(i+2)/2)=1d0-count
!!$     enddo
!!$  enddo



!=========================PPP-PART=========================
  ham=0

!!!! potential term - Ohno parametrization 
  do i=2,nso,2
     temp=0
     do j=2,nso,2
        temp=temp+(14.397)/((28.794/(2*uc)**2+((cord(i/2,1)-cord(j/2,1))**2+(cord(i/2,2)-cord(j/2,2))**2)+(cord(i/2,3)-cord(j/2,3))**2))**0.5
     enddo
     ham(i,i)=0.5*temp
     ham(i-1,i-1)=0.5*temp
  enddo




!!! off-diagonal part
  do i=1,nso-3,2
     ham(i,i+2)=-t
     ham(i+1,i+3)=-t
  enddo
  do i=1,nso
     do j=1,nso
        if(i.ne.j)then
           ham(j,i)=ham(i,j)
        endif
     enddo
  enddo


!============================================================

  !!! supporting operator - distances between different sites for each dimension (x, y, z)
  do i=1,4
     do j=1,4
        do k=1,3
           rad(k,i,j)=cord(i,k)-cord(j,k)
        enddo
     enddo
  enddo

  !pauli matrix
  sx(1,2)=hbar/2
  sx(2,1)=sx(1,2)

  sy(1,2)=-hbar/(2*cplx)
  sy(2,1)=-sy(1,2)
  
  sz(1,1)=hbar/2
  sz(2,2)=hbar/2

  do i=1,2
     do j=1,2
        spin(i,j)=sx(i,j)+sy(i,j)+sz(i,j)
     enddo
  enddo

  do i=1,2
     do j=1,2
        write(4,*) spin(i,j)
     enddo
  enddo
 
!=========================SOC PART=========================

  opera=0

!!Opera is the matrix containing the vectorials product of nearest-neighbors distances written on the atomic orbitals  
  do i=1,nsiti-1
     do j=1,nsiti
        if(i.ne.j)then
           ax=rad(2,i,j)*rad(3,i,i+1)-rad(3,i,j)*rad(2,i,i+1)
           ay=rad(3,i,j)*rad(1,i,i+1)-rad(3,i,i+1)*rad(1,i,j)    
           az=rad(1,i,j)*rad(2,i,i+1)-rad(2,i,j)*rad(1,i,i+1)

           temp=0
           do k=1,3
              temp=temp+rad(k,i,j)**2
           enddo

           opera(i,i+1)=opera(i,i+1)+(dsqrt(temp)**(-3))*dsqrt((ax**2+ay**2+az**2))
           temp2=(dsqrt(temp)**(-3))*dsqrt((ax**2+ay**2+az**2))
        endif
     enddo
  enddo

  do i=2,nsiti
     do j=1,nsiti
        if(i.ne.j)then
           ax=rad(2,i,j)*rad(3,i,i-1)-rad(3,i,j)*rad(2,i,i-1)
           ay=rad(3,i,j)*rad(1,i,i-1)-rad(3,i,i-1)*rad(1,i,j)    
           az=rad(1,i,j)*rad(2,i,i-1)-rad(2,i,j)*rad(1,i,i-1)

           temp=0
           do k=1,3
              temp=temp+rad(k,i,j)**2
           enddo

           opera(i,i-1)=opera(i,i-1)+(dsqrt(temp)**(-3))*dsqrt((ax**2+ay**2+az**2))
           temp2=(dsqrt(temp)**(-3))*dsqrt((ax**2+ay**2+az**2))
           write(6,*) i, i-1, j,  temp2
        endif
     enddo
  enddo

  
  !!Opera is the matrix containing the vectorials product of nearest-neighbors distances written on the atomic spin-orbitals  
  do i=1,nsiti
     do j=1,nsiti
        opera2((2*i)-1,2*j-1)=opera(i,j)
        opera2(2*i,2*j)=opera(i,j)
     enddo
  enddo

!!!=========================VECTORIAL PRODUCT SYMMETRIZED
 !!Mom is opera2 symmetrixed
  do i=1,nso-3,2
     mom(i,i+2)=0.5d0*(opera2(i,i+2)+opera2(i+2,i))
     mom(i+1,i+3)=0.5d0*(opera2(i+1,i+3)+opera2(i+3,i+1))
     mom(i+2,i)=0.5d0*(opera2(i+2,i)+opera2(i,i+2))
     mom(i+3,i+1)=0.5d0*(opera2(i+3,i+1)+opera2(i+1,i+3))
  enddo

  opera3=0
  do i=1,nso-3,2
     opera3(i,i+2)=mom(i,i+2)*spin(1,1)
     opera3(i,i+3)=mom(i,i+3)*spin(1,2)
     opera3(i+1,i+2)=mom(i+1,i+2)*spin(2,1)
     opera3(i+1,i+3)=mom(i+1,i+3)*spin(2,2)
     
  enddo
  do i=3,nso-1,2
     opera3(i,i-2)=mom(i,i-2)*spin(1,1)
     opera3(i,i-1)=mom(i,i-1)*spin(1,2)
     opera3(i+1,i-2)=mom(i+1,i-2)*spin(2,1)
     opera3(i+1,i-1)=mom(i+1,i-1)*spin(2,2)
  enddo

  do i=1,nso
     do j=1,nso
        ham(i,j)=ham(i,j)+pf*opera3(i,j)
     enddo
  enddo

  write(5,*) 'REAL HAMILTONIAN'
  do i=1,nso
     write(5,'(<nso>(2x,f10.5))')(dreal(ham(i,j)),j=1,nso)
  enddo

  write(5,*) 'IMAGINARY HAMILTONIAN'
  do i=1,nso
     write(5,'(<nso>(2x,f10.5))')(dimag(ham(i,j)),j=1,nso)
  enddo
  
  close(4)

  ! =========================START DIAGONALIZATION =========================
  jobvl='V'
  jobvr='V'
  lda=nso
  ldvl=nso
  ldvr=nso
  lwork=max(1,2*nso)
  
  allocate(w(nso),vl(ldvl,nso),vr(ldvr,nso),work(max(1,lwork)),rwork(2*nso))

  call zgeev (jobvl, jobvr, nso, ham,nso,w,vl,ldvl,vr,ldvr,work,lwork, rwork,info)

  if(info.eq.'0')write(*,*) 'gle chaveda'

  write(3,*) 'EIGENVALUES'
  do i=1,nso
     write(3,*) dreal(w(i)), dimag(w(i))
  enddo

 write(3,*) 'REAL EIGENVECTORS'
  do i=1,nso
     write(3,'(<nso>(2x,f10.5))')(dreal(ham(i,j)),j=1,nso)
  enddo

  write(3,*) 'IMAGINARY EIGENVECTORS'
  do i=1,nso
     write(3,'(<nso>(2x,f10.5))')(dimag(ham(i,j)),j=1,nso)
  enddo


end program hamiltonian
