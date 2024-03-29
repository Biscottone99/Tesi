program flash
use module

  implicit none
  character*1:: jobz,uplo
  integer:: lda, lwork,isito,irwork,liwork,info, lrwork, si, sj,nso,conta, count
  integer,allocatable::iwork(:)
  integer::nsiti,dim2, i, n, a, b, c, d, p,j, k, l
  integer::  temp, m, contasito
  integer,allocatable:: vecconfig(:), occupazioni(:), NZ(:), spar(:,:)
  real*8,allocatable:: rwork(:), w(:), dist(:,:,:), hop(:,:), nuclei(:,:), hop2(:,:), spintot(:), carica(:,:), u(:),esite(:), mu(:,:,:), charges(:,:), dipole(:,:), mualpha(:,:), dist2(:,:), mubeta(:,:)
  real*8, allocatable:: muralpha(:,:,:), murbeta(:,:,:)
  complex*16,allocatable::ham(:,:),work(:), hamsoc(:,:),  soc_a(:,:,:), soc_b(:,:,:),soc_mono(:,:,:), pp(:,:),coup(:,:), COUPLING(:,:),pp2(:,:,:,:),  pp2r(:,:,:,:), now(:,:), now2(:,:), ppso(:,:,:), s(:,:,:), ssq(:,:,:),srot(:,:,:), hopalpha(:,:,:), hopbeta(:,:,:), ppa(:,:,:), ppb(:,:,:), ppra(:,:,:), pprb(:,:,:)
  real*8,allocatable:: dsite(:,:), ssite(:), spin3(:), spin2(:), pol(:,:), polr(:,:), tt(:,:), pot(:), energy(:), sqrot(:,:)
  real*8:: Uc, t, PPP, me, gs, e, e0, pi, cl, radius
  logical:: bool, bool1, bool2, bool3, is_hermitian
  real*8::strenght,hbar,hbarev,echarge,emass,dpg, bubu, balu, delta, gamma, gamma2, tollerance
  complex*16::sy(2,2), spin(2,2,3),vec1(3),vec2(3), cp(3)
  complex*16::cplx,pf
  character*1,allocatable::state(:)
  
 nsiti=4
 
  
  Uc=6d0
  t=0.6d0
  nso=nsiti*2
  me=9.1093837015d-31
  gs=2.00231930436256
  e=1.602176634d-19
  e0=8.8541878128d-12
  pi=dacos(-1.d0)
!  write(*,*) pi
  cl=299792458 
  cplx=cmplx(0.d0, 1.d0)


  pf=((gs*e**2)/(8*pi*e0*me*cl**2))*10d10
  write(*,*) pf

  open(1,file='basis.dat')
  open(2,file='geom.dat')
  open(3,file='hamiltonian.dat')
  open(4,file='results.dat')
 ! open(5,file='input.dat')
  open(6,file='check.dat')
  open(10,file='dim2.dat')
  read(10,*) dim2
  close(10)
  

  allocate(vecconfig(dim2), nuclei(nsiti,3),ham(dim2,dim2), occupazioni(nsiti),dist(nsiti,nsiti,3), soc_a(3,nso,nso), soc_b(3,nso, nso), hop(nsiti,nsiti),hop2(nso,nso),soc_mono(3,nso,nso),hamsoc(nso,nso), spintot(dim2), carica(dim2,nsiti),nz(nsiti), u(nsiti), esite(nsiti), coup(dim2,dim2), COUPLING(DIM2, dim2), spar(dim2,dim2), spin3(dim2), spin2(dim2), state(dim2), dipole(dim2,3), MUalpha(DIM2,3),ppso(3,nso,nso),tt(dim2,dim2), pot(dim2), energy(dim2),s(3,nso,nso), mubeta(dim2,3),ppa(3,dim2,dim2), ppb(3,dim2,dim2))
!=========================LETTURA INPUT=========================
  do i=1,dim2
     read(1,*) vecconfig(i)
  enddo
  spin=0
  spin(1,2,1)=1d0
  spin(2,1,1)=1d0

  spin(1,2,2)=cplx
  spin(2,1,2)=-cplx

  spin(1,1,3)=1d0
  spin(2,2,3)=-1d0
  do i=1,nsiti
     u(i)=uc
  enddo
  esite=0
  esite(2)=0!2d0
  esite(3)=esite(2)
  esite(1)=-4d0
  esite(4)=+4d0

  do i=2,3
     nz(i)=1
  enddo
  nz(4)=0
  nz(1)=2
 
  do i=1,nsiti
     read(2,*) nuclei(i,1), nuclei(i,2), nuclei(i,3)
  enddo

  do i=1,nsiti
     write(6,*) i, nz(i), u(i), esite(i)
  enddo
  allocate(dist2(nsiti, 3))
  call calcola_distanze(nuclei, dist2, nsiti)
  write(6,*) 'DISTANCES'
  call write_matrix(dist2, 6, nsiti,nsiti, nsiti)

  
  !=========================CHARGES AND MOMENTUM=========================
  call charge(carica, vecconfig, nz, dim2, nso)
  call dipole_moment(dipole, carica, nuclei, dim2, nsiti)
  

  mualpha=0
  do n=1,dim2
     do i=0,nso-2,2
        isito=(i+2)/2
        bool=btest(vecconfig(n), i)
        do k=1,3
           if(bool)mualpha(n,k)=mualpha(n,k)-nuclei(isito,k)
        enddo
     enddo
  enddo
  
  mubeta=0
  do n=1,dim2
     do i=1,nso-1,2
        isito=(i+1)/2
        bool=btest(vecconfig(n), i)
        do k=1,3
           if(bool)mubeta(n,k)=mubeta(n,k)-nuclei(isito,k)
        enddo
     enddo
  enddo
  !=========================Spin operator=========================
  s=0
  do k=1,3
     do i=1,nso-1,2
        s(k,i,i)=spin(1,1,k)
        s(k,i+1,i+1)=spin(2,2,k)
        s(k,i,i+1)=spin(1,2,k)
        s(k,i+1,i)=spin(2,1,k)
     enddo
  enddo
!!$  do k=1,3
!!$     write(6,*) 'REAL,k=',k
!!$     do i=1,nso
!!$        write(6,'(<nso>(2x,f10.5))') (dreal(s(k,i,j)), j=1,nso)
!!$     enddo
!!$     write(6,*) 'IMAG,k=',k
!!$     do i=1,nso
!!$        write(6,'(<nso>(2x,f10.5))') (dimag(s(k,i,j)), j=1,nso)
!!$     enddo
!!$  enddo
  
  allocate (ssq(3,dim2,dim2))
  do k=1,3
     call sq_oe_op_compl(nso,dim2,s(k,:,:),ssq(k,:,:),vecconfig)
     call check_hermitian(ssq(k,:,:), dim2, is_hermitian)
     if(.not.is_hermitian)write(*,*) 'Not hermitian k=',k
  enddo
  hop2=0
  do i=1,nso-2
     hop2(i,i+2)=t
  enddo
  call copy_upper_to_lower(hop2, nso)
  call  sq_oe_op_real(nso,dim2,hop2,tt,vecconfig) 

  hop=0
  do i=1,nsiti-1
     hop(i,i+1)=t
  enddo
  call copy_upper_to_lower(hop, nsiti)
  call momentum_so (nsiti, nso, hop, nuclei, ppso)

  allocate(hopalpha(3,nso,nso), hopbeta(3,nso,nso))
  hopalpha=0
  do k=1,3
     do i=1,nso-3,2
        hopalpha(k,i,i+2)=ppso(k,i,i+2)
     enddo
  enddo
  call copia_comp_conj(hopalpha, 3, nso)


  hopbeta=0
  do k=1,3
     do i=2,nso-2,2
        hopbeta(k,i,i+2)=ppso(k,i,i+2)
     enddo
  enddo
  call copia_comp_conj(hopbeta, 3, nso)

  do k=1,3
     call check_hermitian(hopalpha(k,:,:), nso, is_hermitian)
     call check_hermitian(hopbeta(k,:,:), nso, bool1)
     if(.not.is_hermitian)write(*,*) 'alpha k=',k
     if(.not.bool1)write(*,*) 'beta k=', k
  enddo
  do k=1,3
     call sq_oe_op_compl(nso, dim2, hopalpha(k,:,:),ppa(k,:,:),vecconfig)
     call check_hermitian(ppa(k,:,:), dim2, bool)
     if(.not.bool)write(*,*) 'ERROR PPA K=',k
      call sq_oe_op_compl(nso, dim2, hopbeta(k,:,:),ppb(k,:,:),vecconfig)
     call check_hermitian(ppb(k,:,:), dim2, bool1)
     if(.not.bool1)write(*,*) 'ERROR PPb K=',k
  enddo
  
!=========================SCRITTURA HAM=========================

  ham=0
  
  call ohno (dim2,nsiti,nuclei,vecconfig, u,nz, pot)
  call site_energy(nso,dim2,esite,vecconfig,energy)
 
  !=========================SPIN-ORBIT COUPLING=========================
  soc_a=0
  SOC_B=0
  soc_mono=0
  do i=1,nso
     do j=1,nso
        do isito = 1,nsiti
           if (isito.ne.(i+1)/2) then 
              ! r x p cross product
              vec1 = (0.d0, 0.d0)
              vec2 = (0.d0, 0.d0)
              cp=0
              do k = 1,3
                 vec1(k) = nuclei((i+1)/2,k) - nuclei(isito,k) !position
                 vec2(k) = ppso(k,(i+1)/2,(j+1)/2) !momentum
              end do
              cp  = cross_product(vec1, vec2)

              ! 1/|r_aA|^3 term
              radius = 0.d0
              do k = 1,3
                 radius = radius + dreal(vec1(k))**2
              end do
              radius = (dsqrt(radius))**3

              cp = cp / radius
              if(i/2*2.eq.i)then
                 si = 2
              else
                 si=1
              endif

              if(j/2*2.eq.j)then
                 sj = 2
              else
                 sj=1
              endif

              do k = 1,3
                 soc_a(k,i,j) = soc_a(k,i,j) + cp(k)  * spin(si,sj,k)
              end do
           end if
        end do

     end do
  end do

  do i=1,nso
     do j=1,nso
        do isito = 1,nsiti
           if(isito.ne.(j+1)/2)then 
              ! r x p cross product
              vec1 = (0.d0, 0.d0)
              vec2 = (0.d0, 0.d0)
              cp=0
              do k = 1,3
                 vec1(k) = nuclei((j+1)/2,k) - nuclei(isito,k) !position
                 vec2(k) = ppso(k,(i+1)/2,(j+1)/2) !momentum
                ! write(*,*) vec2             
              end do
              cp  = cross_product(vec1, vec2)
              !write(*,*) cp

              ! 1/|r_aA|^3 term
              radius = 0.d0
              do k = 1,3
                 radius = radius + dreal(vec1(k))**2
              end do
              radius = (dsqrt(radius))**3

              cp = cp / radius
              
              if(i/2*2.eq.i)then
                 si = 2
              else
                 si=1
              endif

              if(j/2*2.eq.j)then
                 sj = 2
              else
                 sj=1
              endif

              do k = 1,3
                 soc_b(k,i,j) = soc_b(k,i,j) + cp(k)  * spin(si,sj,k)
                
              end do
           end if
        end do

     end do
  end do

  do i = 1,nso
     do j = 1,nso
        do k = 1,3
           soc_mono(k,i,j)    = 0.5d0 * (soc_a(k,i,j) + soc_b(k,i,j))
        end do
     end do
  end do
  

  hamsoc=0d0
  do i=1,nso
     do j=1,nso
        do k=1,3
           hamsoc(i,j)=hamsoc(i,j)+soc_mono(k,i,j)
        enddo
     enddo
  enddo
  write(6,*) 'SOC MONO'
  do n=1,nso
     do m=1,nso
        if(hamsoc(n,m).ne.0)write(6,*) hamsoc(n,m)
     enddo
  enddo
        

  
  call  check_hermitian(soc_mono, nso, is_hermitian)
  if(is_hermitian) write(*,*) 'soc_mono HERMITIANA'
  call  check_hermitian(hamsoc, nso, is_hermitian)
  if(is_hermitian) write(*,*) 'hamsoc_HERMITIANA'
  if(.not.is_hermitian) write(*,*) 'PROBLEMI'

  coup=0
  call sq_oe_op_compl(nso,dim2,hamsoc,coup,vecconfig)
  call  check_hermitian(coup, dim2, is_hermitian)
  if(is_hermitian) write(*,*) 'COUP HERMITIANA'
  
  coup=pf*coup
 
  !============================================================

  !! Adding to hamiltonian
  do n=1,dim2
     ham(n,n)=pot(n)+energy(n)
     do m=1,dim2
        ham(n,m)=ham(n,m)+coup(n,m)-tt(n,m)
     enddo
  enddo
  call  check_hermitian(ham, dim2, is_hermitian)
  if(is_hermitian) write(*,*) 'HAM HERMITIANA'
  jobz ='V'
  uplo='U'
  lrwork=(1+5*dim2+2*dim2**2)
  liwork=(3+5*dim2)
  lwork=(2*dim2+dim2**2)
  allocate(w(dim2),work(max(1,lwork)),rwork(lrwork),iwork(max(1,liwork)))
!==============================DIAGONALIZATION=========================
  call zheevd (jobz, uplo, dim2, ham, dim2, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  call eigenvalues(dim2,1d-10,w,state)
  
  allocate(srot(3,dim2,dim2),sqrot(dim2,dim2))
  do n=1,dim2
     do m=1,dim2
        do k=1,3
           ssq(k,n,m)=ssq(k,n,m)*dconjg(ssq(k,n,m))
        enddo
     enddo

  enddo
  do n=1,dim2
     do m=1,dim2
        do k=1,3
           sqrot(n,m)=sqrot(n,m)+ssq(k,n,m)
        enddo
     enddo
  enddo
  call rotate_real_2x2(dim2,sqrot,sqrot,ham)


  
  write(4,*) 'EIGENVALUES'
  do i=1,12
     write(4,*) w(i)-w(1), state(i)
  enddo

  write(6,*) 'S^2'

    do n=1,dim2
     do m=1,dim2
        if(sqrot(n,m).ne.0)write(4,*) n, m, sqrot(n,m)
     enddo
  enddo
 
  
  !=========================OUTPUT=========================
  allocate(charges(dim2,nsiti))
  charges=0
  call rot_diag(dim2,charges,carica,nsiti,ham)
  
  write(4,*) 'CARICHE'
  call write_matrix(charges, 4, dim2, nsiti, 10)

  allocate(mu(dim2,dim2,3), pp2r(dim2,dim2,2,3), muralpha(dim2,dim2,3), murbeta(dim2,dim2,3))
  murbeta=0
  muralpha=0

  do i=1,dim2 !a
     do l=1,dim2 !b
        do j=1,dim2 !alfa
           do k=1,3
              muralpha(i,l,k)=muralpha(i,l,k)+dconjg(ham(j,i))*ham(j,l)*mualpha(j,k)
              murbeta(i,l,k)=murbeta(i,l,k)+dconjg(ham(j,i))*ham(j,l)*mubeta(j,k)
           enddo
        enddo
     enddo
  enddo


  
  write(4,*) 'DIPOLE ALPHA-BETA'
  write(4,'(<3>(2x,f10.5))') (muralpha(1,1,k)-murbeta(1,1,k), k=1,3)
  write(4,*) 'DIPOLE ALPHA+BETA'
  write(4,'(<3>(2x,f10.5))') (muralpha(1,1,k)+murbeta(1,1,k), k=1,3)
  mu=0
  do i=1,dim2 !a
     do l=1,dim2 !b
        do j=1,dim2 !alfa
           do k=1,3
              mu(i,l,k)=mu(i,l,k)+dconjg(ham(j,i))*ham(j,l)*dipole(j,k)
           enddo
        enddo
     enddo
  enddo


  
!!$  do k=1,3
!!$     call rotate_real_1x2(dim2,mu(:,:,k),dipole(:,k),ham)
!!$  enddo
  
  allocate(ppra(3,dim2,dim2),pprb(3,dim2,dim2))
  do k=1,3
     call rotate_cplx_2x2(dim2,ppra(k,:,:),ppa(k,:,:),ham)
     call rotate_cplx_2x2(dim2,pprb(k,:,:), ppb(k,:,:), ham)
  enddo
  
  write(6,*) "DIPOLE MOMENT"
  do a=1,3
     write(6,*) 'DIPOLE k=',a
     call write_matrix (mu(:,:,a), 6, dim2,dim2, dim2)
  enddo

  call rotate_cplx_2x2(dim2,coupling,coup,ham)
  write(4,*) 'COUPLING'

  do n=1,dim2
     do m=1,dim2
        if(coup(n,m).ne.0)write(4,*) coup(n,m), n, m
     enddo
  enddo

  do k=1,3
     call check_hermitian (ppra(k,:,:), dim2, bool)
     if(.not.bool)write(*,*) 'ERROR PPRA K=',k
     call check_hermitian (pprb(k,:,:), dim2, bool)
     if(.not.bool)write(*,*) 'ERROR PPRB K=',k
     
  enddo
 
  end program
