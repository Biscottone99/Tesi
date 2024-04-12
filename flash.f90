program flash
use module

  implicit none
  character*1:: jobz,uplo
  integer:: lda, lwork,isito,irwork,liwork,info, lrwork, si, sj,nso,conta, count,perm,phase,bbbb
  integer,allocatable::iwork(:)
  integer::nsiti,dim2, i, n, a, b, c, d, p,j, k, l
  integer::  temp, m, contasito,sitoi, sitoj
  integer,allocatable:: vecconfig(:), occupazioni(:), NZ(:), spar(:,:)
  real*8,allocatable:: rwork(:), w(:), dist(:,:,:), hop(:,:), nuclei(:,:), hop2(:,:), spintot(:), carica(:,:), u(:),esite(:), mu(:,:,:), charges(:,:), dipole(:,:), mualpha(:,:), dist2(:,:), mubeta(:,:)
  real*8, allocatable:: muralpha(:,:,:), murbeta(:,:,:), singlet(:), triplet(:),quintet(:),w2(:), sr(:,:), tr(:,:), qr(:,:), eig(:,:), spindensity(:,:), sdr(:,:)
  complex*16,allocatable::ham(:,:),work(:), hamsoc(:,:),  soc_a(:,:,:), soc_b(:,:,:),soc_mono(:,:,:), pp(:,:),coup(:,:), COUPLING(:,:),pp2(:,:,:,:),  pp2r(:,:,:,:), now(:,:), now2(:,:), ppso(:,:,:), s(:,:,:), ssq(:,:,:),srot(:,:,:), hopalpha(:,:,:), hopbeta(:,:,:), ppa(:,:,:), ppb(:,:,:), ppra(:,:,:), pprb(:,:,:), hssotb(:,:,:,:,:),ssotb(:,:,:,:),sso(:,:),mom(:,:,:), ham2(:,:), sqrot(:,:),sqrot2(:,:), hsootb(:,:,:,:,:), sootb(:,:,:,:), soo(:,:), soc(:,:), socr(:,:), mono_coup(:,:,:), bi_coup(:,:,:), temporary(:,:,:,:,:), mcrot(:,:,:), bcrot(:,:,:)
  real*8,allocatable:: dsite(:,:), ssite(:), spin3(:), spin2(:), pol(:,:), polr(:,:), tt(:,:), pot(:), energy(:)
  real*8:: Uc, t, PPP, me, gs, e, e0, pi, cl, radius
  logical:: bool, bool1, bool2, bool3, is_hermitian
  real*8::strenght,hbar,hbarev,echarge,emass,dpg, bubu, balu, delta, gamma, gamma2, tollerance, bibi
  complex*16::sy(2,2), spin(2,2,3),vec1(3),vec2(3), cp(3), cp1(3)
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
  open(11,file='work.dat')
  read(10,*) dim2
  close(10)
  

  allocate(vecconfig(dim2), nuclei(nsiti,3),ham(dim2,dim2), occupazioni(nsiti), soc_a(3,nso,nso), soc_b(3,nso, nso), hop(nsiti,nsiti),hop2(nso,nso),soc_mono(3,nso,nso),hamsoc(nso,nso), spintot(dim2), carica(dim2,nsiti),nz(nsiti), u(nsiti), esite(nsiti), coup(dim2,dim2), COUPLING(DIM2, dim2), spar(dim2,dim2), spin3(dim2), spin2(dim2), state(dim2), dipole(dim2,3), MUalpha(DIM2,3),ppso(3,nso,nso),tt(dim2,dim2), pot(dim2), energy(dim2),s(3,nso,nso), mubeta(dim2,3),ppa(3,dim2,dim2), ppb(3,dim2,dim2),sqrot(dim2,dim2),sqrot2(dim2,dim2),soc(dim2,dim2),socr(dim2,dim2))
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

  do k=1,3
     write(6,*) 'REAL, K=',k
     call write_matrix(dreal(s(k,:,:)), 6, nso, nso, nso)
      write(6,*) 'IMAG, K=',k
     call write_matrix(dimag(s(k,:,:)), 6, nso, nso, nso)
  enddo
  

  allocate (ssq(3,dim2,dim2))
  do k=1,3
     bbbb=70+k
     call sq_oe_op_compl(nso,dim2,s(k,:,:),ssq(k,:,:),vecconfig)
     call check_hermitian(ssq(k,:,:), dim2, is_hermitian)
     call write_matrix(dreal(ssq(k,:,:)), bbbb, dim2, dim2, dim2)
     call write_matrix(dimag(ssq(k,:,:)), bbbb+10, dim2, dim2, dim2)
     if(.not.is_hermitian)write(*,*) 'Not hermitian k=',k
  enddo
  ssq=0.5*ssq

  do k=1,3
     call square_complex_matrix(dim2, ssq(k,:,:))
     call check_hermitian(ssq(k,:,:), dim2, is_hermitian)
     if(.not.is_hermitian)write(*,*) 'Not hermitian k=',k
  enddo

  sqrot=0
  do n=1,dim2
     do m=1,dim2
        do k=1,3
           sqrot(n,m)=sqrot(n,m)+ssq(k,n,m)
        enddo
     enddo
  enddo

  call check_hermitian(sqrot, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Not hermitian'
  
  hop2=0
  do i=1,nso-2
     hop2(i,i+2)=t
  enddo
  call copy_upper_to_lower(hop2, nso)
  call  sq_oe_op_real(nso,dim2,hop2,tt,vecconfig) 
  write(6,*) 'HOP'
  hop=0
  do i=1,nsiti-1
     hop(i,i+1)=-t
  enddo
  call copy_upper_to_lower(hop, nsiti)
  call momentum_so (nsiti, nso, hop, nuclei, ppso)
  call  write_matrix (hop, 6, nsiti, nsiti, nsiti)

 
  do k=1,3
    !!!! write(6,*) 'PPSO K=',k
    ! call  write_matrix (dimag(ppso(k,:,:)), 6, nso, nso, nso)
     !call check_hermitian(ppso(k,:,:), nso, bool)
     if(.not.bool)write(*,*) 'PROBLEMI PPSO K=',k
  enddo

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
  !=========================SPIN DENSITY=========================
  allocate(spindensity(dim2,nsiti),sdr(dim2,nsiti))
  spindensity=0
  do n=1,dim2
     do i=0,nso-1,2
        isito=(i+2)/2
        bool=btest(vecconfig(n),i)
        if(bool)spindensity(n,isito)=spindensity(n,isito)+1
        bool1=btest(vecconfig(n),i+1)
        if(bool1)spindensity(n,isito)=spindensity(n,isito)-1
     enddo
  enddo
  
!=========================SCRITTURA HAM=========================

  ham=0
  
 ! call ohno (dim2,nsiti,nuclei,vecconfig, u,nz, pot)
  !call site_energy(nso,dim2,esite,vecconfig,energy)
  call site_energy_u(nso,dim2,esite,u,vecconfig,energy)
 
  !=========================SPIN-ORBIT COUPLING=========================
  allocate(mom(nsiti,nsiti,3))
  mom=0
  do i = 1,nsiti
     do j = 1,nsiti
        do k = 1,3
           bubu = nuclei(i,k)-nuclei(j,k)
           mom(i,j,k) = cplx * bubu * hop(i,j)
        end do
     end do
  end do

  do k=1,3
     write(6,*) 'MOM K=',k
     call  write_matrix (dimag(mom(:,:,k)), 6, nsiti, nsiti, nsiti)
     call check_hermitian(mom(:,:,k), nsiti, bool)
     if(.not.bool)write(*,*) 'PROBLEMI MOM K=',k
  enddo
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
              sitoi=(i+1)/2
              sitoj=(j+1)/2
              do k = 1,3
                 vec1(k) = nuclei(sitoi,k) - nuclei(isito,k) !position
                 vec2(k) = mom(sitoi,sitoj,k) !momentum
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
                 !if(soc_a(k,i,j).ne.0)write(77,*) soc_a(k,i,j)
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
              sitoi=(i+1)/2
              sitoj=(j+1)/2
                         
              do k = 1,3
                 vec1(k) = nuclei(sitoj,k) - nuclei(isito,k) !position
                 vec2(k)= mom(sitoi,sitoj,k)  !momentum
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
                 si= 1
              endif

              if(j/2*2.eq.j)then
                 sj = 2
              else
                 sj= 1
              endif

              do k = 1,3
                 soc_b(k,i,j) = soc_b(k,i,j) + cp(k)  * spin(si,sj,k)
                 !if(soc_b(k,i,j).ne.0)write(77,*) soc_a(k,i,j)
              end do
           end if
        end do

     end do
  end do

  do i = 1,nso
     do j = 1,nso
        do k = 1,3
           soc_mono(k,i,j)= 0.5d0 * (soc_a(k,i,j) + soc_b(k,i,j))
        end do
     end do
  end do


  allocate(mono_coup(3,dim2,dim2))
  mono_coup=0
  do k=1,3
     call sq_oe_op_compl(nso,dim2,soc_mono(k,:,:),mono_coup(k,:,:),vecconfig)
  enddo
  

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
  !=========================TWO TERM SS0=========================
  !SSO TERM
  allocate(dist(nsiti,nsiti,k),hssotb(3,nso,nso,nso,nso),ssotb(nso,nso,nso,nso),sso(dim2,dim2))
  dist=0
  do i=1,nsiti
     do j=1,nsiti
        do k=1,3
           dist(i,j,k)=nuclei(i,k)-nuclei(j,k)
        enddo
     enddo
  enddo
  hssotb=0
 
  do a=1,nso
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if((b.eq.d).and.((a+1)/2.ne.(b+1)/2))then
                 vec1=0
                 vec2=0
                 cp=0
                 do k=1,3
                    vec1(k)=dist((a+1)/2,(b+1)/2,k)
                    vec2(k)=ppso(k,a,c)
                 enddo
                 cp=cross_product(vec1,vec2)
                 ! write(*,*) cp
                 ! 1/|r_aA|^3 term
                 radius = 0.d0
                 do k = 1,3
                    radius = radius + dreal(vec1(k))**2
                 end do
                 radius = (dsqrt(radius))**3
                ! write(*,*) radius
                 cp = cp / radius
                 if(mod(a,2).eq.0)then
                    si=2
                 else
                    si=1
                 endif

                 if(mod(c,2).eq.0)then
                    sj=2
                 else
                    sj=1
                 endif
                
                 do k=1,3
                    hssotb(k,a,b,c,d)=hssotb(k,a,b,c,d)+cp(k)*spin(si,sj,k)
                    ! write(*,*) hssotb(k,a,b,c,d)
                 enddo
                 
              endif
           enddo
        enddo
     enddo
  enddo

!!$
  ssotb=0
  allocate(temporary(3,nso,nso,nso,nso))
  temporary=0
  
  do k=1,3
     do a=1,nso
        do b=1,nso
           do c=1,nso
              do d=1,nso
                 ssotb(a,b,c,d)=ssotb(a,b,c,d)+0.5d0*(hssotb(k,a,b,c,d)+dconjg(hssotb(k,c,d,a,b)))
                 temporary(k,a,b,c,d)=temporary(k,a,b,c,d)+0.5d0*(hssotb(k,a,b,c,d)+dconjg(hssotb(k,c,d,a,b)))
              enddo
           enddo
        enddo
     enddo
  enddo

  
  do a=1,nso
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if(zabs(ssotb(a,b,c,d)-dconjg(ssotb(c,d,a,b))).ge.1d-10)write(*,*) a, b, c, d, zabs(ssotb(a,b,c,d)-dconjg(ssotb(c,d,a,b)))
           enddo
        enddo
     enddo
  enddo
!!$
  
  !Inizio a passare dal tb al real space
  sso=0
  call bielectron(dim2,nso,pf,vecconfig,ssotb,sso)
  call check_hermitian(sso, dim2, bool)
  if(bool)write(*,*)'sso bene'
  !SOO TERM
  allocate(hsootb(3,nso,nso,nso,nso),sootb(nso,nso,nso,nso),soo(dim2,dim2))
  hsootb=0
 
  do a=1,nso
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if(((b+1)/2.eq.(d+1)/2).and.((a+1)/2.ne.(b+1)/2))then
                 vec1=0
                 vec2=0
                 cp=0
                 do k=1,3
                    vec1(k)=dist((a+1)/2,(b+1)/2,k)
                    vec2(k)=ppso(k,a,c)
                 enddo
                 cp=cross_product(vec1,vec2)
                ! write(*,*) cp
!!$               ! 1/|r_aA|^3 term
                 radius = 0.d0
                 do k = 1,3
                    radius = radius + dreal(vec1(k))**2
                 end do
                 radius = (dsqrt(radius))**3
                ! write(*,*) radius
                 cp = cp / radius
                 if(mod(b,2).eq.0)then
                    si=2
                 else
                    si=1
                 endif

                 if(mod(d,2).eq.0)then
                    sj=2
                 else
                    sj=1
                 endif

                 do k=1,3
                    hsootb(k,a,b,c,d)=hsootb(k,a,b,c,d)+cp(k)*spin(si,sj,k)
                   ! write(*,*) hssotb(k,a,b,c,d)
                 enddo
                
              endif
           enddo
        enddo
     enddo
  enddo

!!$
  sootb=0
 
  do k=1,3
     do a=1,nso
        do b=1,nso
           do c=1,nso
              do d=1,nso
                 sootb(a,b,c,d)=sootb(a,b,c,d)+0.5d0*(hsootb(k,a,b,c,d)+dconjg(hsootb(k,c,d,a,b)))
                 temporary(k,a,b,c,d)=temporary(k,a,b,c,d)+0.5d0*(hsootb(k,a,b,c,d)+dconjg(hsootb(k,c,d,a,b)))
              enddo
           enddo
        enddo
     enddo
  enddo

  
  do a=1,nso
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if(zabs(sootb(a,b,c,d)-dconjg(sootb(c,d,a,b))).ge.1d-10)write(*,*) a, b, c, d, zabs(sootb(a,b,c,d)-dconjg(sootb(c,d,a,b)))
              if(zabs(sootb(a,b,c,d)).ge.1d-10)write(77,*) a,b,c,d
           enddo
        enddo
     enddo
  enddo
!!$
  
  !Inizio a passare dal tb al real space
  soo=0

  call  bielectron(dim2,nso,pf,vecconfig,sootb,soo)
  call check_hermitian(soo,dim2,bool)
  if(bool)write(*,*)'soo bene'

  allocate(bi_coup(3,dim2,dim2))
  bi_coup=0
  do k=1,3
     call  bielectron(dim2,nso,pf,vecconfig,temporary(k,:,:,:,:),bi_coup(k,:,:))
  enddo

  
  !============================================================
  allocate(ham2(dim2,dim2))
  !! Adding to hamiltonian
  do n=1,dim2
     ham(n,n)=pot(n)+energy(n)
     do m=1,dim2
        ham(n,m)=ham(n,m)-tt(n,m)+coup(n,m)-soo(n,m)-sso(n,m)
     enddo
  enddo
  soc=0
  do n=1,dim2
     do m=1,dim2
        soc(n,m)=soc(n,m)-sso(n,m)-soo(n,m)+coup(n,m)
     enddo
  enddo
  
  do n=1,dim2
     ham2(n,n)=pot(n)+energy(n)
     do m=1,dim2
        ham2(n,m)=ham2(n,m)-tt(n,m)
     enddo
  enddo

  !HAM2 HUBBARD HAMILTONIAN


  
  call  check_hermitian(ham, dim2, is_hermitian)
  if(is_hermitian) write(*,*) 'HAM HERMITIANA'
  jobz ='V'
  uplo='U'
  lrwork=(1+5*dim2+2*dim2**2)
  liwork=(3+5*dim2)
  lwork=(2*dim2+dim2**2)
  allocate(w(dim2),work(max(1,lwork)),rwork(lrwork),iwork(max(1,liwork)),w2(dim2), eig(dim2,3))
!==============================DIAGONALIZATION=========================
  call zheevd (jobz, uplo, dim2, ham, dim2, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  call eigenvalues(dim2,1d-8,w,state)

  call rotate_cplx_2x2(dim2,sqrot2,sqrot,ham)
  call check_hermitian(sqrot2, dim2, bool)
  if(.not.bool)write(*,*) 'AAAAH'
  
  write(4,*) 'EIGENVALUES'
  do i=1,12
     write(4,*) w(i)-w(1), state(i), dreal(sqrot2(i,i))
  enddo
  eig=0
  do i=1,dim2
     eig(i,3)=w(i)-w(1)
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
  
  write(4,*) 'DIPOLE ALPHA'
  write(4,'(<3>(2x,f10.5))') (muralpha(1,1,k), k=1,3)
  write(4,*) 'DIPOLE BETA'
  write(4,'(<3>(2x,f10.5))') (murbeta(1,1,k), k=1,3)

!!$  write(4,*) 'DIPOLE ALPHA+BETA'
!!$  write(4,'(<3>(2x,f10.5))') (muralpha(1,1,k)+murbeta(1,1,k), k=1,3)
  mu=0
  
  do k=1,3
     call rotate_real_1x2(dim2,mu(:,:,k),dipole(:,k),ham)
  enddo
  
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

  do k=1,3
     call check_hermitian (ppra(k,:,:), dim2, bool)
     if(.not.bool)write(*,*) 'ERROR PPRA K=',k
     call check_hermitian (pprb(k,:,:), dim2, bool)
     if(.not.bool)write(*,*) 'ERROR PPRB K=',k

  enddo
  write(4,*) 'LINEAR MOMENT ALPHA'
  write(4,'(<3>(2x,f15.8))') (dimag(ppra(k,1,1)), k=1,3)
  write(4,*) 'LINEAR MOMENT BETA'
  write(4,'(<3>(2x,f15.8))') (dimag(pprb(k,1,1)), k=1,3)
  
  !=========================CHECK=========================
  w=0
  call zheevd (jobz, uplo, dim2, ham2, dim2, w, work, lwork,rwork,lrwork,iwork,liwork,info)

  do i=1,dim2
     eig(i,1)=w(i)-w(1)
  enddo
 
  w2=w

  coupling=0
  call rotate_cplx_2x2(dim2,coupling,coup,ham2)
  call check_hermitian(coupling, dim2,bool)
  if(.not.bool)write(*,*) 'Problemi coupling'

  socr=0
  call rotate_cplx_2x2(dim2,socr,soc,ham2)
  call check_hermitian(socr,dim2,bool)
  if(.not.bool)write(*,*) 'problemi socr'
  
  sqrot2=0
  call  rotate_cplx_2x2(dim2,sqrot2,sqrot,ham2)
  call check_hermitian(sqrot2,dim2,bool)
  if(.not.bool)write(*,*) 'problemi sqrot2'

  allocate(mcrot(3,dim2,dim2),bcrot(3,dim2,dim2))
  mcrot=0
  bcrot=0

  do k=1,3
     call  rotate_cplx_2x2(dim2,mcrot(k,:,:),mono_coup(k,:,:),ham2)
     call check_hermitian(mcrot(k,:,:),dim2,bool)
     if(.not.bool)write(*,*)'Problemi mcrot'
     call  rotate_cplx_2x2(dim2,bcrot(k,:,:),bi_coup(k,:,:),ham2)
     call check_hermitian(bcrot(k,:,:),dim2,bool)
     if(.not.bool)write(*,*)'Problemi bcrot'
  enddo

  singlet=0
  triplet=0
  quintet=0
  allocate(singlet(dim2), triplet(dim2), quintet(dim2))
  do n=1,dim2
     if(dreal(sqrot2(n,n)).le.1d-15)singlet(n)=1d0
     if(dabs(dreal(sqrot2(n,n))-2d0).le.1d-12)triplet(n)=1d0
     if(dabs(dreal(sqrot2(n,n))-6d0).le.1d-12)quintet(n)=1d0
  enddo

  ham2=0
  do n=1,dim2
     ham2(n,n)=w(n)
     do m=1,dim2
       ham2(n,m)=ham2(n,m)+coupling(n,m)+socr(n,m)
     enddo
  enddo
  call check_hermitian(ham2, dim2,bool)
  if(.not.bool)write(*,*) 'Problemi ham2'
 
  w=0
  call zheevd (jobz, uplo, dim2, ham2, dim2, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  do i=1,dim2
     eig(i,2)=w(i)-w(1)
  enddo
  allocate(sr(dim2,dim2), tr(dim2,dim2), qr(dim2,dim2))
  sr=0
  tr=0
  qr=0
  call rotate_real_1x2(dim2,sr,singlet,ham2)
  call rotate_real_1x2(dim2,tr,triplet,ham2)
  call rotate_real_1x2(dim2,qr,quintet,ham2)
  write(4,*) 'COMPOSITION STATE' 
  do n=1,14
     write(4,'(I10, 4F15.10)')n, sr(n,n)*100, tr(n,n)*100, qr(n,n)*100
  enddo
  do i=1,12
     write(6,'(I10, 5F20.15)') i, w2(i)-w2(1), w(i)-w(1), sr(i,i)*100, tr(i,i)*100, qr(i,i)*100
  enddo
!!$  
  write(4,*) 'EVOLUTION'
  do i=1,12
     write(4,'(I10, 5F20.15)') i, eig(i,1), eig(i,2), eig(i,3)
  enddo

  do i=1,12
     write(99,*) i,(eig(i,2)-eig(i,1))*1d10, (eig(i,3)-eig(i,1))*1d10
  enddo

 
  do i=1,50
     write(100,*) eig(i,1), dsqrt(zabs(coupling(1,i))**2), dsqrt(zabs(socr(1,i))**2)
  enddo
  do i=1,50
     write(102,*) eig(i,1), dsqrt(zabs(coupling(3,i))**2), dsqrt(zabs(socr(3,i))**2)
  enddo

  open(9876,file='post.dat')
  open(9877,file='post1.dat')
  open(9878,file='mcrot_r.dat')
  open(9879,file='mcrot_i.dat')
  open(9883,file='bcrot_r.dat')
  open(9882,file='bcrot_i.dat')
  
  do i=1,dim2
     write(9876,*) dreal(coupling(1,i)),  dimag(coupling(1,i)), dreal(socr(1,i)), dimag(socr(1,i))
     write(9877,*) dreal(coupling(5,i)),  dimag(coupling(5,i)), dreal(socr(5,i)), dimag(socr(5,i))
     !write(9876,*) dsqrt(zabs(coupling(5,i))**2), dsqrt(zabs(socr(5,i))**2)
  enddo

  do i=1,dim2
     write(9878,*) ( dreal(mcrot(k,5,i)), k=1,3)
     write(9879,*) (  dimag(mcrot(k,5,i)), k=1,3)
     write(9883,*)( dreal(bcrot(k,5,i)), k=1,3)
     write(9882,*)(  dimag(bcrot(k,5,i)), k=1,3)
  enddo


  
  !==========================SPIN DENSITY ROTATION=========================
  sdr=0
  do n=1,dim2
     do m=1,dim2
        do l=1,nsiti
           sdr(n,l)=sdr(n,l)+spindensity(m,l)*zabs(ham(m,n))**2
        enddo
     enddo
  enddo
  write(4,*) 'SPIN DENSITY'
  do n=1,dim2
    write(4,'(I2,4F15.8,A2)') n, sdr(n,1), sdr(n,2), sdr(n,3), sdr(n,4), state(n)
 enddo

 !=========================REDFILL=========================
 open(55,file='mu-bin.dat')
 write(55)(mu-bin(:10,:10,k), k=1,3)
 open(66,file='eigen-bin.dat')
 write(66) w(:10)
 open(77, file='singlet-bin.dat')
 write(77)sr(:10)
 open(88, file='triplet-bin.dat')
 write(88)tr(:10)
 

 
end program flash
