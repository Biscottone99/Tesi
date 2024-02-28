program breit

  implicit none
  character*1:: jobz,uplo
  integer:: lda, lwork,isito,irwork,liwork,info, lrwork, si, sj,nso,conta, count
  integer,allocatable::iwork(:)
  integer::nsiti,dim2, i, n, a, b, c, d, ab, cd, p, nf,j, k, l, perm, nao, ne
  integer::  temp, m, binarysearch, contasito
  integer,allocatable:: vecconfig(:), occupazioni(:), NZ(:), spar(:,:), search(:)
  real*8,allocatable:: cord(:,:),rwork(:), w(:), dist(:,:,:), hop(:,:), nuclei(:,:), hop2(:,:), spintot(:), carica(:,:), q(:), u(:),esite(:)
  complex*16,allocatable::ham(:,:),work(:), hamsoc(:,:),  soc_a(:,:,:), soc_b(:,:,:),soc_mono(:,:,:), pp(:,:,:), coup(:,:), COUPLING(:,:)
  real*8,allocatable:: umat(:,:), emat(:,:), corrint(:,:), correx(:,:), vmat(:,:,:),double(:), double2(:,:), dsite(:,:), ssite(:), spin3(:), spin2(:), sz(:), gigi(:)
  real*8:: Uc, t, PPP, me, gs, e, e0, pi, cl, radius
  logical:: bool, bool1, bool2, bool3
  real*8::strenght,hbar,hbarev,echarge,emass,dpg, bubu, balu, delta, gamma, gamma2, tollerance
  complex*16::sy(2,2), spin(2,2,3),vec1(3),vec2(3), cp(3)
  complex*16::cplx,pf
  character*1,allocatable::state(:)
  
 
  nsiti=4
  ne=2
  
  Uc=4d0
  t=1.d0
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

  open(1,file='basis.dat')
  open(2,file='geom.dat')
  open(3,file='hamiltonian.dat')
  open(4,file='results.dat')
 ! open(5,file='input.dat')
  open(6,file='check.dat')
  open(7,file='eig.dat')
  open(8,FILE='charges.dat')
  open(9,file='work.dat')
  open(11,file='spin.dat')
  open(10,file='dim2.dat')
  read(10,*) dim2
  close(10)
  nf=dim2

  allocate(vecconfig(dim2), cord(nsiti,3), nuclei(nsiti,3),ham(dim2,dim2), occupazioni(nsiti),dist(nsiti,nsiti,3), soc_a(3,nso,nso), soc_b(3,nso, nso), pp(3, nsiti,nsiti),hop(nsiti,nsiti),hop2(nso,nso),soc_mono(3,nso,nso),hamsoc(nso,nso), spintot(dim2), carica(nsiti,dim2), q(dim2), nz(nsiti), u(nsiti), esite(nsiti), coup(dim2,dim2), COUPLING(DIM2, dim2), spar(dim2,dim2), spin3(dim2), spin2(dim2), state(dim2),search(dim2), sz(dim2),gigi(dim2))

  do i=1,dim2
     read(1,*) vecconfig(i)
     read(11,*) sz(i)
  enddo
  
  do i=1,dim2
     search(i)=vecconfig(i)
  enddo

  
  u=uc
  esite=0
!!$  U(2)=uc
!!$  U(3)=u(2)
!!$  u(1)=-delta
!!$  u(4)=+delta
!!$

  esite(1)=-8d0
  esite(4)=-2d0
  
!!$  do i=2,4
!!$     u(i)=0
!!$  enddo

  nz(4)=0
  nz(1)=2
  

  do i=1,nsiti
     read(2,*) cord(i,1), cord(i,2), cord(i,3)
  enddo
  rewind(2)
  do i=1,nsiti
     read(2,*) nuclei(i,1), nuclei(i,2), nuclei(i,3)
  enddo


  do i=1,nsiti
     do p=1,nsiti
        do  k=1,3
           dist(i,p,3)=(cord(i,1)-cord(p,1))
        enddo
     enddo
  enddo

  do n=1,dim2
     contasito=0
     do i=0,2*nsiti-2,2
        contasito=contasito+1
        bool=btest(vecconfig(n),i)
        bool1=btest(vecconfig(n),i+1)
        if(bool)then
           a=1
        else
           a=0
        endif
        if(bool1)then
           b=1
        else
           b=0
        endif
        occupazioni(contasito)=a+b
     enddo


     PPP=0.d0
     do i=1,nsiti
        do p=1,nsiti
           PPP=PPP+(14.397)/((28.794/(u(i)-u(p)))**2+((cord(i,1)-cord(p,1))**2+(cord(i,2)-cord(p,2))**2+(cord(i,3)-cord(p,3))))**0.5*(nz(i)-occupazioni(i))*(nz(p)-occupazioni(p))
        enddo
     enddo
     ham(n,n)=0.5*PPP
  enddo

  do n=1,dim2
     do i=0,nso-1,2
        bool=btest(vecconfig(n),i)
        bool1=btest(vecconfig(n),i+1)
        if(bool)ham(n,n)=ham(n,n)+esite((i+2)/2)
        if(bool1)ham(n,n)=ham(n,n)+esite((i+2)/2)
        if(bool.and.bool1)ham(n,n)=ham(n,n)+u((i+2)/2)
     enddo
  enddo
        
        
  !=========================SUPPORTING=========================
  hop2=0
  do i=1,nso-2
     hop2(i,i+2)=t
  enddo
!============================================================

  do n=1,dim2
     do i=0,nso-3
        bool=btest(vecconfig(n),i)
        bool1=btest(vecconfig(n),i+2)
        if(bool.and.(.not.bool1))then
           temp=ibclr(vecconfig(n),i)
           temp=ibset(temp,i+2)
           m=binary_search(vecconfig,temp,1,dim2)
           if(m.ne.0)then
              perm=0
              do a=i+1,nso-1
                 bool2=btest(vecconfig(n),a)
                 if(bool2)perm=perm+1
              enddo
              
              if(i+2.ne.nso-1)then
                 do a=i+3,nso-1
                    bool2=btest(vecconfig(n),a)
                    if(bool2)perm=perm+1
                 enddo
              endif
              if(perm/2*2.eq.perm)then
                 ham(n,m)=-hop2(i+1,i+3)
                 write(89,*) n, m, i, i+2,dreal( ham(n,m))
              endif
              if(perm/2*2.ne.perm)then
                 ham(n,m)=+hop2(i+1,i+3)
                 write(89,*) n, m, i, i+2, dreal(ham(n,m))
              endif
           endif
        endif
     enddo
  enddo

  !=========================SUPPORTING=========================
    
  hop=0
  do i=1,nsiti-1
     hop(i,i+1)=t
  enddo

  do i=1,nsiti
     do j=1,nsiti
        hop(j,i)=hop(i,j)
     enddo
  enddo

   do k=1,3
     do i=1,nsiti
        do j=1,nsiti
           pp(k,i,j)=cplx*hop(i,j)*(nuclei(i,k)-nuclei(j,k))
        enddo
     enddo
  enddo

  do k=1,3
     write(6,*) 'momentum dimension', k
     do i=1,nsiti
        write(6,'(<nsiti>(2x,f10.5))') (dimag(pp(k,i,j)), j=1,nsiti)
     enddo
  enddo
  
  !pauli matrix
  spin(1,2,1)=1d0
  spin(2,1,1)=1d0

  spin(1,2,2)=cplx
  spin(2,1,2)=-cplx
  
  spin(1,1,3)=1d0
  spin(2,2,3)=-1d0

  do i=1,nso
     do j=1,nso
        do isito = 1,nsiti
           if (isito.ne.(i+1)/2) then 
              ! r x p cross product
              vec1 = (0.d0, 0.d0)
              vec2 = (0.d0, 0.d0)
              do k = 1,3
                 vec1(k) = nuclei((i+1)/2,k) - nuclei(isito,k) !position
                 vec2(k) = pp(k,(i+1)/2,(j+1)/2) !momentum
              end do
             ! write(*,*) dimag(vec2)
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
              do k = 1,3
                 vec1(k) = nuclei((j+1)/2,k) - nuclei(isito,k) !position
                 vec2(k) = pp(k,(i+1)/2,(j+1)/2) !momentum
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
                 soc_b(k,i,j) = soc_b(k,i,j) + cp(k)  * spin(si,sj,k)
                 !write(*,*) cp(k)!soc_b(k,i,j)
              end do
           end if
        end do

     end do
  end do

 
  soc_mono = (0.d0, 0.d0)
  do i = 1,nso
     do j = 1,nso
        do k = 1,3
           soc_mono(k,i,j)    = 0.5d0 * (soc_a(k,i,j) + soc_b(k,i,j))
        end do
     end do
  end do
  nao=nso
 
  hamsoc=0d0
  do i=1,nso
     do j=1,nso
        do k=1,3
           hamsoc(i,j)=hamsoc(i,j)+soc_mono(k,i,j)
        enddo
     enddo
  enddo
!!$  write(6,*) 'REAL HAMSOC'
!!$  do i=1,nso
!!$     write(6,'(<dim2>(2x,f10.5))') (dreal(hamsoc(i,j)),j = 1,nso)
!!$  enddo
!!$
!!$  write(6,*) 'IMAG HAMSOC'
!!$  do i=1,nso
!!$     write(6,'(<dim2>(2x,f10.5))') (dimag(hamsoc(i,j)),j = 1,nso)
!!$  enddo

  write(6,*) 'HAMSOC'
  do i=1,nso
     do j=1,nso
        if(zabs(hamsoc(i,j)).ge.1d-10)write(6,*) i, j, hamsoc(i,j)
     enddo
  enddo
  
  write(6,*) 'REAL HAM'
  do i=1,dim2
     write(6,'(<dim2>(2x,f10.5))') (dreal(ham(i,j)),j = 1,dim2)
  enddo

  write(6,*) 'IMAG HAM'
  do i=1,dim2
     write(6,'(<dim2>(2x,f10.5))') (dimag(ham(i,j)),j = 1,dim2)
  enddo


  !==========================SOC=========================
  do n=1,dim2
     do i=0,nso-3
        bool=btest(vecconfig(n),i)
        if(bool)then
           if(mod(i,2).eq.0)then
              do j=i+2,nso-1
                 bool1=btest(vecconfig(n),j)
                 if(.not.bool1)then
                    temp=ibclr(vecconfig(n),i)
                    temp=ibset(temp,j)

                    m=binary_search(vecconfig,temp,1,dim2)
                    if(m.ne.0)then
                       perm=0
                       do a=i+1,nso-1 !! calcolo elettroni tra sito di partenza e fine
                          bool2=btest(vecconfig(n),a)
                          if(bool2)perm=perm+1
                       enddo
                       if(j.ne.7)then
                          do a=j+1,nso-1 !!calcolo elettroni tra sito di arrivo e fine
                             bool3=btest(vecconfig(n),a)
                             if(bool3)perm=perm+1
                          enddo
                       endif

                       if(mod(perm, 2) .ne. 0)then
                          coup(n,m)=coup(n,m)-pf*hamsoc(i+1,j+1)
                          write(88,*) n,m, i,j,i+1, j+1, coup(n,m)
                       else
                          coup(n,m)=coup(n,m)+pf*hamsoc(i+1,j+1)
                          write(88,*) n,m, i,j, i+1, j+1, coup(n,m)
                       endif
                    endif
                 endif
              enddo
           else
              do j=i+1,nso-1
                 bool1=btest(vecconfig(n),j)
                 if(.not.bool1)then
                    if(mod(i, 2) .ne. 0)then
                       temp=ibclr(vecconfig(n),i)
                       temp=ibset(temp,j)
                       m=binary_search(vecconfig,temp,1,dim2)
                       if(m.ne.0)then
                          perm=0
                          do a=i+1,nso-1 !! calcolo elettroni tra sito di partenza e fine
                             bool2=btest(vecconfig(n),a)
                             if(bool2)perm=perm+1
                          enddo
                          if(j.ne.7)then
                             do a=j+1,nso-1 !!calcolo elettroni tra sito di arrivo e fine
                                bool3=btest(vecconfig(n),a)
                                if(bool3)perm=perm+1
                             enddo
                          endif

                          if(mod(perm, 2) .ne. 0)then
                             coup(n,m)=coup(n,m)-pf*hamsoc(i+1,j+1)
                             write(88,*) n,m, i,j,i+1, j+1, coup(n,m)
                          else
                             coup(n,m)=coup(n,m)+pf*hamsoc(i+1,j+1)
                             write(88,*) n,m, i,j, i+1, j+1, coup(n,m)
                          endif
                       endif
                    endif
                 endif
              enddo
           endif
        endif
     enddo
  enddo
            

write(3,*) 'REAL HAMILTONIAN'
do i = 1,dim2
   write(3,'(<dim2>(2x,f10.5))') (dreal(ham(i,j)),j = 1,dim2)
enddo
write(3,*) 'IMAGINARY HAMILTONIAN'
do i=1,dim2
   write(3,'(<dim2>(2x,f10.5))') (dimag(ham(i,j)),j = 1,dim2)
end do
!!$
!!$do i=1,dim2
!!$   do j=1,dim2
!!$      if(dabs(dimag(ham(i,j))).ge.1d-17)then
!!$         write(*,*) i, j, dimag(ham(i,j))
!!$      endif
!!$   enddo
!!$enddo





!=========================DIAGONALIZATION=========================
jobz ='V'
uplo='U'
lrwork=(1+5*nf+2*nf**2)
liwork=(3+5*nf)
lwork=(2*nf+nf**2)
allocate(w(nf),work(max(1,lwork)),rwork(lrwork),iwork(max(1,liwork)))

call zheevd (jobz, uplo, nf, ham, nf, w, work, lwork,rwork,lrwork,iwork,liwork,info)

write(4,*) 'EIGENVALUES'
do i=1,nf
   write(7,*) w(i)-w(1)
   write(4,*) i, w(i)
enddo

write(4,*) 'REAL EIGENVECTORS'
do i=1,nf
   write(4,'(<nf>(2x,f10.5))')(dreal(ham(i,j)),j=1,nf)
enddo

write(4,*) 'IMAGINARY EIGENVECTORS'
do i=1,nf
   write(4,'(<nf>(2x,f10.5))')(dimag(ham(i,j)),j=1,nf)
enddo
!=========================COUPLING=========================
do i=1,dim2 !a
   do l=1,dim2 !b
      do j=1,dim2 !alfa 
         do k=1,dim2 !beta
            coupling(i,l)=coupling(i,l)+dconjg(ham(j,i))*ham(k,l)*coup(j,k)
         enddo
      enddo
   enddo
enddo

write(77,*) 'IMAG COUP'
do i=1,dim2
   write(77,'(<dim2>(2x,f10.5))') (dimag(coup(i,j)),j = 1,dim2)
enddo


write(77,*) 'real COUP'
do i=1,dim2
   write(77,'(<dim2>(2x,f10.5))') (dreal(coup(i,j)),j = 1,dim2)
enddo


!=========================CHARGES=========================

do n=1,dim2
   do i=0,nso-1,2
      do m=1,dim2
         bool=btest(vecconfig(m),i)
         bool1=btest(vecconfig(m),i+1)
         if(bool)then
            a=1
         else
            a=0
         endif
         if(bool1)then
            b=1
         else
            b=0
         endif
         carica((i+2)/2,n)=carica((i+2)/2,n)+(nz((i+2)/2)-a-b)*(zabs(ham(m,n)))**2
      enddo
   enddo
enddo

do n=1,dim2
   do i=1,nsiti
      q(n)=q(n)+carica(i,n)
   enddo
enddo


write(4,*) 'CARICA'

do i=1,nsiti
   write(4,'(<dim2>(2x,f10.5))') (carica(i,j), j=1,dim2 )
   
enddo

do i=1,nsiti
   write(8,'(<2>(2x,f10.5))') (carica(i,j), j=1,2)
enddo

do i=1,2
   do j=1,nsiti
      write(9,*) carica(j,i)
   enddo
enddo

!==========================Coupling================================
spin3=0
do n=1,dim2
   a=0
   b=0
   do i=0,nso-1
      bool=btest(vecconfig(n),i)
      if(bool)then
         if(i/2*2.eq.i)then
            a=a+1
         else
            b=b+1
         endif
      else
         spin3(n)=(a-b)*0.5d0
      endif
   enddo
enddo
spin2=0
do i=1,dim2 !a
   do j=1,dim2 !alfa 
      spin2(i)=spin2(i)+dconjg(ham(j,i))*ham(j,i)*spin3(j)
   enddo
enddo

tollerance=1d-8
do i = 1, dim2
   if (count(abs(w - w(i)) < tollerance) == 3) then
      state(i) = 'T'
   else if (count(abs(w - w(i)) < tollerance) == 2) then
      state(i) = 'D'
   elseif (count(abs(w - w(i)) < tollerance) == 5) then
      state(i) = 'Q'
   else
      state(i) = 'S'
   endif
end do

do i=1,dim2
   write(77,*) w(i), state(i)
enddo

spar=0
do i=1,dim2
   do j=1,dim2
      if(zabs(coupling(i,j)).ge.1d-5)then
         spar(i,j)=j
        
      endif
   enddo
enddo

do i=1,6
   write(4,*) 'Configuration=', i, state(i), gigi(i)
   do j=1,dim2
      if(spar(i,j).ne.0)then
         write(4,*) j, state(j), gigi(j), zabs(coupling(i,j))
      endif
   enddo
enddo
         
!============================================================



contains
function cross_product(v1, v2) result(result_vector)
  implicit none
  complex*16, dimension(3), intent(in) :: v1, v2
  complex*16, dimension(3) :: result_vector

  result_vector(1) = v1(2) * v2(3) - v1(3) * v2(2)
  result_vector(2) = v1(3) * v2(1) - v1(1) * v2(3)
  result_vector(3) = v1(1) * v2(2) - v1(2) * v2(1)
end function cross_product




recursive function binary_search(arr, target, low, high) result(index)
    integer, intent(in) :: arr(:), target, low, high
    integer :: index, mid

    if (low > high) then
      index = 0  ! Target not found
    else
      mid = (low + high) / 2

      if (arr(mid) == target) then
        index = mid  ! Target found
      else if (arr(mid) < target) then
        index = binary_search(arr, target, mid + 1, high)
      else
        index = binary_search(arr, target, low, mid - 1)
      end if
    end if
   end function binary_search



   
end program

!__________________________________________________
