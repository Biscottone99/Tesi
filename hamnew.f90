program hoppingtotal

  implicit none
  character*1:: jobz,uplo
  integer:: lda, lwork,isito,irwork,liwork,info, lrwork, si, sj,nso,conta
  integer,allocatable::iwork(:)
  integer::nsiti,dim2, i, n, a, b, c, d, ab, cd, p, nf,j, k, l, perm, nao, ne
  integer::  temp, m, binarysearch, contasito
  integer,allocatable:: vecconfig(:), occupazioni(:)
  real*8,allocatable:: cord(:,:),rwork(:), w(:), dist(:,:,:), hop(:,:), nuclei(:,:), hop2(:,:)
  complex*16,allocatable::ham(:,:),work(:), hamsoc(:,:),  soc_a(:,:,:), soc_b(:,:,:),soc_mono(:,:,:), pp(:,:,:)
  real*8,allocatable:: umat(:,:), emat(:,:), corrint(:,:), correx(:,:), vmat(:,:,:),double(:), double2(:,:), dsite(:,:), ssite(:)
  real*8:: Uc, t, PPP, me, gs, e, e0, pi, cl, radius
  logical:: bool, bool1, bool2, bool3
  real*8::strenght,hbar,hbarev,echarge,emass,dpg, bubu, balu
  complex*16::sy(2,2), spin(2,2,3),vec1(3),vec2(3), cp(3)
  complex*16::cplx,pf
  
  external binarysearch
  nsiti=4
  ne=2
  
  Uc=11.26d0
  dim2=28
  nf=dim2
  t=2.5d0
  nso=nsiti*2
  me=9.1093837015d-31
  gs=2.00231930436256
  e=1.602176634d-19
  e0=8.8541878128d-12
  pi=dacos(-1.d0)
  write(*,*) pi
  cl=299792458 
  cplx=cmplx(0.d0, 1.d0)


  pf=((gs*e**2)/(8*pi*e0*me*cl**2))*10d10

  open(1,file='basis.dat')
  open(2,file='geom.dat')
  open(3,file='hamiltonian.dat')
  open(4,file='results.dat')
  open(6,file='check.dat')
  open(7,file='eig.dat')
  

  allocate(vecconfig(dim2), cord(nsiti,3), nuclei(nsiti,3),ham(dim2,dim2), occupazioni(nsiti),dist(nsiti,nsiti,3), soc_a(3,nso,nso), soc_b(3,nso, nso), pp(3, nsiti,nsiti),hop(nsiti,nsiti),hop2(nso,nso),soc_mono(3,nso,nso),hamsoc(nso,nso))

  do i=1,dim2
     read(1,*) vecconfig(i)
  enddo

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
           PPP=PPP+(14.397)/((28.794/(2*uc))**2+((cord(i,1)-cord(p,1))**2+(cord(i,2)-cord(p,2))**2+(cord(i,3)-cord(p,3))))**0.5*(1-occupazioni(i))*(1-occupazioni(p))
        enddo
     enddo
     ham(n,n)=0.5*PPP
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
           m=binarysearch(1,dim2,vecconfig,temp)
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
              if(perm/2*2.eq.perm)ham(n,m)=hop2(i+1,i+3)
              if(perm/2*2.ne.perm)ham(n,m)=-hop2(i+1,i+3)

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
  write(6,*) 'REAL HAMSOC'
  do i=1,nso
     write(6,'(<dim2>(2x,f10.5))') (dreal(hamsoc(i,j)),j = 1,nso)
  enddo

  write(6,*) 'IMAG HAMSOC'
  do i=1,nso
     write(6,'(<dim2>(2x,f10.5))') (dimag(hamsoc(i,j)),j = 1,nso)
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
           if(2*i/2.eq.i)then
              do j=i+2,i+3
                 bool1=btest(vecconfig(n),j)
                 if(.not.bool1)then
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

                    temp=ibclr(vecconfig(n),i)
                    temp=ibset(temp,j)
  
                    m=binarysearch(1,dim2,vecconfig,temp)
                              
                 if(2*perm/2.ne.perm)then
                    ham(n,m)=ham(n,m)-pf*hamsoc(i+1,j+1)
                 else
                    ham(n,m)=ham(n,m)+pf*hamsoc(i+1,j+1)
                 endif
                 if(hamsoc(i+1,j+1).ne.0d0)write(77,*) n, m, i+1,j+1
              endif
           enddo
        else
           do j=i+1,i+2
              bool1=btest(vecconfig(n),j)
              if(.not.bool1)then
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
                 temp=ibclr(vecconfig(n),i)
                 temp=ibset(temp,j)

                 m=binarysearch(1,dim2,vecconfig,temp)

                 if(2*perm/2.ne.perm)then
                    ham(n,m)=ham(n,m)-pf*hamsoc(i+1,j+1)
                 else
                    ham(n,m)=ham(n,m)+pf*hamsoc(i+1,j+1)
                 endif
                 if(hamsoc(i+1,j+1).ne.0d0)write(77,*) n, m, i+1,j+1
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







contains
function cross_product(v1, v2) result(result_vector)
  implicit none
  complex*16, dimension(3), intent(in) :: v1, v2
  complex*16, dimension(3) :: result_vector

  result_vector(1) = v1(2) * v2(3) - v1(3) * v2(2)
  result_vector(2) = v1(3) * v2(1) - v1(1) * v2(3)
  result_vector(3) = v1(1) * v2(2) - v1(2) * v2(1)
end function cross_product

end program

!__________________________________________________
integer function binarysearch(i, length, array, val)
  ! Given an array and a value, returns the index of the element that
  ! is closest to, but less than, the given value.
  ! Uses a binary search algorithm.
  ! "delta" is the tolerance used to determine if two values are equal
  ! if ( abs(x1 - x2) <= delta) then
  ! assume x1 = x2
  ! endif

implicit none

integer, intent(in) :: length, i
integer, dimension(length), intent(in) :: array
integer, intent(in) :: val

!integer :: binarysearch

integer :: left, middle, right

left = i
right = length
binarysearch=0


if (val.lt.array(left) .or. val.gt.array(right)) go to 10

do

   if (left .gt. right) then
      exit
      !write(*,*) 'ERRORE!!!'
   endif

   !divisione=((left+right) / 2.0)
   !middle = jnint(divisione)
   middle=(left+right)/2

   if ( array(middle) .eq. val ) then
      binarySearch = middle
      return
   else 
      if (array(middle) .gt. val) then
         right = middle - 1
      else
         left = middle + 1
      end if
   end if
end do

binarysearch = right
10 continue
end function binarysearch



