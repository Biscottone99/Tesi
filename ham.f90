program hamiltonian
  implicit none
  integer:: a, b, c, i, j, k, nsiti,nso, count, ne, n, lda, lwork, info, nf
  integer,allocatable:: config(:)
  real*8::ax, ay, az, temp, uc,energy, t, temp2
  real*8, allocatable:: cord(:,:),rad(:,:,:), pot(:,:), charge(:,:), opera(:,:), ham(:,:), w(:), work(:),eig(:,:)
  logical::bool, bool1
  character*1:: jobz, uplo

  nsiti=4
  ne=2
  uc=11.26d0
  nso=nsiti*2
  nf=28
  t=2.50d0
  allocate (cord(nsiti,3),config(nf), charge(nf,nsiti), pot(nsiti,nsiti), rad(3,nsiti,nsiti),opera(nsiti,nsiti), ham(nsiti,nsiti), eig(nsiti,nsiti) )


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



!!!! PPP part
  ham=0

  !!! Interelectronic interaction
  do i=1,4
     ham(i,i)=uc
  enddo
!!!! potential term - Ohno parametrization 
  do i=1,4
     temp=0
     do j=1,4
        temp=temp+(14.397)/((28.794/(2*uc)**2+((cord(i,1)-cord(j,1))**2+(cord(i,2)-cord(j,2))**2)+(cord(i,3)-cord(j,3))**2))**0.5
     enddo
     ham(i,i)=ham(i,i)+0.5*temp
  enddo


  !!! off-diagonal part
  do i=1,nsiti-1
     ham(i,i+1)=-t
  enddo
  do i=2,4
     ham(i,i-1)=-t
  enddo


  do i=1,nsiti
     write(5,'(<nsiti>(2x,f10.5))')(ham(i,j),j=1,nsiti)
  enddo
  close(5)


  !!! supporting operator - distances between different sites for each dimension (x, y, z)
  do i=1,4
     do j=1,4
        do k=1,3
           rad(k,i,j)=cord(i,k)-cord(j,k)
        enddo
     enddo
  enddo
  do k=1,3
     write(4,*) 'DIMENSIONE K=',k
     do i=1,nsiti
        write(4,'(<nsiti>(2x,f10.5))')(rad(k,i,j),j=1,nsiti)
     enddo
  enddo



!!!! start writing SOC part

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
           write(5,*) i, i+1, j, temp, temp2
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
           write(6,*) i, i+1, j, temp, temp2
        endif
     enddo
  enddo


  write(4,*) 'MATRIX OF CONFIGURATION', 1
  do i=1,nsiti
     write(4,'(<nsiti>(2x,f10.5))')(opera(i,j),j=1,nsiti)
  enddo
  close(4)

!!! DIAGONALIZATION
  
  call jacobi (ham,eig,1d-10,nsiti)
  write(3,*) 'EIGENVALUES'
  do i=1,nsiti
     write(3,*) i, ham(i,i)
  enddo
  energy=0.d0
  do i=1,2
     energy=energy+2*ham(i,i)
  enddo
  energy=energy

  write(3,*) 'EIGENVECTORS'
  do i=1,nsiti
     write(3,'(<nsiti>(2x,f10.5))')(eig(i,j),j=1,nsiti)
  enddo


endprogram hamiltonian

subroutine Jacobi(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer i, j, k, n
double precision a(n,n),x(n,n)
double precision abserr, b2, bar
double precision beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do
return
end subroutine Jacobi
