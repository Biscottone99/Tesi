program ham
  integer:: a, b, c, i, j, k, nsiti,nso, count, ne
  integer,allocatable:: config(:)
  real*8::ax, ay, az, temp
  real*8, allocatable:: cord(:,:),rad(:,:,:), pot(:,:), charge(:,:), opera(:,:)
  logical::bool, bool1
    
  nsiti=4
  ne=2
  
  nso=nsiti*2
  nf=28
  allocate (cord(nsiti,3),config(nf), charge(nf,nsiti), pot(nsiti,nsiti), rad(3,nsiti,nsiti),opera(nsiti,nsiti) )


  open(1,file='geom.dat')
  do i=1,4
     read(1,*) cord(i,1), cord(i,2), cord(i,3)
  enddo
  close(1)
  open(2,file='configurations.dat')
  do i=1, nf
     read(2,*) config(i)
  enddo
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


  open(5,file='check.dat')
  do i=1,nsiti-1
   !  opera(i,i+1)=0
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
   !  opera(i,i+1)=0
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
           write(5,*) i, i+1, j, temp, temp2
        endif
     enddo
  enddo
  
  open(4,file='operatore.dat')
  write(4,*) 'MATRIX OF CONFIGURATION', 1
  do i=1,nsiti
     write(4,'(<nsiti>(2x,f10.5))')(opera(i,j),j=1,nsiti)
  enddo
close(4)
endprogram ham
