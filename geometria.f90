program geometria
  implicit none
  real*8:: l, alpha
  integer:: i
  real*8:: x, y,z, pi
  real*8,allocatable::cord(:,:)
  l=1.50
  pi=dacos(-1.d0)
  alpha=pi/4
  write(*,*) dsin(alpha), dcos(alpha)
  open(1,file='geometria.dat')
  open(2,file='geom.dat')
  allocate(cord(4,3))
  do i=1,4
     if(i==1)then
        cord(i,1)=l*dcos(alpha)
        cord(i,2)=0.d0
        cord(i,3)=0.d0
     endif
     if(i==2)then
        cord(i,1)=0.d0
        cord(i,2)=l*dcos(alpha)
        cord(i,3)=l*dsin(alpha)
     endif
     if(i==3)then
        cord(i,1)=-l*dcos(alpha)
        cord(i,2)=0.d0
        cord(i,3)=2*l*dsin(alpha)
     endif
     if(i==4)then
        cord(i,1)=0.d0
        cord(i,2)=-l*dsin(alpha)
        cord(i,3)=3*l*dsin(alpha)
     endif
     write(1,*) i, cord(i,1), cord(i,2), cord(i,3)
     write(2,*) cord(i,1), cord(i,2), cord(i,3)
  enddo
  
  close(1)
end program geometria
