program basis
  implicit none
  integer:: nsiti,i,max,min, check1, check2, nso, count, config, nf, a,b, n, spin, ne
 ! real*8::spin
  character:: array(8)
  logical::bool
  nsiti=4 !siti di disposizione elettroni
  nso=nsiti*2
  ne=2
  open(1,file='basis.dat')
  open(2,file='configurations.dat')
  max=0
  do i=2*nsiti-1,2*nsiti-2,-1
     max=max+2**i
  enddo
  min=0
  do i=0,1
     min=min+2**i
  enddo

!  write(1,*)min,max
!!$  check1=0
!!$  do i=4,7
!!$     check1=check1+2**i
!!$  enddo
!!$  check2=0
!!$  do i=0,3
!!$     check2=check2+2**i
!!$  enddo
!!$


!!!! Starting whit writing of configurations
  nf=0
  do n=min,max
     count=0
     config=0
     a=0
     b=0
     spin=0
     do i=0,nso-1
        bool=btest(n,i)
        if(bool)then
           array(i+1)='1'
           count=count+1
           if(i/2*2.eq.i)then
              a=a+1
           else
              b=b+1
           endif
        else
           array(i+1)='0'          
        endif
     enddo
    
     if(count.eq.ne)then
        nf=nf+1
        spin=(a-b)*0.5d0
        do i=0,nso-1
           if(array(i+1).eq.'1')then
              config=config+2**i
           endif
        enddo
        write(1,*) config
        write(2,*) (array(i),i=nso,1,-1), spin
     endif
  enddo
  write(*,*)'numero di funzioni di base uguale', nf
  write(2,*)'numero di funzioni di base uguale=', nf

endprogram basis
