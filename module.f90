module module
  implicit none
contains
function cross_product(v1, v2) result(result_vector)
  implicit none
  complex*16, dimension(3), intent(in) :: v1, v2
  complex*16, dimension(3) :: result_vector

  result_vector(1) = v1(2) * v2(3) - v1(3) * v2(2)
  result_vector(2) = v1(3) * v2(1) - v1(1) * v2(3)
  result_vector(3) = v1(1) * v2(2) - v1(2) * v2(1)
end function cross_product

subroutine copy_upper_to_lower(matrix, n)
  implicit none
  integer, intent(in) :: n
  real*8, dimension(n, n), intent(inout) :: matrix

  integer :: i, j

  do i = 1, n
    do j = i+1, n
      matrix(j, i) = matrix(i, j)
    end do
  end do
end subroutine copy_upper_to_lower



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

    !generates a DOUBLE PRECISION one-electron operator in second quantization
    !given a real space basis and one-electron operator in first quantization
    !nso    number of spin-orbitals
    !dim    real-space dimension
    !op_tb  operator in the first quantization (tight binding basis)
    !op     operator in the real space basis
    !basis  1D array containing the integers that describe the real space basis in bit representaion
    !NB: the same spin-orbital ordering used to describe the real space configurations
    !    must be used for the first quantization basis
    subroutine sq_oe_op_real(nso,dim,op_tb,op,basis)
        implicit none
        integer iso,jso,conta,i,j,istate,jstate,step,a
        integer, intent (in) :: dim,nso,basis(dim)
        double precision, intent (in) :: op_tb(nso,nso)
        double precision, intent (out) :: op(dim,dim)
        double precision phase

    op = 0.d0
    do j = 1,dim !col index
        do iso = 0,nso-1 ! creation op index
            do jso = 0,nso-1 ! annihilation op index
                jstate = basis(j)
                if (btest(jstate,jso)) then
                    istate = ibclr(jstate,jso)
                    if (.not.btest(istate,iso)) then

                        istate = ibset(istate,iso)

                        i = binary_search(basis, istate, 1, dim) !row index
                        if (i/=0) then
                            !determine the phase
                            !get direction from iso to jso
                            if (jso>iso) step = -1
                            if (iso>jso) step = 1

                            if (iso==jso) then
                                phase = 1.d0
                                goto 1000
                            end if

                            conta = 0
                            do a = jso+step, iso-step, step
                                if (btest(istate,a)) conta = conta + 1
                            end do

                            if (conta/2*2==conta) then
                                phase = 1.d0
                            else
                                phase = -1.d0
                            end if

                            1000 continue

                            !write(*,*) phase,i,j,(btest(istate,a),a=0,nso-1), 0,0, (btest(jstate,a),a=0,nso-1)
                            op(i,j) = op(i,j) + phase * op_tb(iso+1,jso+1)
                        end if
                    end if
                end if

            end do
        end do
        !write(*,*) ''

    end do
    end subroutine sq_oe_op_real

    !generates a DOUBLE COMPLEX one-electron operator in second quantization
    !given a real space basis and one-electron operator in first quantization
    !nso    number of spin-orbitals
    !dim    real-space dimension
    !op_tb  operator in the first quantization (tight binding basis)
    !op     operator in the real space basis
    !basis  1D array containing the integers that describe the real space basis in bit representaion
    !NB: the same spin-orbital ordering used to describe the real space configurations
    !    must be used for the first quantization basis
    subroutine sq_oe_op_compl(nso,dim,op_tb,op,basis)
    implicit none
    integer iso,jso,conta,i,j,istate,jstate,step,a
    integer, intent (in) :: dim,nso,basis(dim)
    double complex, intent (in) :: op_tb(nso,nso)
    double complex, intent (out) :: op(dim,dim)
    double precision phase

    op = 0.d0
    do j = 1,dim !col index
        do iso = 0,nso-1 ! creation op index
            do jso = 0,nso-1 ! annihilation op index
                jstate = basis(j)
                if (btest(jstate,jso)) then
                    istate = ibclr(jstate,jso)
                    if (.not.btest(istate,iso)) then

                        istate = ibset(istate,iso)

                        i = binary_search(basis, istate, 1, dim) !row index
                        if (i/=0) then
                            !determine the phase
                            !get direction from iso to jso
                            if (jso>iso) step = -1
                            if (iso>jso) step = 1

                            if (iso==jso) then
                                phase = 1.d0
                                goto 1000
                            end if

                            conta = 0
                            do a = jso+step, iso-step, step
                                if (btest(istate,a)) conta = conta + 1
                            end do

                            if (conta/2*2==conta) then
                                phase = 1.d0
                            else
                                phase = -1.d0
                            end if

                            1000 continue

                            !write(*,*) phase,i,j,(btest(istate,a),a=0,nso-1), 0,0, (btest(jstate,a),a=0,nso-1)
                            op(i,j) = op(i,j) + phase * op_tb(iso+1,jso+1)
                        end if
                    end if
                end if

            end do
        end do
        !write(*,*) ''

    end do
    end subroutine sq_oe_op_compl

    subroutine rotate_real_2x2(dim2,coupling,coup,ham)
      implicit none
      integer::i,l,j,k
      integer, intent(in)::dim2
      real*8,intent(out)::coupling(dim2,dim2)
      real*8,intent(in)::coup(dim2,dim2)
      complex*16,intent(in)::ham(dim2,dim2)

      !$omp parallel do default(none)&
      !$omp private(i,l,j,k)&
      !$omp shared(coup,ham,dim2)&
      !$omp reduction(+:coupling)
      do i=1,dim2 !a
         do l=1,dim2 !b
            do j=1,dim2 !alfa 
               do k=1,dim2 !beta
                  coupling(i,l)=coupling(i,l)+dconjg(ham(j,i))*ham(k,l)*coup(j,k)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
    end subroutine rotate_real_2x2



    subroutine rotate_cplx_2x2(dim2,coupling,coup,ham)
      implicit none
      integer::i,l,j,k
      integer, intent(in)::dim2
      complex*16,intent(out)::coupling(dim2,dim2)
      complex*16,intent(in)::coup(dim2,dim2)
      complex*16,intent(in)::ham(dim2,dim2)

      !$omp parallel do default(none)&
      !$omp private(i,l,j,k)&
      !$omp shared(coup,ham,dim2)&
      !$omp reduction(+:coupling)
      do i=1,dim2 !a
         do l=1,dim2 !b
            do j=1,dim2 !alfa 
               do k=1,dim2 !beta
                  coupling(i,l)=coupling(i,l)+dconjg(ham(j,i))*ham(k,l)*coup(j,k)
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
    end subroutine rotate_cplx_2x2

    !ruota una matrice reale non quadrata sullo spazio degli autostati
    !out è la matrice che ottieni
    !in è la matrice che ruoti
    !rot è la matrice complessa con gli autovettori
    !col è il numero di colonne
    subroutine rot_diag(dim2,out,in,col,rot)
      integer::i,j,k
      integer, intent(in)::dim2,col
      complex*16, intent(in)::rot(dim2,dim2)
      real*8, intent(in)::in(dim2,col)
      real*8,intent(out)::out(dim2,col)
      do i=1,dim2 !a
         do j=1,dim2 !alfa
            do k=1,col
               out(i,k)=out(i,k)+dconjg(rot(j,i))*rot(j,i)*in(j,k)
            enddo
         enddo
      enddo
    end subroutine rot_diag

    subroutine copia_matrice_complessa_multi(matrice_origine, matrice_destinazione, righe, colonne, riga_fissa, colonna_fissa, dim3, dim4)
      integer, intent(in) :: righe, colonne, riga_fissa, colonna_fissa, dim3, dim4
      complex*16, intent(in) :: matrice_origine(righe, colonne, dim3, dim4)
      complex*16, intent(out) :: matrice_destinazione(righe, colonne)
      integer :: i, j

      do i = 1, righe
         do j = 1, colonne
            matrice_destinazione(i, j) = matrice_origine(i, j, riga_fissa, colonna_fissa)
         end do
      end do
    end subroutine copia_matrice_complessa_multi

    subroutine copia_matrice_complessa_cont(matrice_origine, matrice_destinazione, righe, colonne, riga_fissa, colonna_fissa, dim3,dim4)
      integer, intent(in) :: righe, colonne, riga_fissa, colonna_fissa, dim3, dim4
      complex*16, intent(in) :: matrice_origine(righe, colonne)
      complex*16, intent(out) :: matrice_destinazione(righe, colonne,dim3,dim4)
      integer :: i, j

      do i = 1, righe
         do j = 1, colonne
            matrice_destinazione(i, j, riga_fissa, colonna_fissa) = matrice_origine(i, j)
         end do
      end do
    end subroutine copia_matrice_complessa_cont


    subroutine check_hermitian(matrix, n, is_hermitian)
      implicit none
      integer, intent(in) :: n
      complex*16, intent(in) :: matrix(n, n)
      logical, intent(out) :: is_hermitian

      integer :: i, j
      complex*16 :: conj_transpose_matrix(n, n)

      ! Calcola la trasposta coniugata della matrice
      do i = 1, n
         do j = 1, n
            conj_transpose_matrix(i, j) = dconjg(matrix(j, i))
         end do
      end do

      ! Verifica se la matrice è Hermitiana
      is_hermitian = all(matrix == conj_transpose_matrix)

    end subroutine check_hermitian

    subroutine ohno(dim2,nsiti,nuclei,vecconfig,u, nz, pot)
      implicit none
      integer, intent(in)::dim2, vecconfig(dim2), nsiti, nz(nsiti)
      integer::sito, i, p, occupazioni(nsiti), a,b, k,n, j
      real*8::PPP, r(nsiti,nsiti), dx, dy, dz
      logical::bool, bool1
      real*8,intent(in):: nuclei(nsiti,3), u(nsiti)

      real*8, intent(out)::pot(dim2)

      pot=0
      do n=1,dim2
         sito=0
         do i=0,2*nsiti-2,2
            sito=(i+2)/2
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
            occupazioni(sito)=a+b
         enddo

         r = 0.0d0

         do i = 1, nsiti
            do j = i + 1, nsiti
               dx = nuclei(i,1) - nuclei(j,1)
               dy = nuclei(i,2) - nuclei(j,2)
               dz = nuclei(i,3) - nuclei(j,3)

               r(i, j) = dsqrt(dx**2 + dy**2 + dz**2)
               r(j, i) = r(i, j)
            enddo
         enddo

         PPP=0.d0
         do i=1,nsiti
            do p=1,nsiti
               if(i.ne.p) PPP=PPP+(14.397)/dsqrt(r(i,p)**2+(28.794/(u(i)+u(p)))**2)*(nz(i)-occupazioni(i))*(nz(p)-occupazioni(p))
               if((i.eq.p).and.(occupazioni(i).eq.2))PPP=PPP+2*u(i)
            enddo
         enddo
         pot(n)=0.5*PPP
      enddo
    end subroutine ohno

    subroutine site_energy_ohno(nso,dim2,esite,u,vecconfig,energy)
      implicit none
      integer::n,i, sito
      real*8, intent(in)::esite(nso/2), u(nso/2)
      integer,intent(in)::dim2, vecconfig(dim2), nso
      real*8, intent(out)::energy(dim2)
      logical::bool, bool1
      
      
      energy(n)=0
      do n=1,dim2
         do i=0, nso-1
            sito=(i+2)/2
            if(btest(vecconfig(n),i))energy(n)=energy(n)+esite(sito)
         enddo
      enddo

      do n=1,dim2
         do i=0, nso-2,2
            sito=(i+2)/2
            bool=btest(vecconfig(n),i)
            bool1=btest(vecconfig(n),i+1)
            if(bool.and.bool1)energy(n)=energy(n)+u(sito)
         enddo
      enddo

    end subroutine site_energy_ohno

      subroutine site_energy(nso,dim2,esite,vecconfig,energy)
      implicit none
      integer::n,i, sito
      real*8, intent(in)::esite(nso/2)
      integer,intent(in)::dim2, vecconfig(dim2), nso
      real*8, intent(out)::energy(dim2)
      logical::bool, bool1
      
      
      energy(n)=0
      do n=1,dim2
         do i=0, nso-1
            sito=(i+2)/2
            if(btest(vecconfig(n),i))energy(n)=energy(n)+esite(sito)
         enddo
      enddo
    end subroutine site_energy
            
            



    
    subroutine momentum_so (nsiti, nso, hop, nuclei, pp_so)
      implicit none
      integer::k, i, j
      real*8, intent(in)::hop(nsiti,nsiti), nuclei(nsiti,3)
      complex*16, intent(out)::pp_so(3,nso,nso)
      complex*16::cplx, p(3,nsiti,nsiti)
      integer, intent(in):: nso, nsiti
      cplx=cmplx(0.d0,1.d0)
      p=0
      pp_so=0
      do k=1,3
         do i=1,nsiti
            do j=1,nsiti
               p(k,i,j)=cplx*hop(i,j)*(nuclei(i,k)-nuclei(j,k))
            enddo
         enddo
      enddo
      
      do K=1,3
         do i=1,nso-2
            pp_so(k,i,i+2)=p(k,(i+1)/2, (i+3)/2)
            pp_so(k,i+2,i)=p(k,(i+3)/2,(i+1)/2)
         enddo
      enddo
    end subroutine momentum_so

    
    subroutine write_matrix (matrix, numfile, righe, colonne)
      integer::i,j
      integer, intent(in):: numfile
      integer, intent(in)::righe, colonne
      real*8,intent(in)::matrix(righe,colonne)
      do i=1,righe
         write(numfile,'(<colonne>(2x,f10.5))') (matrix(i,j), j=1,colonne)
      enddo
    end subroutine write_matrix

    subroutine eigenvalues(dim2,tollerance,w,state)
      
      
      integer, intent(in) :: dim2
      real*8, intent(in) :: w(dim2)
      real*8, intent(in) :: tollerance
      character(len=1), intent(out) :: state(dim2)

      integer :: i, count_result, j

      do i = 1, dim2
         ! Inizializza il contatore
         count_result = 0

         ! Conta gli elementi che soddisfano la condizione
         do j = 1, dim2
            if (abs(w(j) - w(i)) < tollerance) then
               count_result = count_result + 1
            endif
         end do

         ! Determina lo stato in base al risultato del conteggio
         if (count_result == 3) then
            state(i) = 'T'
         elseif (count_result == 5) then
            state(i) = 'Q'
         else
            state(i) = 'S'
         endif
      end do

      
    end subroutine eigenvalues

    subroutine charge(carica, vecconfig, nz, dim2, nso)
      implicit none
      integer, intent(in):: nz(nso/2), dim2, vecconfig(dim2), nso
      real*8, intent(out):: carica(dim2, nso/2)
      integer:: n, i, a, b, sito
      logical::bool, bool1

      carica = 0

      do n = 1, dim2
         do i = 0, nso-1, 2
            bool = btest(vecconfig(n), i)
            bool1 = btest(vecconfig(n), i+1)

            if (bool) then
               a = 1
            else
               a = 0
            endif

            if (bool1) then
               b = 1
            else
               b = 0
            endif
            sito=(i+2)/2
            carica(n, sito) = nz(sito) - (a + b)
         enddo
      enddo
    end subroutine charge

    subroutine dipole_moment(dipole, carica, nuclei, dim2, nsiti)
      implicit none
      integer, intent(in)::dim2, nsiti
      real*8, intent(in):: carica(dim2,nsiti), nuclei(nsiti,3)
      real*8, intent(out):: dipole(dim2,3)
      integer::n,k,j

      dipole = 0.0d0

      do n = 1, dim2
         do k = 1, 3
            do j = 1, nsiti
               dipole(n, k) = dipole(n, k) + carica(n, j) * nuclei(j, k)
            enddo
         enddo
      enddo
    end subroutine dipole_moment
!!$
!!$    subroutine copia_quad(matrice, dimensione)
!!$      complex(16), intent(inout) :: matrice(dimensione,dimensione)
!!$      integer, intent(in) :: dimensione
!!$      integer :: i, j
!!$
!!$      do i = 1, dimensione
!!$         do j = i+1, dimensione
!!$            matrice(j, i) = (matrice(i, j))
!!$         end do
!!$      end do
!!$
!!$    end subroutine copia_quad

    subroutine quadrato_matrice_tridimensionale(matrice, dim2, strati,quad)
      integer,intent(in)::strati, dim2
      complex*16,intent(in)::matrice(strati,dim2,dim2)
      complex*16,intent(out)::quad(strati,dim2,dim2)
      integer::i,j,k,m
      complex*16::temp(dim2,dim2)

      ! Calcola il quadrato di ogni "strato" lungo la dimensione k
      do k = 1, strati
         do i = 1, dim2
            do j = 1, dim2
               temp(i, j) = (0.0, 0.0)
               do m = 1, dim2
                  temp(i, j) = temp(i, j) + matrice(k, m, i) * matrice(k, m, j)
               end do
            end do
         end do

         ! Copia il risultato nella matrice originale
         do i = 1, dim2
            do j = 1, dim2
               quad(k, i, j) = temp(i, j)
            end do
         end do
      end do

    end subroutine quadrato_matrice_tridimensionale
  end module module
