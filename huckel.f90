      program main

      implicit none

      integer                                   :: i, n, j, k
      real*8                                    :: beta, beta_plus, beta_minus
      real*8, dimension(:, :), allocatable      :: H, Eigenvectors
      real*8, dimension(:), allocatable         :: Eigenvalues
      integer                                   :: n_cells, atom_cell, n_pi, nvect
      double precision, allocatable             :: alpha(:)
      integer,allocatable                       :: vect(:)
      character(len=20)                         :: dimerization, system_type
      

      open(1,file="huckel.in")
      open(2,file="huckel.out")
      open(3,file="eigenvalues.dat")
      open(4,file="eigenvectors.dat")

      read(1,*)
      read(1,*)
      read(1,*) system_type                                    ! Linear or cyclic
      read(1,*) dimerization                                   ! Type of dimerization 
      read(1,*) n_cells                                        ! Number of unit cells
      read(1,*) atom_cell                                      ! Number of atoms for unit cell
      read(1,*) n_pi                                           ! Number of pi electrons

      allocate(alpha(atom_cell))

      read(1,*) (alpha(i), i=1,atom_cell)
      read(1,*) beta
      read(1,*) beta_plus
      read(1,*) beta_minus

      n = n_cells*atom_cell

      allocate(H(n,n), Eigenvectors(n,n), Eigenvalues(n))
      read(1,*) nvect                                          ! Number of eigenvectors printed in the output
      allocate(vect(nvect))
      read(1,*) (vect(i), i =1, nvect)                         ! Which eigenvectors are printed in the output

      close(1)

      H = 0.d0                                                 ! Initialization of the matrix

! Build the Huckel matrix -> 4 cases: no dimerization, atom dimerization, bond dimerization, atom+bond dimerization

! No dimerization, linear or cyclic system

      if (dimerization=='none') then                           ! No dimerization
         do i=1,n
            H(i,i)=alpha(1)
            H(i,i+1)=beta
            H(i+1,i)=beta
         end do
      end if

      if (system_type=='cyclic') then
         H(1,n) = beta
         H(n,1) = beta
      else
         H(1,n) = 0.d0
         H(1,n) = 0.d0
      end if


! Introduce dimerization (bond or atom)

      if (dimerization=='bond') then                       ! bond dimerization
         if (system_type=='cyclic') then                   
            if (mod(n,2)==0) then
               H(1,n) = beta_minus
               H(n,1) = beta_minus
            else
               write(*,'(a60)') "The number of atoms must be even if the system is cyclic"
            end if
         end if
         do i = 1, n
            H(i, i) = alpha(1)
         end do
         do i=1, n-1
            if (mod(i,2) /= 0) then
               H(i, i+1) = beta_plus
               H(i+1, i) = beta_plus
            else
               H(i, i+1) = beta_minus
               H(i+1, i) = beta_minus
            end if
         end do
      elseif (dimerization=='atom') then                   ! atom dimerizationm
         do i = 1, n
            k = mod(i-1,atom_cell) + 1
            H(i,i) = alpha(k)
            H(i,i+1) = beta
            H(i+1,i) = beta
         end do
         if (system_type=='cyclic') then
            H(1,n) = beta
            H(1,n) = beta
        end if
      end if
      

! Atom and bond dimerization

      if (dimerization=='atom+bond') then
         if (system_type=='cyclic') then
            if (mod(n,2)==0) then
               H(1,n) = beta_minus
               H(n,1) = beta_minus
            else
               write(*,'(a60)') "The number of atoms must be even if the system is cyclic"
            end if
         end if
         do i = 1, n
            k = mod(i-1,atom_cell) + 1
            H(i,i) = alpha(k)
            if (mod(i,2) == 0) then
               H(i, i+1) = beta_minus
               H(i+1, i) = beta_minus
            else
               H(i, i+1) = beta_plus
               H(i+1, i) = beta_plus
            end if
         end do
      end if


! Diagonalize the matrix

      call diagonalization(H,Eigenvalues,Eigenvectors,N)

! Write the results in the output file      

      write(2,'(a60)') "------------------------------------ INPUT -------------------------------------"
      write(2,'(a30,I3)') "Number of atoms: ", n
      write(2,'(a30,A)') "Type of system: ", system_type
      write(2,'(a30,A)') "Type of dimerization: ", dimerization
      write(2,'(a30,I3)') "Number of unit cells: ", n_cells
      write(2,'(a30,I3)') "Number of atoms per unit cell: ", atom_cell
      write(2,'(a30,I3)') "Number of pi electrons: ", n_pi
      if (dimerization == 'atom' .or. dimerization == 'atom+bond') then
         write(2,*) "Alpha values: ", alpha
      else
         write(2,'(a30,F10.6)') "Alpha value: ", alpha(1)
      end if
      if (dimerization == 'bond' .or. dimerization == 'atom+bond') then
          write(2,'(a15,2f10.6)') "Beta values: ", beta_plus, beta_minus
      else
         write(2, '(a15, F10.6)') "Beta value:", beta
      end if   
      write(2,*) ""
      write(2,'(a60)') "------------------------------------ HUCKEL HAMILTONIAN ------------------------------------"
      write(2,*) ""
      do i=1,n
         write(2,*) H(i,:)                                                                        
      end do
      write(2,*) ""
      write(2,'(a60)') "---------------------------------------- EIGENVALUES ----------------------------------------"
      write(2,*) ""
      do i=1, n
         write(2,"(I3,1X,30(F9.6,1X))") i, ((real(i)-1)/(real(n)-1)), Eigenvalues(i)
      end do
      write(2,*) ""
      
      write(2,'(a60)') "--------------------------------------- HOMO-LUMO gap --------------------------------------"
      write(2,*) "" 
      if (mod(n,2)==0) then
         write(2,*) Eigenvalues((n/2)+1)-Eigenvalues(n/2), 'eV'
      else if (mod(n,2)/=0) then
         write(2,*) Eigenvalues(((n+1)/2)+1)-Eigenvalues((n+1)/2), 'eV'
      end if
      write(2,*) ""
      write(2,'(a60)') "--------------------------------------- EIGENVECTORS ---------------------------------------"
      write(2,*) ""
      write(2, "(5X,30(I3,8X))") (vect(i), i=1,nvect) 

      do i = 1, n
         write(2, "(I3,1X,30(F9.6,1X))") i, (Eigenvectors(i, vect(j)), j=1,nvect)
      end do

! Write the eigenvalues and the eigenvectors in different files in order to plot them

      
      write(3,'(a15)') "# Eigenvalues"
      do i=1, n
         write(3,"(I3,1X,30(F9.6,1X))") i, ((real(i)-1)/(real(n)-1)), Eigenvalues(i)
      end do

      write(4,'(a15)') "# Eigenvectors"
      write(4, "(5X,30(I3,8X))") (vect(i), i=1,nvect)
      do i = 1, n
         write(4, "(I3,1X,30(F9.6,1X))") i, (Eigenvectors(i, vect(j)), j=1,nvect)
      end do

      close(2)
      close(3)
      close(4)
     
     end program main
      
      SUBROUTINE Diagonalization(Matrix,Eigenvalues,Eigenvectors,N)
!************************************************************************!
! The subroutine performs the diagonalization of a real symmetryc matrix !
!************************************************************************!

        implicit none
                                                       
        Real*8,dimension(N,N),intent(in)         :: Matrix 
        Real*8,dimension (N,N),intent(inout)     :: Eigenvectors
        Real*8,dimension(:,:),allocatable        :: Matrix_tmp
        Real*8,dimension(N),intent(inout)        :: Eigenvalues
        Real*8,dimension(:),allocatable          :: WORK
        integer,intent(in)                       :: N
        integer                                  :: INFO, LWORK

        allocate (Matrix_tmp(N,N))
        allocate (work(1))

        Matrix_tmp (:,:)=Matrix (:,:)

        call DSYEV('V','U',N,Matrix_tmp,N, & !We set LWORK to -1 in order
                   Eigenvalues,WORK,-1,INFO) !to get back the optimal 
                                             !size of Work.

        if (INFO .NE. 0) stop 

        LWORK=WORK(1) !As Lwork was set to -1
                      !WORK(1) returns the optimal LWORK
        deallocate (WORK)
        allocate (WORK(LWORK))

        Matrix_tmp (:,:)=Matrix (:,:) !The Matrix_tmp could have changed
                                      !and it would not be the entering
                                      !Hamiltonian anymore.

        call DSYEV('V','U',N,Matrix_tmp,N,Eigenvalues,WORK,LWORK,INFO)

        If (INFO .LT. 0) print*, "diagonalization failure: wrong argument"
        If (INFO .GT. 0) print*, "diagonalization failure: convergence not reached"
        if (INFO .NE. 0) stop

        EigenVectors(:,:)=Matrix_tmp(:,:)

        deallocate (work,matrix_tmp)

endsubroutine
