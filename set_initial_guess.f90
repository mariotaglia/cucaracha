! El input m es infile!:D
subroutine set_initial_guess(m)
#   include "control_run.h"
    use globales, only: dimR
    use csys, only: x1, xsolbulk, xflag
    use pore, only: xh, psi
    implicit none
    integer, intent(in) :: m
    integer :: i
    real(kind=8) :: c
    c=1.0
select case (m)
    case ( 0 )  ! Bulk solution
         xh(:dimR) = xsolbulk*c
!        xg1(:dimR) = xsolbulk
         x1(:dimR) = xsolbulk*c
! Psi Initial values
!          xg1(dimR+1:) = 0.0d0
!         psi(0:dimR+1) = 0.0d0 ! No hacen falta los exremos
!                               ! se corrigen en set_pore_distrib
         psi(1:dimR) = 0.0d0
         x1(dimR+1:) = 0.0d0
    case ( 1 ) 
         open(unit=92,file='initial_guess')
         read(92,*);read(92,*) ! ignore first and second lines
        
             do i=1,2*dimR
                 read(92,*) x1(i)
                 !read(92,*) xg1(i)
             end do
        
         close(92)
    ! xh and Psi Initial values are in xg1
          xh(:dimR) = x1(:dimR) 
          psi(1:dimR)= x1(dimR+1:)
    !      psi(0) = 0. ! No hace fatal se resuelve en set_pore_distrib
    !      psi(dimR+1) = 0. ! NO hace falta se resuelve en set_pore_distrib
          !x1(:) = xg1(:)
    case ( 2 )  ! Last solution
         x1(:) = xflag(:)
      !   xh(:dimR) = x1(:dimR) 
      !   psi(1:dimR)= x1(dimR+1:)
end select
! Por unica vez. Luego es llamada desde fkfun
# ifdef fdebug
    print*, "set_initial_guess.f90 L39: entro a set_poro_distrib"
# endif
    call set_pore_distrib
end subroutine set_initial_guess
