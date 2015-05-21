
! En este modulo van los vectores con informacion del pore

module pore
#   include "control_run.h"
    use globales, only: dimR
    implicit none

! Volume fraction
    real(KIND=8), dimension(:), allocatable :: xh          ! fraccion de volumen solvente
    real(KIND=8), dimension(:), allocatable :: avpol       ! fraccion de volumen del polimero 
    real(kind=8), dimension(:), allocatable :: xpot        !, xpot_neg(dimR), xpot_pos(dimR)
    real(KIND=8), dimension(:), allocatable :: qtot        ! Carga total
    real(kind=8), dimension(:), allocatable :: psi         ! psi se define asi por la discretizacion de Poisson y las boundary conditions 
    real(KIND=8), dimension(:), allocatable :: xpos        ! pos ion
!    real(KIND=8), dimension(:), allocatable :: xpos2       ! pos ion
    real(KIND=8), dimension(:), allocatable :: xneg        ! neg ioni
    real(KIND=8), dimension(:), allocatable :: xHplus      ! H+
    real(KIND=8), dimension(:), allocatable :: xOHmin      ! OH-
    real(kind=8), dimension(:), allocatable :: xtotal      ! xtotal poor solvent
# if POL == 0
    real(kind=8), dimension(:), allocatable :: fdis ! PAH weakpol
# elif POL == 1
    real(kind=8), dimension(:), allocatable :: fdis, fdis2 ! PMEP weakpol
# endif
    real(KIND=8), dimension(:), allocatable :: pp 
    real(kind=8) :: fdiswall 
end module pore                                                                                             
                                                                                                            
