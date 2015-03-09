
! En este modulo van los vectores con informacion del pore

module pore
    use globales, only: dimR
    implicit none

! Volume fraction
!!     real(KIND=8), dimension(dimR)   :: xh     ! fraccion de volumen solvente
!!     real(KIND=8), dimension(dimR)   :: avpol  ! fraccion de volumen del polimero 
!!     real(kind=8), dimension(dimR)   :: xpot!, xpot_neg(dimR), xpot_pos(dimR)
!!     real(KIND=8), dimension(dimR)   :: qtot ! Carga total
!!     real(kind=8), dimension(dimR+1) :: psi ! psi se define asi por la discretizacion de Poisson y las boundary conditions 
!!     real(KIND=8), dimension(dimR)   :: xpos ! pos ion
!!     real(KIND=8), dimension(dimR)   :: xpos2 ! pos ion
!!     real(KIND=8), dimension(dimR)   :: xneg ! neg ioni
!!     real(KIND=8), dimension(dimR)   :: xHplus ! H+
!!     real(KIND=8), dimension(dimR)   :: xOHmin ! OH-
!!     real(kind=8), dimension(dimR)   :: xtotal ! xtotal poor solvent
!!     real(kind=8), dimension(dimR)   :: fdis, fdis2 !  weakpol
!!     real(kind=8) :: fdiswall 
!    Free_Energy, Free_Energy2, F_Mix_s
!    real(kind=8) :: F_Mix_pos, F_Mix_neg, F_Mix_Hplus, F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_electro, F_eps
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
    real(kind=8), dimension(:), allocatable :: fdis, fdis2 !  weakpol
    real(KIND=8), dimension(:), allocatable :: pp 
    real(kind=8) :: fdiswall 
end module pore                                                                                             
                                                                                                            
