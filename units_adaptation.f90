subroutine units_adaptation
#   include "control_run.h"
    use globales
    use csys
    implicit none

! Here the main idea is to adapt the units of the input parameters
! there is many ways to define surface units
! Here we use vsol/delta because this is size independent?
! be carefull sigmaq is charge/surface as a function of solvent 
! size and discretization size
   vsigma(:) = vsigma(:) *vsol/ delta ! La superficie de referencia es vsol/delta
   sigmaq = sigmaq *vsol/delta ! Funciona! ver energia pong o 9.0 Electrostatic

! Electrostatic
    constq=(4*pi*lb) ! multiplicative factor in poisson eq  ! SI UNITS

! Chemical Equilibrium
    Ka=10**(-pKa)
    Kb=10**(-pKb)
    K_Cl=10**(-pK_Cl)
    Kwall=10**(-pKawall)
    print*, "Kwall, K_Cl, Ka, Kb : " , Kwall, K_Cl, Ka, Kb
    print*, 'Units adaptation OK!'
end subroutine
