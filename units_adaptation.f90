! Here the main idea is to adapt the units of the input parameters

subroutine units_adaptation
#   include "control_run.h"
    use globales
    use csys
    implicit none

! there is many ways to define surface units
! Here we use vsol/delta because this is size independent?
! be carefull sigmaq is charge/surface as a function of solvent size and discretization size
    sigma = sigma *vsol/ delta
    sigmaq = sigmaq *vsol/delta ! Funciona! ver energia pong o 9.0 Electrostatic
    
    constq=delta*delta*lb ! multiplicative factor in poisson eq  ! SI UNITS

! Chemical Equilibrium
    Ka=10**(-pKa)
    Kb=10**(-pKb)
    Kwall=10**(-pKawall)

! Volume fraction second salt in mol/l, csalt2=0! Check!
    xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   
!      xsalt2=(csalt2*Na/(1.0d24))*(vsalt2*vsol)   
# if CHAIN == 1 && mupol == 1
    xpolbulk =(cpol*Na/(1.0d24))*(vpol*vsol)   
#endif
    print*, 'Units adaptation OK!'
end subroutine
