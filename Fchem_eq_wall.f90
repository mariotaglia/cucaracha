function fchem_eq_wall()
#   include "control_run.h"
    use globales, only: delta, vsol
    use csys, only: expmuHplus, Kwall0, sigmaq 
!    use FreeEnergy, only: checknumber
    use pore, only: fdiswall
    implicit none
    real(kind=8) :: fchem_eq_wall
    
    fchem_eq_wall = 0.0
# if fsigmaq == 1    
    if ( fdiswall  == 0.0 ) then 
        print*, "**********************************************************************************************"
        print*, "WARNING! fdiswall could not be zero! see chemical equilibrium in the wall: Fchem_eq_wall.f90"
        print*, "**********************************************************************************************"
    else 
        fchem_eq_wall = fchem_eq_wall + fdiswall  *(dlog(fdiswall/Kwall0)+dlog(expmuHplus))
    endif 
    if ( (1-fdiswall) == 0.0 ) then 
        print*, "**********************************************************************************************"
        print*, "WARNING! fdiswall could not be one! see chemical equilibrium in the wall: Fchem_eq_wall.f90"
        print*, "**********************************************************************************************"
    else 
        fchem_eq_wall = fchem_eq_wall + (1- fdiswall) * dlog(1-fdiswall)
    endif 
    fchem_eq_wall = fchem_eq_wall * sigmaq*(delta/vsol)
! *************************
!  Falta el potencial quimico standard del estado cero menos el pot.chem. del prontoN! (? que onda?)
! *************************

!    print*, "fchem_eq_wall llama checknumber: fchem_eq_wall", fchem_eq_wall
!    call checknumber(fchem_eq_wall)
#endif
    
    return
    contains 
        subroutine checknumber(var, arg)
            implicit none
            real(kind=8), intent(in) :: var
            character(len=*), intent(in) :: arg
            ! Check if Not a Number
            if ( var /= var ) then
                print*, arg, " real number is NaN"
                call printstate('checknumber real NaN')
                stop
            endif
        
            if ( var-1 == var ) then
                print*, arg, " real number is infinity"
                call printstate('checknumber real infinity')
                stop
            endif
        endsubroutine checknumber
end function 
