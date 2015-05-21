function fchem_eq()
#   include "control_run.h"
    use globales, only: dimR, delta, vpol, vsol, radio
    use csys, only: expmuHplus, expmuOHmin, Ka0 
!    use FreeEnergy, only: checknumber
    use pore, only: avpol, fdis!, fdis2
    implicit none
    real(kind=8) :: fchem_eq
    integer :: i
    
!    print*,  "WARNING! fchem_eq Calc. Be sure that neither f's is equal to zero!"
!    print*,  "Check fdis(), fdis2() and 1-f()-fdis2()" 
!    call printstate("fchem L14") ! Report of State
    
    fchem_eq = 0.0
    
    do i = 1, dimR
# if POL == 0 
! Nueva expresion para el equilibrio quimico
      fchem_eq = fchem_eq + ( fdis(i) *dlog( fdis(i)/Ka0 ) &
                         + (1-fdis(i))*dlog( (1-fdis(i)) ) )  *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio
      fchem_eq = fchem_eq + fdis(i)*dlog(expmuOHmin) *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio ! ojo! expmuHplus ya tiene el signo menos!!
! Original Mario
!      fchem_eq = fchem_eq + fdis(i)*dlog(fdis(i))            *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
!      fchem_eq = fchem_eq + (1.0-fdis(i))*dlog(1.0-fdis(i))  *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
!                   
!      fchem_eq = fchem_eq + (1.0-fdis(i))*dlog(Ka0)            *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
!      fchem_eq = fchem_eq + (1.0-fdis(i))*(-dlog(expmuOHmin)) *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
# elif POL == 1
! Here PMEP - TO CHECK!
      fchem_eq = fchem_eq + ( fdis(iR)*dlog(fdis(iR)/Ka0) &
                            + fdis2(iR)*dlog(fdis2(iR)/(Ka0*Kb0)) &
                            + (1.0-fdis(iR)-fdis2(iR))*dlog(1.0-fdis(iR)-fdis2(iR)) &
                            ) *avpol(iR)/vpol*(dfloat(iR)-0.5)*delta/Radio

      fchem_eq = fchem_eq + ( fdis(iR)*(dlog(expmuHplus)) & ! ojo! expmuHplus ya tiene el signo menos!!
                            + 2*fdis2(iR)*(dlog(expmuHplus)) )*avpol(iR)/vpol*(dfloat(iR)-0.5)*delta/Radio
# endif

    enddo

!    print*, "fchem_eq llama checknumber: fchem_eq", fchem_eq
!    call checknumber(fchem_eq,'fchem_eq')
    
    return
    contains 
!        function checknumber(var, arg) result(bool)
!            implicit none
!            real(kind=8), intent(in) :: var
!            character(len=*), intent(in) :: arg
!            logical :: bool 
!            ! Check if Not a Number
!            if ( var /= var ) then
!                print*, arg, " real number is NaN"
!                call printstate('checknumber real NaN')
!                bool=.true. 
!            endif
!        
!            if ( var-1 == var ) then
!                print*, arg, " real number is infinity"
!                call printstate('checknumber real infinity')
!                bool=.true.
!            endif
!        end function checknumber
!        subroutine checknumber(var, arg)
!            implicit none
!            real(kind=8), intent(in) :: var
!            character(len=*), intent(in) :: arg
!            ! Check if Not a Number
!            if ( var /= var ) then
!                print*, arg, " real number is NaN"
!                call printstate('checknumber real NaN')
!                stop
!            endif
!        
!            if ( var-1 == var ) then
!                print*, arg, " real number is infinity"
!                call printstate('checknumber real infinity')
!                stop
!            endif
!        endsubroutine checknumber
end function 
