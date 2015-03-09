function fchem_eq_wall()
    use csys, only: expmuHplus, Kwall0 
!    use FreeEnergy, only: checknumber
    use pore, only: fdiswall
    implicit none
    real(kind=8) :: fchem_eq_wall
    
    fchem_eq_wall = 0.0
    
    fchem_eq_wall = fchem_eq_wall + fdiswall  *dlog(fdiswall/Kwall0)
    fchem_eq_wall = fchem_eq_wall + (1-fdiswall) *dlog(1.0-fdiswall)
    fchem_eq_wall = fchem_eq_wall + fdiswall*(dlog(expmuHplus))
! *************************
!  Falta el potencial quimico standard del estado cero minos el por chem del prontoN! (? que onda?)
! *************************

!    print*, "fchem_eq_wall llama checknumber: fchem_eq_wall", fchem_eq_wall
!    call checknumber(fchem_eq_wall)
    
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
