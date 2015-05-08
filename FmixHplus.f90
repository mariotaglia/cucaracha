real(kind=8) function fmixHplus()
    use globales, only: dimR, Radio, delta, vsol
    use csys, only: expmuHplus, xHplusbulk
    use pore, only: xHplus 
!    use FreeEnergy
    implicit none
    integer :: iR
    fmixHplus=0
! Siempre se calcula la energia respecto de la de bulk!
    do iR = 1, dimR
        fmixHplus = fmixHplus + (xHplus(iR)/vsol)*(dlog(xHplus(iR))-1.0 -dlog(expmuHplus)) *delta*(dfloat(iR)-0.5)*delta/Radio
        fmixHplus = fmixHplus - (xHplusbulk/vsol)*(dlog(xHplusbulk)-1.0 -dlog(expmuHplus)) *delta*(dfloat(iR)-0.5)*delta/Radio
    enddo
!    print*, "fmixHplus llama checknumber: fmixHplus: ", fmixHplus
    call checknumber(fmixHplus, 'Energia fmixHplus')
    
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
