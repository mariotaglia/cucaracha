real(kind=8) function fmixpos()
    use globales, only: dimR, Radio, delta, vsalt, vsol
    use csys, only: expmupos, xposbulk
    use pore, only: xpos 
!    use FreeEnergy, only: checknumber
    implicit none
    integer :: iR
    fmixpos=0
! Siempre se calcula la energia respecto de la de bulk!
    do iR = 1, dimR
! Con los volumenes camibados (MARIO~)
        fmixpos = fmixpos + (xpos(iR)/(vsalt*vsol)) &
                  *(dlog(xpos(iR)/vsalt) -1.0 - dlog(expmupos/vsalt) ) *delta*(dfloat(iR)-0.5)*delta/Radio
        fmixpos = fmixpos - (xposbulk/(vsalt*vsol)) &
                  *(dlog(xposbulk/vsalt) -1.0 - dlog(expmupos/vsalt) ) *delta*(dfloat(iR)-0.5)*delta/Radio
!        fmixpos = fmixpos + (xpos(iR)/(vsalt*vsol))*(dlog(xpos(iR)) -1.0 - dlog(expmupos) ) *delta*(dfloat(iR)-0.5)*delta/Radio
!        fmixpos = fmixpos - (xposbulk/(vsalt*vsol))*(dlog(xposbulk) -1.0 - dlog(expmupos) ) *delta*(dfloat(iR)-0.5)*delta/Radio
! MI version 
!        fmixpos = fmixpos + (xpos(iR)/vsalt)*(dlog(xpos(iR)) -1.0 - dlog(expmupos) ) *(dfloat(iR)-0.5)*delta/Radio
!        fmixpos = fmixpos - (xposbulk/vsalt)*(dlog(xposbulk) -1.0 - dlog(expmupos) ) *(dfloat(iR)-0.5)*delta/Radio
    enddo

!    print*, "fmixpos llama checknumber: fmixpos = ", fmixpos
    call checknumber(fmixpos,'fmixpos')
    
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
