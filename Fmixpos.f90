real(kind=8) function fmixpos()
    use globales, only: dimR, Radio, delta, vsalt
    use csys, only: expmupos, xposbulk
    use pore, only: xpos 
!    use FreeEnergy, only: checknumber
    implicit none
    integer :: iR
    fmixpos=0
    do iR = 1, dimR
        fmixpos = fmixpos + xpos(iR)*(dlog(xpos(iR)/vsalt)-1.0 - dlog(expmupos) + dlog(vsalt))*(dfloat(iR)-0.5)*delta/Radio
        fmixpos = fmixpos - xposbulk*(dlog(xposbulk/vsalt)-1.0 - dlog(expmupos) + dlog(vsalt))*(dfloat(iR)-0.5)*delta/Radio
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
