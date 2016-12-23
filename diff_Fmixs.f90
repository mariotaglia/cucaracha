function fmixs()
    use globales, only: delta, Radio, dimR
    use csys, only: xsolbulk
    use pore, only: xh
!    use FreeEnergy, only: checknumber
    implicit none
    
    real(kind=8) :: fmixs
    integer :: iR
    fmixs=0
    do iR = 1, dimR
        fmixs = fmixs + xh(iR)*(dlog(xh(iR))-1.0)*(dfloat(iR)-0.5+radio/delta)*delta/Radio 
        fmixs = fmixs - xsolbulk*(dlog(xsolbulk)-1.0)*(dfloat(iR)-0.5+radio/delta)*delta/Radio
    enddo

!    print*, "fmixs llama checknumber: fmixs", fmixs
    call checknumber(fmixs,'fmixs')
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
end function fmixs

