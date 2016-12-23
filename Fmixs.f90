function fmixs()
    use globales, only: dimR, Radio, delta, vsol
    use csys, only: xsolbulk
    use pore, only: xh 
!    use FreeEnergy, only: checknumber
    implicit none
    
    real(kind=8) :: fmixs
    integer :: iR
! El output de la funcion es la suma de todo esto
    fmixs=0
! Siempre se calcula la energia respecto de la de bulk
    do iR = 1, dimR
        fmixs = fmixs + ( xh(iR) /vsol) *( dlog(xh(iR)) -1.0) *delta*(dfloat(iR)-0.5+radio/delta)*delta/Radio 
        fmixs = fmixs - (xsolbulk/vsol)*(dlog(xsolbulk)-1.0) *delta*(dfloat(iR)-0.5+radio/delta)*delta/Radio
    enddo

!    print*, "fmixs llama checknumber: fmixs", fmixs
!    call checknumber(fmixs,'fmixs')
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

