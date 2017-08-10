real(kind=8) function fmixOHmin()
    use globales, only: dimR, Radio, delta, vsol
    use csys, only: expmuOHmin, xOHminbulk
    use pore, only: xOHmin
!    use FreeEnergy, only: checknumber
    implicit none
    integer :: iR
    fmixOHmin=0
! Siempre se calcula la energia respecto de la de bulk!
    do iR = 1, dimR
        fmixOHmin = fmixOHmin + (xOHmin(iR)/vsol)*(dlog(xOHmin(iR))-1.0 &
        -dlog(expmuOHmin)) *delta*(dfloat(iR)-0.5+radio/delta)*delta/Radio
        fmixOHmin = fmixOHmin - (xOHminbulk/vsol)*(dlog(xOHminbulk)-1.0 &
         -dlog(expmuOHmin)) *delta*(dfloat(iR)-0.5+radio/delta)*delta/Radio
! No deberia llevar un delta mas como aca: (?)
!        fmixOHmin = fmixOHmin + xOHmin(iR)*(dlog(xOHmin(iR))-1.0 -dlog(expmuOHmin)) *delta*(dfloat(iR)-0.5)*delta/Radio
!        fmixOHmin = fmixOHmin - xOHminbulk*(dlog(xOHminbulk)-1.0 -dlog(expmuOHmin)) *delta*(dfloat(iR)-0.5)*delta/Radio
    enddo

!    print*, "fmixOHmin llama checknumber: fmixOHmin", fmixOHmin
!    call checknumber(fmixOHmin)
    
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
