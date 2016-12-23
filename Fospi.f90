real(kind=8) function fospi()
    use globales, only: dimR, Radio, delta, vsalt, vpol, vsol
    use csys, only: expmuOHmin, xOHminbulk
    use pore, only: avpol, xOHmin, xHplus, xpos, xneg, xh
!    use FreeEnergy, only: checknumber
    implicit none
    integer :: iR
    fospi=0
! Siempre se calcula la energia respecto de la de bulk!
! Osmotic Preassure Energy (?) constraints! should be zero
    do iR = 1, dimR
        fospi = fospi + ( xOHmin(iR) + xHplus(iR) + xh(iR) &
                      +  (xneg(iR) + xpos(iR) ) &
                      + avpol(iR) - 1.0 ) &
                      *(dlog(xh(iR))/vsol) *delta*(dfloat(iR)-0.5+radio/delta)*delta/Radio
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
