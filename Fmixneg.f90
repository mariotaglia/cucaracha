real(kind=8) function fmixneg()
    use globales, only: dimR, Radio, delta, vsalt
    use csys, only: expmuneg, xnegbulk
    use pore, only: xneg
!    use FreeEnergy, only: checknumber
    implicit none
    integer :: iR
    fmixneg=0
    do iR = 1, dimR
       fmixneg = fmixneg + xneg(iR)*(dlog(xneg(iR)/vsalt) -1.0 - dlog(expmuneg) + dlog(vsalt))*(dfloat(iR)-0.5)*delta/Radio
       fmixneg = fmixneg - xnegbulk*(dlog(xnegbulk/vsalt) -1.0 - dlog(expmuneg) + dlog(vsalt))*(dfloat(iR)-0.5)*delta/Radio
    enddo 

!    print*, "fmixneg llama checknumber: fmixneg: ", fmixneg
!    call checknumber(fmixneg)
    
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
