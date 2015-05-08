function fpol_sup()
    use globales, only: delta, radio, dimR
    use csys, only: eps1
    use pore, only: avpol
!    use FreeEnergy, only: checknumber
    implicit none
    real(kind=8) :: fpol_sup
    integer :: iR
!  2015-04-20 Interaccion polimero superficie (electrostatica(?))
    fpol_sup=0.0

    do iR = dimR-1, dimR
!    iR=dimR
            fpol_sup = fpol_sup - eps1*avpol(iR)*delta*(dfloat(iR)-0.5)*delta/radio
    enddo
    print*, "fpol_sup L16: " , fpol_sup 
!    print*, "f_vdW llama checknumber: f_vdW", f_vdW
    call checknumber(fpol_sup,'fpol_sup')
    
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
