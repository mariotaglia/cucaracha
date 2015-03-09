function fvdW()
    use globales, only: delta, radio, dimR
    use csys, only: Xu
    use pore, only: xtotal
!    use FreeEnergy, only: checknumber
    implicit none
    real(kind=8) :: fvdW
    integer :: iR,j
    
    fvdW=0.0

    do iR = 1, dimR
        do j = 1, dimR
            fvdW = fvdW - 0.5000*delta*xtotal(iR)*xtotal(j)*Xu(iR,j) *(dfloat(iR)-0.5)*delta/radio
        enddo
    enddo

!    print*, "f_vdW llama checknumber: f_vdW", f_vdW
    call checknumber(fvdW, 'fvdW')
    
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
