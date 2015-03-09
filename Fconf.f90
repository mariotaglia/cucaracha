function fconf_pol()
    use mpmodule
    use globales, only: cuantas, delta, vsol
    use csys, only: q, pro, sigma 
!    use FreeEnergy, only: checknumber
    implicit none
    real(kind=8) :: fconf_pol, aux_mp
    integer :: i
    fconf_pol=0

    do i = 1, cuantas
        aux_mp = log(pro(i)/q)
        fconf_pol = fconf_pol + (pro(i)/q)*aux_mp /vsol*delta*sigma
        !fconf_pol = fconf_pol + (pro(i)/q)*dlog((pro(i))/q)  /vsol*delta*sigma
    enddo

!    print*, "fconf_pol llama checknumber: fconf_pol", fconf_pol
    call checknumber(fconf_pol,'fconf_pol')
    
    return
    contains 
        subroutine checknumber(var, arg)
            implicit none
            real(kind=8), intent(in) :: var
            character(len=*), intent(in) :: arg
            ! Check if Not a Number
            
            if ( var /= var ) then
                print*, arg, " real number is NaN"
!                call printstate('checknumber real NaN')
                call printstate(arg)
                stop
            endif
        
            if ( var-1 == var ) then
                print*, arg, " real number is infinity"
                call printstate('checknumber real infinity')
                stop
            endif
        endsubroutine checknumber
end function 
