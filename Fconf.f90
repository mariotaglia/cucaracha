function fconf_pol()
    use mpmodule
    use globales, only: chaintot, delta, vsol
    use csys, only: log_q, pro, sigma 
!    use FreeEnergy, only: checknumber
    implicit none
    real(kind=8) :: fconf_pol, aux_mp
    integer :: i
    fconf_pol=0
! Siempre se calcula la energia respecto de la de bulk!
!    print*, "Fconf10: Valor log_q", log_q
    do i = 1, chaintot
!        aux_mp = log(pro(i)/q)
!        aux_mp = log(pro(i)) - log_q
!        fconf_pol = fconf_pol + (pro(i)/q)*aux_mp /vsol*delta*sigma
!        fconf_pol = fconf_pol + sigma*delta*pro(i)*log(pro(i))/vsol
        fconf_pol = fconf_pol + (pro(i))*log(pro(i)) 
        !fconf_pol = fconf_pol + (pro(i)/q)*dlog((pro(i))/q)  /vsol*delta*sigma
    enddo
! Aca mp deneroa janer ima coreccion? elefante
!        fconf_pol = (sigma*delta/vsol)*fconf_pol  

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
