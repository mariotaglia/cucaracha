function fchem_eq()
    use globales, only: dimR, delta, vpol, radio
    use csys, only: expmuHplus, Ka0 
!    use FreeEnergy, only: checknumber
    use pore, only: avpol, fdis!, fdis2
    implicit none
    real(kind=8) :: fchem_eq
    integer :: i
    
!    print*,  "WARNING! fchem_eq Calc. Be sure that neither f's is equal to zero!"
!    print*,  "Check fdis(), fdis2() and 1-f()-fdis2()" 
!    call printstate("fchem L14") ! Report of State
    
    fchem_eq = 0.0
    
    do i = 1, dimR
      
      fchem_eq = fchem_eq +    fdis(i)  *dlog(fdis(i)*Ka0)  *avpol(i)/vpol*(dfloat(i)-0.5)*delta/Radio
      fchem_eq = fchem_eq + (1-fdis(i)) *dlog(1.0-fdis(i))  *avpol(i)/vpol*(dfloat(i)-0.5)*delta/Radio
      fchem_eq = fchem_eq - fdis(i)*(dlog(expmuHplus)) *avpol(i)/vpol*(dfloat(i)-0.5)*delta/Radio ! ojo! expmuHplus ya tiene el signo menos!!
                  ! Falta el potencial quimico standard del polimero en el estado 0!!! (?)
    enddo

!    print*, "fchem_eq llama checknumber: fchem_eq", fchem_eq
!    call checknumber(fchem_eq,'fchem_eq')
    
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
