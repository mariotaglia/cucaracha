!
!

function factorcurv (i)
    use globales, only: pi, radio, delta 
    use csys, only: longporo
    implicit none
    integer, intent(in) :: i
    double precision :: factorcurv

    factorcurv = radio/((dfloat(i) - 0.5d0)*delta+radio)  
    factorcurv = factorcurv**2  ! elevamos al cuadrado por geometría esférica.
                         ! aparece en la expresion para la densidad de monomeros avpol

    return
end function factorcurv
