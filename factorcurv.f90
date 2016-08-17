!
!

function factorcurv (i)
    use globales, only: pi, radio, delta 
    use csys, only: longporo
    implicit none
    integer, intent(in) :: i
    double precision :: factorcurv

    factorcurv = radio/((dfloat(i) - 0.5d0)*delta) ! /(2*pi*longporo) <- este termino se cancela con el que 
                                                   ! aparece en la expresion para la densidad de monomeros avpol

    return
end function factorcurv
