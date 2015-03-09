!
!

function factorcurv (i)
    use globales 
    implicit none
    integer, intent(in) :: i
    double precision :: factorcurv

    factorcurv = radio/((dfloat(i) - 0.5d0)*delta)

    return
end function factorcurv
