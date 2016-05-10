! This routine create the polymer/chain configurations
! WARNING: It call cadenas72.mr!
!    

subroutine creador ! crea cadenas
#   include "control_run.h"
! Las cadenas se colocan cerca de la posicion
! (x,y) = (0, -radio)
    use globales
    use csys
    implicit none
    interface
        subroutine cadenas72mr(cha, ncha)
            use globales
            integer :: ncha
!            real(kind=8), dimension(3,200,100), intent(out) :: cha
            real(kind=8), dimension(3,200,130), intent(out) :: cha
        endsubroutine
        subroutine imprimir_cadenas()
        endsubroutine
    endinterface
!    REAL(KIND=8)  indax, inday, indaz, rands
!    REAL(KIND=8), dimension(3,200,100) :: chains ! en un modulo especial?
    REAL(KIND=8), dimension(3,200,130) :: chains ! en un modulo especial?
    INTEGER :: i, il, j, ncha=0

    il=0

    do while (il.lt.cuantas)
        call cadenas72mr(chains,ncha)
        do i=1,ncha
            il = il + 1
            if(il.gt.cuantas) exit
            do j=1,long         
                in1(il,j,1)=chains(1,j,i)
                in1(il,j,2)=chains(2,j,i)
                in1(il,j,3)=chains(3,j,i)
            end do           
        end do
    end do

# ifdef fdebug
    print*, "creador.f90 L42: imprimo cadenas"
    call imprimir_cadenas
#endif
 100  return
end subroutine creador
