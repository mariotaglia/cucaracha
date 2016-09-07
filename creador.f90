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
        subroutine cadenas72mr_rosen(cha, ncha, wcha)
            use globales
            integer :: ncha
!            real(kind=8), dimension(3,200,100), intent(out) :: cha
            real(kind=8), dimension(3,200,130), intent(out) :: cha
            real(kind=8), dimension(200), intent(out) :: wcha
        endsubroutine
        subroutine imprimir_cadenas()
        endsubroutine
    endinterface
!    REAL(KIND=8)  indax, inday, indaz, rands
!    REAL(KIND=8), dimension(3,200,100) :: chains ! en un modulo especial?
    REAL(KIND=8), dimension(3,200,130) :: chains ! en un modulo especial?
    REAL(KIND=8), dimension(200) :: wchains ! en un modulo especial?
    INTEGER :: i, il, j, ncha
    real*8 inw, in1(long,3)

    ncha=0
    sumallweight = 0.0
    chaintot = 0
    do while (chaintot.lt.cuantas)
        if(rosen.eq.1) then
        call cadenas72mr_rosen(chains,ncha, wchains) ! ncha: number of generated chains
        else
        call cadenas72mr(chains,ncha,wchains) ! ncha: number of generated chains
        endif
!        print*, ncha, wchains(1)
        do i=1,ncha
            if(chaintot.gt.cuantas) exit
            do j=1,long         
                in1(j,1)=chains(1,j,i)
                in1(j,2)=chains(2,j,i)
                in1(j,3)=chains(3,j,i)
            end do           
                inw = wchains(i)
                call pxs(in1,inw)
        end do
    end do
    if(rosen.eq.-1) pR(chaintot,:) = dimR

    print*, 'sumall weight = ', sumallweight
# ifdef fdebug
    print*, "creador.f90 L42: imprimo cadenas"
    call imprimir_cadenas(in1)
#endif
 100  return
end subroutine creador
