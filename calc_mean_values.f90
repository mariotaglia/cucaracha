! Calcula fmedio y Rmedio
!-----------------------------------------------------------------

subroutine calc_mean_values(pHbulk)
#   include "control_run.h"
    use globales, only: delta, vpol, vsol, zpos, pi, lb, radio, zwall
    use csys, only: sigmaq
    use pore
    implicit none
    real(kind=8), intent(in) :: pHbulk
    integer :: iR
    real(kind=8) :: sumpol, fmedio, fmedio2, fdisw, Rmedio, sumcharge

    print*, "Calculating mean values... "
        fdisw = 0.0
#if CHAIN == 1
       fmedio = 0.0
       sumpol = 0.0
      fmedio2 = 0.0
       Rmedio = 0.0
    sumcharge = 0.0

    do iR = 1, dimR
        sumpol = sumpol + avpol(iR)*(dfloat(iR)-0.5) ! suma de avpol*r (falta constante)
!       sumpolp = sumpolp +
!     & avpol(iR)*(dfloat(iR)-0.5) ! suma de avpol*r (falta constante)
!       sumpoln = sumpoln +
!     & avpol(iR)*(dfloat(iR)-0.5) ! suma de avpol*r (falta constante)

         fmedio = fmedio + fdis(iR)*avpol(iR)*(dfloat(iR)-0.5) ! suma de frdis*avpol*r (falta constante)
        fmedio2 = fmedio2 + fdis2(iR)*avpol(iR)*(dfloat(iR)-0.5) ! suma de frdis*avpol*r (falta constante)
!       fmedio = fmedio +
!     & fdis(iR)*avpoln(iR)*(dfloat(iR)-0.5) ! suma de frdis*avpol*r (falta constante)
!       fmedio2 = fmedio2 +
!     & fdis2(iR)*avpolp(iR)*(dfloat(iR)-0.5) ! suma de frdis*avpol*r (falta constante)
          Rmedio = Rmedio + avpol(iR)*(((dfloat(iR)-0.5)*delta)**2) ! suma de avpol*r^2 (falta constante)
        sumcharge=sumcharge + (dfloat(iR)-0.5)*(avpol(iR)/(vpol*vsol))*zpos *(fdis(iR)+2*fdis2(iR))
!       sumcharge=sumcharge + (dfloat(iR)-0.5)*(avpoln(iR)/vpol)*fdis(iR)*zneg
!     & + (dfloat(iR)-0.5)*avpolp(iR)*fdis2(iR)*zpos ! para switterion
    enddo
! the idea is to calculate the net charge in the wall?
!        fdisw = (sigmaq /(4.0*pi*lb*delta) ) *fdiswall ! ??
!        fdisw = sigmaq*2*pi*radio*longporo /(4.0*pi*lb*delta) ) *fdiswall ! ??
    if ( sumpol /= 0.0 ) then ! si sumpol es distinto de cero
           fmedio = fmedio/sumpol    ! fraccion desprotonada media por monomero
          fmedio2 = fmedio2/sumpol    ! fraccion desprotonada media por monomero
           Rmedio = Rmedio/sumpol
        sumcharge = sumcharge/sumpol ! carga total por unidad de monomero
    end if

    write(318,*) pHbulk, Rmedio, sumcharge
! NOTE: fdiswall was calculated in set_pore_distrib be carefull! ;)
    write(313,*) pHbulk, fmedio, fdiswall, sigmaq*zwall*fdiswall !, fmedio2, fdiswall !, fdisw <- cual es el sentido de fdisw?
#else
    write(313,*) pHbulk, fdiswall, sigmaq*fdiswall*zwall
#endif

end subroutine calc_mean_values
