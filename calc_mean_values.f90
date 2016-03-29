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
#if CHAIN != 0 
   fmedio = 0.0
   sumpol = 0.0
  fmedio2 = 0.0
   Rmedio = 0.0
sumcharge = 0.0

do iR = 1, dimR
    sumpol = sumpol + avpol(iR)*(dfloat(iR)-0.5) ! suma de avpol*r (falta constante)
!sumpolp = sumpolp + avpolp(iR)*(dfloat(iR)-0.5) ! suma de avpolp*r (falta constante)
!sumpoln = sumpoln + avpoln(iR)*(dfloat(iR)-0.5) ! suma de avpoln*r (falta constante)

!   suma de avpol*r^2 (falta constante)
    Rmedio = Rmedio + avpol(iR)*(((dfloat(iR)-0.5)*delta)**2) 
!   suma de frdis*avpol*r (falta constante)
    fmedio = fmedio + fdis(iR)*avpol(iR)*(dfloat(iR)-0.5) 

! suma de frdis*avpoln*r (falta constante)
!    fmedio = fmedio + fdis(iR)*avpoln(iR)*(dfloat(iR)-0.5) 
#ifdef PAHCL
    fmedio2 = fmedio2 + fdis2(iR)*avpol(iR)*(dfloat(iR)-0.5) 
#endif
# if POL == 0 || POL == 2/* PAH */
    sumcharge=sumcharge + (dfloat(iR)-0.5)*(avpol(iR)/(vpol*vsol))*zpos *(fdis(iR))
# elif POL == 1 /* PMEP */
  ! suma de frdis*avpol*r (falta constante)
    fmedio2 = fmedio2 + fdis2(iR)*avpol(iR)*(dfloat(iR)-0.5) 
!    fmedio2 = fmedio2 + fdis2(iR)*avpolp(iR)*(dfloat(iR)-0.5) 
    sumcharge=sumcharge + (dfloat(iR)-0.5)*(avpol(iR)/(vpol*vsol))*zpos &
                                          *(fdis(iR)+2*fdis2(iR))
!  sumcharge=sumcharge + (dfloat(iR)-0.5)*(avpoln(iR)/vpol)*fdis(iR)*zneg &
!                      + (dfloat(iR)-0.5)*avpolp(iR)*fdis2(iR)*zpos ! para switterion
# endif /* PAH || PMEP */
! Estaria faltando una constante porque se calcula el numero de cargas.
! Pero todas se normalizan por el numero de cargas totales (sumpol)
enddo
! the idea is to calculate the net charge in the wall?
!    fdisw = (sigmaq /(4.0*pi*lb*delta) ) *fdiswall ! ??
!    fdisw = sigmaq*2*pi*radio*longporo /(4.0*pi*lb*delta) ) *fdiswall ! ??
if ( sumpol /= 0.0 ) then ! si sumpol es distinto de cero
      fmedio = fmedio/sumpol    ! fraccion desprotonada media por monomero
      Rmedio = Rmedio/sumpol
! carga polimerica total por unidad de superficie 
   sumcharge = (delta/radio)*sumcharge 
# if POL == 1
     fmedio2 = fmedio2/sumpol    ! fraccion desprotonada media por monomero
# endif
end if

write(318,*) pHbulk, sigmaq*(delta/vsol)*zwall*fdiswall +sumcharge, &
                     sigmaq*(delta/vsol)*zwall*fdiswall, sumcharge
!    write(318,*) pHbulk, Rmedio, sigmaq*(delta/vsol)*zwall*fdiswall +sumcharge, &
!                                 sigmaq*(delta/vsol)*zwall*fdiswall, sumcharge
! NOTE: fdiswall was calculated in set_pore_distrib be carefull! ;)
#   if POL == 0 /* PAH */
#ifdef PAHCL
! writing fmedio.dat:
    write(313,*) pHbulk, fmedio, fmedio2, fdiswall !, fmedio2, fdiswall
#else
! writing fmedio.dat:
    write(313,*) pHbulk, fmedio, fdiswall !, fmedio2, fdiswall
#endif
#   elif POL == 1 /* PMEP */
! writing fmedio.dat:
    write(313,*) pHbulk, fmedio, fmedio2, fdiswall
#   elif POL == 2 /* Neutral Polymer */
    write(313,*) pHbulk, fdiswall !, fmedio2, fdiswall
#   endif /* PAH || PMEP */

#else
! writing fmedio.dat:
write(313,*) pHbulk, fdiswall, (sigmaq*delta/vsol)*fdiswall*zwall, sumcharge
!write(313,*) pHbulk, fdiswall, (sigmaq)*fdiswall*zwall, sumcharge
#endif

end subroutine calc_mean_values
