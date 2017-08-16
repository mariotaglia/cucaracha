function pong_energy()
#   include "control_run.h"
    use globales, only: delta, radio, dimR, vsalt, vsol, vpol, pi, zwall, long
    use csys, only: xsolbulk, xHplusbulk, xOHminbulk, xnegbulk, xposbulk, xpolbulk, sigma, sigmaq, K_CL0, cpol
#ifdef PAHCL
    use pore, only: xh, xHplus, xOHmin, xpos, xneg, psi, qtot, avpol, xpot, fdis, fdis2,avpol!, fdiswall
#else  
    use pore, only: xh, xHplus, xOHmin, xpos, xneg, psi, qtot, avpol, xpot, fdis, avpol!, fdiswall
#endif
    implicit none

    real(kind=8) :: pong_energy, sumpi, sumrho, sumel, sumpol!, sum_disos
    integer :: iR
    real*8 :: rdrR

! Energía segun Pong:
! Esta energia se obtiene lugeo de reemplazar las expresiones 
! que se obtienen luego de minimizar la energía libre. 

    sumpi  = 0.0 ! sumpi :: osmotic pressure contribution (solvent expression)
    sumrho = 0.0  ! sumrho:: density of free ions contribution
    sumel  = 0.0   ! sumel :: No considered superficial charge distributio
    sumpol = 0.0  ! sumpol :: comes from polymer chain when they are free (not used when grafted chains)

!**********************************************************************************
! El calculo de esta energia es la diferencia entre la energia total del sistema 
! y la energia del bulk. Por eso aparecen las magnitudes en bulk. 
!**********************************************************************************
    do iR=1,dimR
        rdrR = (delta**3)*((dfloat(iR)-0.5+radio/delta)**2)/(Radio**2)
        sumpi = sumpi+dlog(xh(iR))/vsol *rdrR
        sumpi = sumpi-dlog(xsolbulk)/vsol *rdrR  
        
       sumrho = sumrho + ( - xh(iR) -xHplus(iR) -xOHmin(iR) &
                            -(xpos(iR)+xneg(iR))/vsalt &
!                            - xpos2(iR)/vsalt2 & ! sum over  rho_i i=+,-,si
                          ) /vsol *rdrR
! Bulk sumrho
        sumrho = sumrho - ( - xsolbulk -xHplusbulk -xOHminbulk &
                            -(xposbulk+xnegbulk)/vsalt &
!                            - xposbulk2/vsalt2 & ! sum over  rho_i i=+,-,si
                          ) /vsol *rdrR
! electrostatic part free energy
        sumel = sumel - qtot(iR) *psi(iR)/2.0 *rdrR

! Bulk polymer chain density
!        sumpol = sumpol - (- xpolbulk / (long*vpol) )/vsol *delta*(dfloat(iR)-0.5)*delta/Radio 
#if CHAIN != 0
! Writing output in std_mupol.dat ! Not clear that this line work! 
  !      write(202,*) iR, dlog(avpol(iR)/vpol)
#endif

    enddo
!  free polymer chain density
!        sumpol = sumpol - sigma*delta/vsol

!    print*, "surface charge: ", sigmaq*psi(dimR)/2.0
!    print*, " sumpi ", sumpi, " sumrho ", sumrho, " sumel ", sumel        
! output.aux
    write(324,*) cpol, sumpi, sumrho, sumel !, sumpol 
    pong_energy = sumpi + sumrho + sumel ! + sumpol 

end function pong_energy
