!--------------------------------------------------------------------
!        Calcula Energia libre
!--------------------------------------------------------------------

subroutine calc_energy(pHbulk)
#   include "control_run.h"
    use globales, only: delta, vsol, vsalt, radio, dimR , pi, lb, zwall
    use pore, only: qtot, psi, fdiswall
    use csys, only: eps1, q, sigma, sigmaq
    use FreeEnergy
    use mpmodule
    implicit none
    real(kind=8), intent(in) :: pHbulk
    real(kind=8) :: suma_pong, aux_mp
    real(kind=8) :: Free_Energy, Free_Energy2, F_Mix_s
!    real(kind=8) :: F_Mix_pos, F_Mix_neg, F_Mix_Hplus, F_Mix_OHmin, &
!                    F_Conf, F_Eq, F_Eq_wall, F_vdW, F_electro, F_eps
    integer :: iR,j

    print*, "Calculating energies..." 
    Free_Energy = 0.0
    Free_Energy2 = 0.0

! 1. Mezcla solvente
    F_Mix_s = 0.0 
    F_Mix_s = fmixs() ! 1. Mezcla solvente
    F_Mix_s = F_Mix_s * delta/vsol
    Free_Energy = Free_Energy + F_Mix_s

! 2. Mezcla ion positivo
    F_Mix_pos = 0.0 
    F_Mix_pos = fmixpos() ! 2. Mezcla ion positivo
    F_Mix_pos = F_Mix_pos * delta/vsol/vsalt
    Free_Energy = Free_Energy + F_Mix_pos

! 3. Mezcla ion negativo
    F_Mix_neg = 0.0
    F_Mix_neg = fmixneg()
    F_Mix_neg = F_Mix_neg * delta/vsol/vsalt
    Free_Energy = Free_Energy + F_Mix_neg

! 4. Mezcla protones
    F_Mix_Hplus = 0.0
    F_Mix_Hplus = fmixHplus()
    F_Mix_Hplus = F_Mix_Hplus * delta/vsol
    Free_Energy = Free_Energy + F_Mix_Hplus

! 5. Mezcla hidroxilos
    F_Mix_OHmin = 0.0
    F_Mix_OHmin = fmixOHmin()
    F_Mix_OHmin = F_Mix_OHmin * delta/vsol
    Free_Energy = Free_Energy + F_Mix_OHmin

! 6. Entropia interna polimero
    F_Conf = 0.0
# if CHAIN == 1 
    F_Conf = fconf_pol()
# endif
    Free_Energy = Free_Energy + F_Conf
!    print*, "E + F_Fconf" , Free_energy

! 7. Chemical Equilibrium
    F_Eq = 0.0 
    F_Eq_wall = 0.0
# if CHAIN == 1 
      F_Eq = fchem_eq()

      if ( sigmaq /= 0.0 ) then 
            print*, "sigmaq es distinto de cero!", sigmaq 
   !     F_Eq_wall = fchem_eq_wall()* sigmaq
      end if 
      F_Eq = F_Eq *delta/vsol
# endif
      F_Eq_wall = F_Eq_wall *delta/vsol
      Free_Energy = Free_Energy + F_Eq + F_Eq_wall
!    print*, "E + F_Eq" , Free_energy

! 8.vdW ! Ojo, los  son negativos => atraccion
    F_vdW = 0.0
#ifdef VDW
    F_vdW = fvdW()
    Free_Energy = Free_Energy + F_vdW
    print*, "E + F_vdW" , Free_energy
#endif

! 9. Electrostatic ! VER ESTO...
    F_electro = 0.0    
    do iR  = 1, dimR
        F_electro = F_electro + delta*psi(iR)*qtot(iR)/2.0/vsol * (dfloat(iR)-0.5)*delta/Radio
    enddo
       ! F_electro = F_electro + psi(dimR+1)*sigmaq*zwall
    Free_Energy = Free_Energy + F_electro
!    print*, "E + F_electro" , Free_energy
!*******************************************************************
! Observacion:
! La energia correspondiente a la distribucion de carga superficial 
! sigmaq NO es tenida en cuenta. Esto es por que la integral corespondiente
! a la energia electrostatica en Free_energy llega hasta R-delta/2
!*******************************************************************

! 10. Pol-Sup
    F_eps = 0.0
#if CHAIN == 1
    F_eps = fpol_sup() 
!    if(eps1.ne.0.0) then
!        print*, 'EPS should be 0 or correct free energy!'
!        stop
!    endif
    Free_Energy = Free_Energy + F_eps
#endif
    print*, "E " , Free_energy
 
! Segun Pong -  Estas son las condiciones de contorno del sistema?
! sumpi - 
! sumrho - 
! sumel :: No considered superficial charge distributio
    Free_Energy2 = 0.0
    suma_pong = pong_energy()
    
    aux_mp = 0.0
# if CHAIN == 1
    aux_mp = log(q)
# endif
    Free_Energy2 = suma_pong - (delta/vsol)*sigma*aux_mp !&- F_vdW !&
   !                 + (delta/vsol)*zwall*sigmaq*fdiswall/(4*pi*lb) * psi(dimR+1)/2  

!*******************************************************************
! Observacion:
! La energia correspondiente a la distribucion de carga superficial 
! sigmaq NO es tenida en cuenta. Esto es por que la integral corespondiente
! a la energia electrostatica en Free_energy llega hasta R-delta/2
!*******************************************************************

    Print*, "E2 ", Free_Energy2
    Print*, "diff", Free_Energy -  Free_Energy2

! Guarda energia libre
         write(301,*) pHbulk, Free_energy
         write(302,*) pHbulk, F_Mix_s 
         write(303,*) pHbulk, F_Mix_pos
         write(304,*) pHbulk, F_Mix_neg
         write(305,*) pHbulk, F_Mix_Hplus
         write(306,*) pHbulk, F_Mix_OHmin
         write(308,*) pHbulk, F_Eq, F_Eq_wall
#ifdef VDW
         write(309,*) pHbulk, F_vdW
#endif
# if CHAIN == 1
         write(307,*) pHbulk, F_Conf
         write(310,*) pHbulk, F_eps
# endif
         write(311,*) pHbulk, F_electro
         write(312,*) pHbulk, Free_energy2

!c--------------------------- FIN DE ENERGIA LIBRE -----------------

end subroutine calc_energy
