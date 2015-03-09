function pong_energy()
    use globales, only: delta, radio, dimR, vsalt, vsol
    use csys, only: sigmaq, xsolbulk, xHplusbulk, xOHminbulk, xnegbulk, xposbulk
    use pore, only: xh, xHplus, xOHmin, xpos, xneg, psi, qtot
    implicit none

    real(kind=8) :: pong_energy, sumpi, sumrho, sumel!, sum_disos
!    real(kind=8) :: aux_sumpi, aux_sumrho, aux_sumel!, sum_disos
    integer :: iR

    sumpi = 0.0
    sumrho=0.0
    sumel=0.0

!    aux_sumpi = 0.0
!    aux_sumrho=0.0
!    aux_sumel=0.0
!    sum_disos=0.0
! El calculo de esta energia es la diferencia entre la energia total del sistema 
! y la energia del bulk. Por eso aparecen las magnitudes en bulk. 
    do iR=1,dimR
        sumpi = sumpi+dlog(xh(iR))*(dfloat(iR)-0.5)*delta/Radio     
        sumpi = sumpi-dlog(xsolbulk)*(dfloat(iR)-0.5)*delta/Radio  
!        aux_sumpi = aux_sumpi+dlog(xh(iR))*(dfloat(iR)-0.5)*delta/Radio     
! x*** Son fracciones de volumen si queres las densidades de particulas entonces
! (xpos/vsalt)/vsol
        sumrho = sumrho + ( - xh(iR) -xHplus(iR) -xOHmin(iR) &
                            -(xpos(iR)+xneg(iR))/vsalt &
!                            - xpos2(iR)/vsalt2 & ! sum over  rho_i i=+,-,si
                          ) *(dfloat(iR)-0.5)*delta/Radio

        sumrho = sumrho - ( - xsolbulk -xHplusbulk -xOHminbulk &
                            -(xposbulk+xnegbulk)/vsalt &
!                            - xposbulk2/vsalt2 & ! sum over  rho_i i=+,-,si
                          ) *(dfloat(iR)-0.5)*delta/Radio
! electrostatic part free energy
        sumel = sumel - qtot(iR)*psi(iR)/2.0 *(dfloat(iR)-0.5)*delta/Radio  
! To consider Chemical equilibrium (facundo 12-2-2015)
!        sum_disos = sum_disos + ( dlog( fdis(iR)*Ka0 / (1-fdis(iR)) ) - 1 &
!                                ) * fdis(iR) * avpol(iR)/vpol *(dfloat(iR)-0.5)*delta/Radio

    enddo
!        aux_sumpi = aux_sumpi- radio* dlog(xsolbulk)/2.0  

!    sumel = sumel - sigmaq*psi(dimR)/2.0
!    print*, "surface charge: ", sigmaq*psi(dimR)/2.0

!    print*, " sumpi ", sumpi, " sumrho ", sumrho, " sumel ", sumel        

    sumpi = (delta/vsol)*sumpi
    sumrho = (delta/vsol)*sumrho
    sumel = (delta/vsol)*sumel

    pong_energy = sumpi + sumrho + sumel 

end function pong_energy
