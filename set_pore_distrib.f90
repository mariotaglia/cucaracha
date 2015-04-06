subroutine set_pore_distrib
#   include "control_run.h" ! To control the code to be included in the compilation
    use mpmodule
    use globales
    use csys
    use pore
    implicit none
    interface    
        function factorcurv(i)
            integer, intent(in) :: i
            double precision :: factorcurv
        end function factorcurv
    end interface 
    integer :: iR, aR, i,j 
    real(kind=8) :: shift ! it is  multiplicative factor in probability pro(cuantas)
    type (mp_real) q ! este valor puede ser MUY grande
    real(kind=8) :: aux_mp
! Dissociation in the inner wall of the pore 
! + interaction with the electrostatic Potential    
!    fdiswall = Kwall0 / (expmuHplus*dexp(psi(dimR+1)*sigmaq*zwall) + Kwall0) 
!    fdiswall = Kwall0 / (expmuHplus*dexp(psi(dimR)*sigmaq*zwall*2*pi*radio) + Kwall0) 
    fdiswall = Kwall0 / (expmuHplus*dexp(psi(dimR)*zwall) + Kwall0) 

! Boundary Conditions: Electrical Potential
!! Estas ecuaciones son suplementarias a la eq. de Poisson discretizada
!   *   La derivada en r=0 tiene que ser cero. Se utiliza orden uno en la derivada en r=0.
!   *   La derivada en r=R considera la carga superficial sigmaq. Se utiliza orden 2 en la derivada.
    psi(0) = psi(1) ! La derivada en r = 0 es igual a cero
!    psi(dimR+1) = psi(dimR-1) + 2*delta*sigmaq*zwall*fdiswall ! La derivada en r=R es el salto de la carga superficial, ver notas
!    psi(dimR+1) = psi(dimR-1) + 2*delta*sigmaq*zwall*fdiswall*2*pi*radio ! La derivada en r=R es el salto de la carga superficial, ver notas
    psi(dimR+1) = psi(dimR) + (4*pi*lb*delta)*sigmaq*zwall*fdiswall ! La derivada en r=R es el salto de la carga superficial, ver notas
    
!    call printstate("set_pore_distrib L23")

! Volume fractions and degrees of dissociations: fdis
    do iR=1,dimR
! Mis expresiones:
! expmu**** = son las expresiones de los pot. quimicos con las concentraciones de BULK,
! * xh : concentracion de solvente, tiene una expresion en terminos de pi(r).
        xpos(iR) = expmupos*(xh(iR)**vsalt) *dexp(-psi(iR)*zpos) ! ion plus volume fraction

        xneg(iR) = expmuneg*(xh(iR)**vsalt) *dexp(-psi(iR)*zneg) ! ion neg volume fraction
     
        xHplus(iR) = expmuHplus*(xh(iR)**vHplus) *dexp(-psi(iR)*zH)           ! H+ volume fraction
     
        xOHmin(iR) = expmuOHmin*(xh(iR)**vOHmin) *dexp(-psi(iR)*zOH)           ! OH-  volume fraction
!        fdis(iR) = dissos_degree(1,iR)
! Para estas expresiones es necesario definir las constantes K1 y K2
!          denon = ( expmuOHmin*dexp(psi(iR)*zpos) )**2
!     &          + expmuOHmin*dexp(psi(iR)*zpos)*Ka0
!     &          + Ka0*Kb0 ! Kb0 = Ka20 for dobule protonation
#   if CHAIN == 1
!         fdis(iR) = expmuOHmin*dexp(psi(iR)*zpos)*Ka0 / denon 
!         fdis2(iR) = Ka0*Kb0 / denon ! Kb0 = Ka20 for dobule protonation
        fdis(iR) = 1.0d0/(1.0d0 + (Ka0*dexp(psi(iR)*zpos)/expmuHplus) )
!       fdis2(iR) = 1.0d0/(1.0d0 + (dexp(psi(iR)*zpos)*expmuOHmin/Kb0) )
        fdis2(iR) = 0.0
#   endif
    end do

# if CHAIN == 1
!!!!!! AQUI FALTA ACTUALIZAR xH! debe tomar el valor de  x1
!----------------------------------------------------------------------------- 
! To calculate VdW force coefficient
!----------------------------------------------------------------------------- 
    do iR = 1, dimR
! El termino: (xh(iR)**vpol) , viene de reemplazar 
! la presion osmotica por la expresion para el solvente
        xpot(iR) = (xh(iR)**vpol) /(1-fdis(iR)-fdis2(iR) )
# ifdef VDW
! Taking in to account the probability with the vdW force coefficient for every layer.
        do j = 1, dimR
!        xpot(iR) = xpot(iR) * dexp(vsol*vpol* Xu(iR,j)*xtotal(j+1)) ! Porque j+1 ??? (elefante)
            xpot(iR) = xpot(iR) * dexp(vsol*vpol* Xu(iR,j)*xtotal(j))
        end do
# endif
    enddo
!          q = 0.0 ! Normalization of P(alpha)
          q = '0.0' ! Normalization of P(alpha) remember that: type (mp_real) q !!
      shift = 1.0 ! (?)
    avpol(:)= 0.0 ! line important to probability calculus
! Estos Do's estan bien! siempre dejar la coma a la izquierda
    do i=1,cuantas ! i enumerate configurations (configurations ensamble)
        pro(i)=1.0!*shift
!         do j=1,long, 2 ! (long=28) ! Be carefull! Here you choose the type of the first segment
        do j=1,long ! (long=28)
            aR = pR(i, j) ! pR()'s output is the layer where is the segment j of configuration i.
!            bR = pR(i, j+1)   
! The configuration's probability is the product of layer's probability
! many times as particles in that layer. 
            pro(i) = pro(i) * xpot(aR)
! ATENTTION: Now you insert structural detail about alternating type of 
! monomers, A - B - A - B - etc
!            pro(i) = pro(i) * xpot_neg(aR) * xpot_pos(bR)
        enddo
! q es la suma de todas las probabilidades
        q=q+pro(i)
        do j = 1,long
            aR = pR(i, j) ! pR devuelve en que layer se encuentra el monómero j del polímero i.
            avpol(aR) = avpol(aR) + pro(i)*sigma*vpol*factorcurv(aR) ! cilindro, ver notas...
! ver eq:factorcurv in mis_apuntes.lyx
!            bR = pR(i, j+1)
!            avpol(aR) = avpol(aR) + pro(i)*sigma*vpol ! plano
! sigma*vpol*Factorcurv(aR) es el factor en la densidad de polimero.
!            avpoln(aR) = avpoln(aR) + pro(i)*sigma*vpol*Factorcurv(aR) ! cilindro, ver notas...
!            avpolp(bR) = avpolp(bR) + pro(i)*sigma*vpol*Factorcurv(bR) ! cilindro, ver notas...
        end do
    enddo ! End loop over chains/configurations
    call mpwrite(11,q) 
    log_q = log(q) ! Variable clave en el calculo de energias!
    print*, "Se normalizan todas las probabilidades por q!"
    do i=1,cuantas 
        pro(i) = pro(i)/q
    enddo
    print*, "pro(:): ", pro(:)
!    call printstate("L105 set_pore_distrib")
    do iR=1, dimR            ! norma avpol
        aux_mp =  avpol(iR)/ q  ! q_mp_real
        avpol(iR) = aux_mp
!        avpol(iR) = avpol(iR) / q
!        avpoln(iR) = avpoln(iR) / q
!        avpolp(iR) = avpolp(iR) / q 
!         avpol(iR) = avpoln(iR) + avpolp(iR)
    enddo
! Chequeo q
!    call checknum(q,'q en set_pore_distrib')
!    call printstate("L115 set_pore_distrib")
# endif

! Local charge
    do iR=1,dimR  
        ! Ojo! volver aca si no hay carga superfical!
!    do iR=1,dimR-1
!         qtot(iR) = 
!     &   (zpos*xpos(iR)+zneg*xneg(iR))/vsalt + zpos2*xpos2(iR)/vsalt2+
!     &   avpol(iR)*(zpol/vpol)*fdis(iR) + xHplus(iR)-xOHmin(iR) 

        qtot(iR) = (zpos*xpos(iR)+zneg*xneg(iR))/vsalt &
# if CHAIN == 1
            + (fdis(iR) + 2*fdis2(iR))*zpos*avpol(iR)/vpol &
# endif
            + xHplus(iR)-xOHmin(iR) 
!         + zpos2*xpos2(iR)/vsalt2 &
!     &   + (avpolp(iR)*fdis2(iR)*zpos + avpoln(iR)*fdis(iR)*zneg)/vpol
    enddo

end subroutine set_pore_distrib
