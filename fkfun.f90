!*************************************************************
! This subroutine define the function to minimize

subroutine fkfun(x,f,ier)
! fkfun define the function to minimize, f = f(x)
! No me gusta que setea los parametros del poro aca dentro.
! Esta funcion no deberia operar sobre los parametros del poro. 
! Su unica funcion deberia armar la funcion f = f(x) a minimizar por kinsol
!
#    include "control_run.h"
    use globales 
    use csys
    use pore
    use mpmodule
!    use interface_checknumber
    implicit none 

    real(KIND=8), intent(in), dimension(2*dimR) :: x ! son variables de salida también
    real(KIND=8), intent(out), dimension(2*dimR) :: f ! son variables de salida también
    integer(KIND=8), intent(out) :: ier
!    real(KIND=8) :: norma
    type (mp_real) norma ! este valor puede ser MUY grande
    integer :: i, iR

! Recupera xh y psi desde x()
! x(:) es parametro de entrada de fkfun es lo que devuelve kinsol (?)
    do iR=1,dimR
        xh(iR)=x(iR) ! Solvent volume fraction
        psi(iR)=x(iR+dimR) ! Electrical potential
    enddo
!print*, xh(:)
!print*, psi(:)
!    call printstate('fkfun 32')
!**********************************
! Boundary conditions are inside the function set_pore_distrib
    call set_pore_distrib ! Esto debe hacerse aca porque fkfun recibe el input!
!**********************************
        
# ifdef VDW
!**********************************
! Calculo de xtotal 
! PARA poor solvent en el lattice
    do iR=1,dimR
!        xtotal(iR) = 1.0 - xpos(iR) - xneg(iR) - xpos2(iR) - xh(iR) - xHplus(iR) - xOHmin(iR) ! xtotal es todo menos solvente e iones
        xtotal(iR) = 1.0 - xh(iR) - xpos(iR) - xneg(iR) - xHplus(iR) - xOHmin(iR) ! xtotal es todo menos solvente e iones (polimero?)
    enddo
!**********************************
#endif /* VDW */

!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------
! Volume fraction
    do iR=1,dimR
        f(iR)= xh(iR) & 
#if CHAIN != 0
             + avpol(iR) & 
#endif
             + xneg(iR) + xpos(iR) &
!             + xpos2(iR) &
             + xHplus(iR) + xOHmin(iR) - 1.000000d0
    enddo

! Poisson eq. ! poner una funcion (independizarme de la geometria!)
    do iR=1,dimR
! Plano
!               f(iR+dimR)=
!     & psi(iR+1) -2*psi(iR) + psi(iR-1) +
!     & qtot(iR)*constq

! Cilindro (derivada al centro), ver notas mis_apuntes.lyx eq.discret_poisson, f() debe ser cero
        f(iR+dimR)= psi(iR+1) -2*psi(iR) + psi(iR-1) &
                  + (0.5/(dfloat(iR)-0.5))*(psi(iR+1)-psi(iR-1)) & ! termino de curvatura
                  + qtot(iR)*delta*delta*constq
!*******************************************************************
! Observacion:
! La energia correspondiente a la distribucion de carga superficial 
! sigmaq NO es tenida en cuenta. Esto es por que la integral corespondiente
! a la energia electrostatica en Free_energy llega hasta R-delta/2
! Ver calculo de energias y energia de Pong.
!*******************************************************************
     
! elefante(7): porqué dividir por -2? 
        f(iR+dimR)=f(iR+dimR)/(-2.0) ! mejora kinsol...
    enddo
    
    iter = iter + 1
! Kinsol minimizes norma. ideally norma = 0.
!    norma = 0.0
    norma = '0.0'
    do i = 1, 2*dimR
        norma = norma +(f(i))**2
!        call checknumber(norma,'norma en fkfun') ! if NaN or Infinity checknumber stops the program
!        call Rchecknumber(norma, 'norma en fkfun') ! if NaN or Infinity checknumber stops the program
    enddo
!    print*, "fkfun67: iter, norma, q", iter, norma
    write(10,'(A14,I4,A8)', advance='no') "fkfun95: iter ", iter, " norma: "
    call mpwrite(10,norma)
    ier = 0
   
    return
    contains
        subroutine Rchecknumber(var,arg)
            implicit none
            real(kind=8), intent(in) :: var
            character(len=*), intent(in) :: arg
            ! Check if Not a Number
            if ( var /= var ) then
                print*, arg, " real number is NaN"
                call printstate(arg)
                !call printstate('checknumber real NaN')
                stop
            endif
        
            if ( var-1 == var ) then
                print*, arg, " real number is infinity"
                call printstate(arg)
                !call printstate('checknumber real infinity')
                stop
            endif
        
        endsubroutine Rchecknumber
end subroutine fkfun
