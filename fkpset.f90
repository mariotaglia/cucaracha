!*******************************************************************************
! The routine kpreco is the preconditioner setup routine. It must have
! that specific name be used in order that the c code can find and link
! to it.  The argument list must also be as illustrated below:
!

subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)
!          FKPSET(  U  , USCALE,  FVAL, FSCALE,VTEMP1,VTEMP2, IER)
! It must perform any evaluation of Jacobian-related data and preprocessing needed for the solution
! of the preconditioned linear systems by FKPSOL. The variables U through FSCALE are for use in
! the preconditioning setup process. Typically, the system function FKFUN is called before any calls
! to FKPSET, so that FVAL will have been updated. U is the current solution iterate. The arrays
! VTEMP1 and VTEMP2 are available for work space. If scaling is being used, USCALE and FSCALE are
! available for those operations requiring scaling
!
    use globales, only: neq
    use pore, only: pp
    implicit none
    real(KIND=8) :: udata(*), uscale(*), fdata(*), fscale(*)
    real(KIND=8) :: vtemp1(*), vtemp2(*)
    integer,intent(out):: ier
    
    integer*8 :: i
!    real(KIND=8),dimension(2*dimR) :: pp ! pp difined in kinsl_param
! double precision == real(KIND=8)

!    common /pcom/ pp(2*dimR)
! neq = 2*dimR in module globales
    do i = 1, neq/2
       pp(i) = 0.1 / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
    enddo

    do i = neq/2+1, neq
       pp(i) = 1.0
    enddo

    ier = 0

    return
end subroutine fkpset
