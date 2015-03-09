!*****************************************************************************
!     The routine kpsol is the preconditioner solve routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:
!

subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
!
! Typically this routine will use only U, FVAL, VTEM and FTEM. It must solve the preconditioned linear
! system P z = r, where r = VTEM is input, and store the solution z in VTEM as well. Here P is the
! right preconditioner. If scaling is being used, the routine supplied must also account for scaling
! on either coordinate or function value, as given in the arrays USCALE and FSCALE, respectively.
!
    use globales, only: neq, dimR
!    use kinsol_param
    use pore, only: pp
    implicit none

    integer, intent(out) :: ier
    integer(kind=8) :: i
!    real(KIND=8) udata(*), uscale(*), fdata(*), fscale(*)
    real(KIND=8), dimension(2*dimR) :: udata, uscale, fdata, fscale
    real(KIND=8), dimension(2*dimR) :: vv, ftem

    do  i = 1, neq
       vv(i) = vv(i) * pp(i)
    enddo
    ier = 0

    return
end subroutine fkpsol
