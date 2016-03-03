![KINSOL ERROR]  KINSetMaxSetupCalls
!  Illegal msbset < 0.
!
! ******************************************************
!
! Subrutina que llama a kinsol
!
! ******************************************************

subroutine call_kinsol(x1, ier)
!subroutine call_kinsol(x1, xg1, ier)
#   include "control_run.h"
    use globales
    implicit none
    
interface
! Paso 1 - Se especifica la funcion a resolver mediante fkfun. fval = f(u)
!    SUBROUTINE FKFUN (U, FVAL, IER)
!    DIMENSION U(*), FVAL(*)
    subroutine fkfun(x,f,ier)
    real(KIND=8), intent(out), DIMENSION(:) :: x,f ! son variables de salida también
    integer(KIND=8), intent(out) :: ier
    end subroutine fkfun
    ! fkpset: Se especifica el jacobiano a utilizar. fdata es la f de fkfun. 
    subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)
    real(KIND=8) udata(*), uscale(*), fdata(*), fscale(*)
    real(KIND=8) vtemp1(*), vtemp2(*)
    integer,intent(out):: ier
    end subroutine fkpset
end interface

real(KIND=8), dimension(dimR*2), intent(inout) :: x1!, xg1

! Kinsol error flag !! Integer of 8 bits between -128 and 127
integer(KIND=1), intent(inout) :: ier 

real(KIND=8), dimension(dimR*2) :: sscale, constr
real(KIND=8) :: rout(2) ! Kinsol additional out information
real(KIND=8) :: fnormtol, scsteptol

integer(KIND=8), dimension(15) :: iout(15) ! Kinsol additional output information

! max. numb. of nonlinear iterations without a call to the preconditioner 
!setup function. integer of 8 bits [-128,127]
integer(KIND=1) :: msbpre 

integer :: i
integer :: globalstrat, maxl, maxlrst, maxiter

    print*, "call_kinsol: Variable neq desde call_kinsol: ", neq    
! INICIA KINSOL
! msbpre:  Cuantas iteraciones pasan antes de actualizar la matriz 
!          del precondicionador (operacion expensive operation) 

! maximum number of iterations without prec. setup (consultar kinsol_guide p.56) 
!    msbpre  = 5 ! (MARIO)
    msbpre  = 10 ! Default value: elefante

    maxiter = 300 ! Maximum number of non-linear iterations (Default=200)
!    maxiter = 200 ! Maximum number of non-linear iterations (Default=200)

    fnormtol = 1.0d-8 ! Function-norm stopping tolerance
    scsteptol = 1.0d-8 ! Scaled-step stopping tolerance
! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
!    maxl = 10 este estaba antes! elefante!
    maxl = neq ! Rikkert Comment
    maxlrst = 2 ! maximum number of restarts
    globalstrat = 0 ! this should be 1 for Inexact Newton or 2 for line search. elefante decia 0!
! Paso 2 - Inicializa the serial nvector module fnvinits(key,neq,ier); 
!            (key=3 for kinsol solver)
! neq=2*dimR  is declared in module_globales, is the size of the vectors x1 and xg1, defined in read_input.f90
! **********
! fnvinits inits NVECTOR module. 3 - for kinsol, neq equation number,
    call fnvinits(3, neq, ier) 
! EXISTE la version en paralelo fnvinitp
    if (ier .ne. 0) then       !  ier error flag (0 is OK)
        print*, 'SUNDIALS_ERROR: FNVINITS returned IER = ', ier
        stop
    endif
! Paso 3 -  Problem Specification
! iout - integer outputs
! rout - real outputs
    call fkinmalloc(iout, rout, ier) ! Allocates memory and output additional info
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
        stop
    endif

! maximum number of iterations without prec. setup 
    call fkinsetiin('MAX_SETUPS', msbpre, ier)  
! maximum number of Nonlinear iterations (Default=200)
    call fkinsetiin('MAX_NITER', maxiter, ier) 

    call fkinsetrin('FNORM_TOL', fnormtol, ier) ! Function-norm stopping tolerance
    call fkinsetrin('SSTEP_TOL', scsteptol, ier)! Scaled-step stopping tolerance
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSET*IN returned IER = ', ier
        stop
    endif

! Next definition of constraints
! Constraints (N Vector) vector of constraint flags. If constraints[i] is
!  0.0 then no constraint is imposed   on ui .
!  1.0 then u(i) will be constrained to be u(i) ≥ 0.0
! -1.0 then u(i) will be constrained to be u(i) ≤ 0.0
!  2.0 then u(i) will be constrained to be u(i) > 0.0
! -2.0 then u(i) will be constrained to be u(i) < 0.0
    do i = 1, neq  !constraint vector
! Indicate that no inequality constraints should be imposed on the solution vector.
        constr(i) = 0.0  
    enddo
    do i = 1, neq/2  !constraint vector
! The solvent volume fraction should be positive.
! this part of the solution vector (x1) will be constrained to be ui > 0.0.
        constr(i) = 2.0  
! PROBAR ESTA MODIFICACION PARA EVITAR DOS CICLOS!
! this part of the solution vector (x1) will be constrained to be ui > 0.0.
!        constr(neq/2 + i) = 2.0  
    enddo

    call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSETVIN returned IER = ', ier
        stop
    endif

!     CALL FKINSPTFQMR (MAXL, IER)
! The Generalized Minimal Residual Method (GMRES)

!  Scale Preconditioned GMRES solution of linear system (???)
!     call fkinspgmr(maxl, maxlrst, ier) 

! Biconjugate gradient stabilized method - (BCG)
! MAXL is the maximum Krylov subspace dimension.
    call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
        call fkinfree ! libera memoria
        stop
    endif

! Esta funcion indica a Kinsol que debe buscar las funciones fkpsol y fkpset 
! (precondiciones dadas por el usuario)
! fkinspilssetprec(flag,ier) flag=1.
    call fkinspilssetprec(1, ier) ! preconditiones
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSPILSSETPREC returned IER = ', ier
        call fkinfree ! libera memoria
        stop
    endif

    do i = 1, neq ! scaling vector
        sscale(i) = 1.0 ! all ==1 mean: no scaling is done
    enddo

!    do i = 1, neq ! Initial guess
!        x1(i) = xg1(i)
!    enddo

! Paso 6- Problem Solution
! El input es x1 y la Solucion vuelve en la variable x1
    call fkinsol(x1, globalstrat, sscale, sscale, ier)         ! Llama a kinsol
! U = x1 is an array containing the initial guess on input, 
!        and the solution on return. 
! GLOBALSTRAT is an integer (type INTEGER) defining the global strategy choice 
!             = 1 specifies Inexact Newton, 
!             = 2 indicates line search). 
! USCALE = sscale is an array of scaling factors for the U vector.
! FSCALE = sscale is an array of scaling factors for the FVAL vector. 
! IER is an integer completion flag and will have one of the following values: 
!           0 to indicate success, 
!           1 to indicate that the initial guess satisfies F (u) = 0 
!                         within tolerances, 
!           2 to indicate apparent stalling (small step), or a negative value 
!                         to indicate an error or failure.
    if (ier < 0) then
        print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
        print*, 'Linear Solver returned IER = ', iout(9)
        call fkinfree ! libera memoria
        stop
    endif

!    if (ier .lt. 0) then
!        print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
!        call fkinfree
!    endif

    return
endsubroutine call_kinsol
