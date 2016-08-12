! *****************************************************************************
! 
! Este es la version f90 del programa de nanocanales
! version actual: 1.0.1 
! *****************************************************************************

program nanochannel
#   include "control_run.h"
    use globales
    use csys
    use pore 
    
    implicit none
    interface
        function check_run(ipH,ier)
        logical :: check_run
        integer, intent(inout) :: ipH
        integer(KIND=1), intent(in) :: ier
        end function 
        subroutine save_data(ipH,isigma)
        integer, intent(in) :: ipH
        end subroutine

    end interface

    integer (KIND=1) :: ier = 0 ! kinsol error flag
    integer :: i 
    print*, " This is the monolayer program: "
    print*, 'GIT version = ', _VERSION 

! *****************************************************************************
! Estas subrutinas deberian estar declaradas en algun lugar en un contains o un 
! interface sin embargo funciona asi por que no reciben argumentos(creo) y por 
! que estan linkeadas con el Makefile.   
! *****************************************************************************
! Input data and units adaptation 
    call read_input 
    call units_adaptation ! units and variables adaptation

    call allocating(1) ! Allocating memory

! If chain ==1 then prepare the memory and variables
#if CHAIN !=0   
    print*, "Entro a creador!"
    call creador ! Creating the chains
    print*, "Salgo de creador!"
    call pxs ! Chequea que todos los segmentos esten dentro del slab.
#ifdef  VDW
    call kai ! Calcula los parametros de L-J 
#endif
#endif
    call mpinit(15) ! Initial working precision, number of digits =15
    

!    write(11,*) "set_pore_distrib 109: q "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Falta una etapa donde se definir si se continua una corrida anterior
! En este caso habria una variable 'infile' por ahora es infile = 0.
! VER
!         if (ipH.eq.1) then 
!         endif
! IMPORTANT: DEFINE INITIAL GUESS -> SUBROUTINE!
    ipH=1;icsalt=1;icpol=1
! Setup bulk properties before setup initial guess ocurre en fkfun! (?)
! bulk properties without polymer
    print*, "vcpol(icpol): ", vcpol(icpol)
    call set_bulk_properties(pHs(ipH),vcsalt(icsalt),vcpol(icpol)) 
! Inside nanochannel set x1
    call set_initial_guess(0) ! 0 - bulk solution as initial guess
    
    call printstate("Aloop L69") ! Report of State
    
    print*, '********************************************'
    write(0,'(23a)'), '*** Comienza el loop principal ***' ! stderr
! *****************************************************************************
! Principal Loop loop
! *****************************************************************************
    ipH=1 
    do while (ipH <= npH)  ! Start principal Loop
     ! Actualizo las condiciones de bulk se repite solo para ipH=1
        print*, 'pH bulk =', pHs(ipH), 'ipH/npH:', ipH, '/', npH
        call open_files(1,ipH) ! Open files to save data? how to do that?
     !   call set_bulk_properties(pHs(ipH)) 

!        call set_pore_distrib !Not necesarry this function is inside fkfun

! Resolution of the equations
        ier=0
        !print*, "Solucion inicial x1: ", x1(:)
        icpol = 1
        do while (icpol <= ncpol)  ! Start principal Loop
            cpol=vcpol(icpol)
            print*, 'Cpol bulk =', vcpol(icpol), 'ipH/npH:', icpol, '/', ncpol
            call set_bulk_properties(pHs(ipH),vcsalt(icsalt),vcpol(icpol)) 
# ifdef fdebug
            call printstate("fdebug Aloop L94") ! Report of State
# endif
            call call_kinsol(x1, ier)
             
            write(10,*), " " ! para formatear fort.10
            write(11,*), " " ! para formatear fort.11
            xh(:) = x1(:dimR)    ! Solvent volume fraction
            psi(1:dimR) = x1(dimR+1:) ! Electrostatic Potential
#ifdef fdebug
        print*, 'monolayer L91,  xh(:)= ' , xh(:)
        print*, 'monolayer L91, psi(:)= ' , psi(:)
#endif 

! FUNCION: check_run == true if error.
        if ( check_run(ipH,ier) ) then
!       Si hubo un error entonces imprimir alguna informacion y terminar
            print*, "funcion check run: Mal run!"
            print*, "cpol: ", vcpol(icpol)
            stop
        ! Mean value between the last good value and the actual wrong value        
        !        goto 257 (abajo de 255(no calculapHbulk))
!            aux_del = 2
        else
!            print*, " funcion check run: " // "RUN OK!"
!  Si no exploto guardar xflag! (input para la proxima iteracion)
            if ( icpol==1  ) then
! Esto es solo para la primer iteracion por que se supone 
! esta solución es "cercana" a la proxima iteracion de ph
                select case (2)
                    case ( 0 )  ! Bulk solution
                    case ( 1 )  ! Save to file 
                         open(unit=92,file='initial_guess')
                            write(92,*) "# sigma= ", sigma*delta/vsol, " pH= ", pHs(ipH)
                            write(92,*) ! ignore first and second lines
                        
                            do i=1,2*dimR
                                write(92,*) x1(i)
                            end do
                         close(92)
                    case ( 2 )  ! Last solution
                            xflag(:) = x1(:) ! xflag sirve como input para la proxima iteracion
                end select
            endif
! Se escribe el output 
            call save_data(ipH,icpol) ! Saving data
            call calc_energy(pHs(ipH)) ! CALCULO DE ENERGIAS!
            call calc_mean_values(pHs(ipH)) ! Rmedio
!            call calc_adsorvedchains(pHs(ipH)) !Nro de cadenas adsorvidas en el poro
!            call calc_pkas() 
! Calculo magnitudes derivadas: Gporo, Gneg, Gpos, fmedio, Rmedio,etc.
!            call calc_conductance(pHs(ipH))
        endif
        icpol= icpol +1
        enddo ! loop over sigma
!       Inside nanochannel set x1
        call set_initial_guess(2) ! 2 - read from xflag
        ipH=ipH+1
   enddo  ! End principal Loop 

    call open_files(0) ! Closing all files
    call allocating(0) !DeAllocating
    write(0,'(23a)'), '*** Fin loop principal ***' ! stderr
end program nanochannel
