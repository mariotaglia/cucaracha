! *****************************************************************************
! 
! Este es la version f90 del programa de nanocanales
! 
! *****************************************************************************

program nanochannel
!    include"control_run.h"
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
        subroutine save_data(ipH)
        integer, intent(in) :: ipH
        end subroutine

    end interface

    integer (KIND=1) :: ier = 0 ! kinsol error flag
    
    print*, " This is the monolayer program: "

! *****************************************************************************
! Estas subrutinas deberian estar declaradas en algun lugar en un contains o un 
! interface sin embargo funciona asi por que no reciben argumentos(creo) y por 
! que estan linkeadas con el Makefile.   
! *****************************************************************************
! Input data and units adaptation 
    call read_input 
    call units_adaptation ! units and variables adaptation

    call allocating(1) ! Allocating memory
    call creador ! Creating the chains
    call pxs ! Chequea que todos los segmentos esten dentro del slab.
!    call kai ! Calcula los parametros de L-J Xu(i,j) Solo para fvdW
    call mpinit(15) ! Initial working precision, number of digits
    call open_files(1) ! Open files to save data? how to do that?

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Falta una etapa donde se definir si se continua una corrida anterior
! En este caso habria una variable 'infile' por ahora es infile = 0.
! VER
!         if (ipH.eq.1) then 
!         endif
! IMPORTAN: DEFINE INITIAL GUESS -> SUBROUTINE!

    call set_bulk_properties(pHs(1)) ! Setup bulk properties before setup initial guess ocurre en fkfun! (?)
    call set_initial_guess(0) ! 0 - bulk solution
!    call set_pore_distrib Not necesarry this function is inside fkfun

    call printstate("Aloop L57") ! Report of State
! *****************************************************************************
! Principal Loop loop
! *****************************************************************************
    ipH=1 
    do while (ipH <= npH)  ! Start principal Loop
        print*, 'pH bulk =', pHs(ipH), 'ipH/npH:', ipH, '/', npH
!        print*, 'q', q
! Actualizo las condiciones de bulk se repite solo para ipH=1
        call set_bulk_properties(pHs(ipH))
! Resolution of the equations
        ier=0
!        print*, "Solucion inicial x1: ", x1(:)
        call call_kinsol(x1, ier) 

        xh(:) = x1(:dimR)    ! Solvent volume fraction
        psi(1:dimR) = x1(dimR+1:) ! Electrostatic Potential

! FUNCION: check_run == true if error.
!        if ( ipH == 4 .and. aux_del==0 ) aux_del = 1
!        if ( aux_del == 1 ) ier = -1
        if ( check_run(ipH,ier) ) then
!       Si hubo un error entonces imprimir alguna informacion y terminar
            print*, "funcion check run: Mal run!"
        ! Mean value between the last good value and the actual wrong value        
        !        goto 257 (abajo de 255(no calculapHbulk))
!            aux_del = 2
        else
!            print*, " funcion check run: " // "RUN OK!"
!           Si no exploto guardar xflag! (input para la proxima iteracion)
!            xflag(:) = x1(:) ! xflag sirve como input para la proxima iteracion
!           Preparar el siguiente paso: infile = 2
! Se escribe el output 
            call save_data(ipH) ! Saving data
            call calc_energy(pHs(ipH)) ! CALCULO DE ENERGIAS! 
            call calc_mean_values(pHs(ipH)) ! Rmedio 
! Calculo magnitudes derivadas: Gporo, Gneg, Gpos, fmedio, Rmedio,etc.
            call calc_conductance(pHs(ipH))
        endif

!        call printstate("main_L99") ! Report of State

        ipH=ipH+1
!        call set_pore_distrib() ! pore properties not necessary occur inside fkfun

   enddo  ! End principal Loop 

    call open_files(0) ! Closing all files

    !DeAllocating
    call allocating(0)
end program nanochannel
