function check_run(j_pHbulk,ier)
    use globales, only: error
    use csys
    implicit none
    logical :: check_run
    integer, intent(inout) :: j_pHbulk
    integer(kind=1), intent(in) :: ier
! Chequea si exploto... => Sistema anti-crash
!    logical :: check_run
    real(kind=8) :: mean_pH!, pHbulk
        
!    pHbulk=pHs(j_pHbulk)
!    if( (ier < 0) .or. (norma > error) ) then ! exploto...
    if( (ier < 0) ) then ! exploto...
        call printstate('Check Run BUM!')
        print*, 'Check Run says: BUM! see: printstate.txt '
        print*, 'Error en solver: ', ier
!        print*, 'norma ', norma
!        print*, 'q ', q
        print*, 'st de error', st
        print*, 'pH de error', pHs(j_pHbulk)
       
!        mean_pH = (pHbulk + pHs(j_pHbulk-1))/2
        if ( j_pHbulk ==1 ) then 
            print*, "*********************************"
            print*, "Falla para el primer valor de pH!"
            print*, "*********************************"
            stop
        endif

        mean_pH = (pHs(j_pHbulk) + pHs(j_pHbulk-1))/2
        write(1010,*) 'Fallo pH ', pHs(j_pHbulk), ' Paso a ', mean_pH
        pHs(j_pHbulk) = mean_pH
        flag = 1
        j_pHbulk=j_pHbulk-1 ! one step backward. 
        
        check_run = .true.
    else
        infile = 2 ! no vuelve a leer infile
        check_run = .false.
    endif

    if(flag.eq.1) then  ! habia un error...
        print*, 'Recupero del error'
        print*, 'OK', pHs(j_pHbulk)
        flag = 0
    endif

end function check_run
