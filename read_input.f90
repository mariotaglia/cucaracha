subroutine read_input
    use globales, only: cuantas, dimR, long, radio, delta, neq
    use csys
    implicit none
    integer :: i
    character :: basura
!    print*, "Numero de avogadro", Na
    open(8,status='old', action='read' )

    print*, 'Reading input parameters'

! Input variables from fort.8

      read(8, *), basura
      read(8, *), sigma   ! Surface coverage

      read(8, *), basura
      read(8, *), sigmaq   ! Surface chargeable sites

      read(8, *), basura
      read(8, *), eps1     ! Interaccion Pol-Sup 
! Salt concentration
      read(8, *), basura
      read(8, *), csalt   ! Bulk Salt concentration

      read(8, *), basura
      read(8, *), pKawall ! Eq. Const. wall charge

      read(8, *), basura
      read(8, *), npKa
      read(8, *), pKa, pKb    ! Eq. Const. Polymer
     
      read(8, *), basura
      read(8, *), pKb     ! Eq. Const. Polymer
 
      read(8, *), basura
      read(8, *), npH
      allocate(pHs(npH))  ! OJO! esta memoria debe ser liberada en el main! becareful!
      read(8, *), (pHs(i), i=1, npH)

      read(8, *), basura
      read(8, *), infile
! st parameter is for vdW interaction on/off ?
      read(8, *), basura
!      read(8, *), st
      read(8, *), cuantas, dimR, long     ! radio size
        radio=dimR*delta
        neq=2*dimR

      read(8, *), basura
      read(8, *), movpos, movneg, movHplus, movOHmin, longporo
    close(unit=8)

end subroutine 
