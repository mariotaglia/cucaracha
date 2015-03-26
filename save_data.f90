!
!
!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine save_data(ipH)
#   include "control_run.h"
    use globales, only: delta, dimR
    use csys, only: pHs!,x1
    use pore, only: avpol, qtot, xh, xpos, xneg, xHplus, xOHmin, fdis, psi!, fdis!, fdis2
    integer, intent(in) :: ipH
    real(kind=8) :: pH
    character(len=5) :: title
!    character(len=17) :: outfile
!    character(len=30) :: format_string, filename

!    integer :: i
    
!    write(outfile,'(A14,I3.3)') 'initial_guess_',ipH
!        open(unit=92,file=outfile)
!            write(92,*) '# This lines are omitted'
!            write(92,*) '# Printed with save_data function'
!            do i=1,2*dimR
!                write(92,*) x1(i)
!            enddo
!        close(92)

    print*, "Saving data..."
! GUARDAR! usar savetodisk
        pH=pHs(ipH)
#if CHAIN == 1
! Polimero
      title = 'avpol'
      call savetodisk(avpol, title, pH ,ipH)
! fdis
      title = 'fdis1'
      call savetodisk(fdis, title, pH, ipH)

!! ! fdis2
!!       title = 'fdis2'
!!       call savetodisk(fdis2, title, pH, ipH)

#endif
! Polimero
      title = 'qtodo'
      call savetodisk(qtot, title, pH ,ipH)

! Solvente
      title = 'avsol'
      call savetodisk(xh, title, pH, ipH)

! Cationes
      title = 'avpos'
      call savetodisk(xpos, title, pH, ipH)

! Cationes2
!      title = 'avpo2'
!      call savetodisk(xpos2, title, pH, ipH)

! Aniones
      title = 'avneg'
      call savetodisk(xneg, title, pH, ipH)

! H+
      title = 'avHpl'
      call savetodisk(xHplus, title, pH, ipH)

! OH-
      title = 'avOHm'
      call savetodisk(xOHmin, title, pH,ipH)

! Potencial electrostatico
!        call printstate('save_data L71')
      title = 'poten'
      call savetodisk(psi, title, pH, ipH)

!!     format_string='(A5,A1,I2.2,F0.2,A1,I3.3,A4)'
!!     ! Construye el array filename
!!     write(filename,format_string) title,'_', int(pHs(ipH)) , pHs(ipH)-int(pHs(ipH)), '_', ipH, '.dat'
!!     ! Abre el archivo y guarda la informacion 
!!     open(unit=45, file=filename)
!! !        do iR=1, dimR
!! !            write(45,*)(-dimR + iR - 0.5)*delta,array(dimR-iR+1)
!! !        enddo
!!         do i=0,dimR+1
!!             write(45,*)(i-0.5)*delta,psi(i)
!!         enddo
!!     close(45)


end subroutine save_data

