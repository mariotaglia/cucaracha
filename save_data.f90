!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine save_data(ipH,icpol)
#   include "control_run.h"
    use globales, only: delta, dimR, vsalt, vsol
    use csys, only: pHs, expmuHplus, expmuOHmin, expmupos, expmuneg, vcpol, sigma,csalt!,x1
    use pore, only: avpol, qtot, xh, xpos, xneg, xHplus, xOHmin, fdis, psi!, fdis!, fdis2
    integer, intent(in) :: ipH, icpol
    real(kind=8) :: pH, cpol
    character(len=5) :: title
!    character(len=17) :: outfile
    character(len=30) :: format_string, filename

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
        cpol=vcpol(icpol)
!    print*, 'SIGMA: ', pH, cpol, csalt, dimR, sigma
    write(325,*) cpol, sigma
# if CHAIN != 0
! Polimero
      title = 'avpol'
      call savetodisk(avpol, title, sigma, pH ,ipH)
! fdis
      title = 'fdis1'
      call savetodisk(fdis, title, sigma, pH, ipH)

!! ! fdis2
!!       title = 'fdis2'
!!       call savetodisk(fdis2, title, sigma, pH, ipH)

#endif /* CHAIN != 0 */
! Total Charge
      title = 'qtodo'
      call savetodisk(qtot, title, sigma, pH ,ipH)

! Solvente
      title = 'avsol'
      call savetodisk(xh, title, sigma, pH, ipH)

! Cationes
      title = 'avpos'
      call savetodisk(xpos, title, sigma, pH, ipH)

! Cationes2
!      title = 'avpo2'
!      call savetodisk(xpos2, title, sigma, pH, ipH)

! Aniones
      title = 'avneg'
      call savetodisk(xneg, title, sigma, pH, ipH)

! H+
      title = 'avHpl'
      call savetodisk(xHplus, title, sigma, pH, ipH)

! OH-
      title = 'avOHm'
      call savetodisk(xOHmin, title, sigma, pH,ipH)

! Potencial electrostatico
!        call printstate('save_data L71')
      title = 'poten'
      call savetodisk(psi, title, sigma, pH, ipH)

# ifdef debug
    format_string = '(A5,A1,I2.2,F0.2,A1,I3.3,A4)'
    ! Construye el array filename
    title='b_mus'
    !write(filename,format_string) title,'_', int(pHs(ipH)) , pHs(ipH)-int(pHs(ipH)), '_', ipH, '.dat'
    write(filename,format_string) title,'_', int(pH) , pH-int(pH), '_', ipH, '.dat'
    ! Abre el archivo y guarda la informacion 
    open(unit=45, file=filename)
        do iR=1, dimR
            write(45,*) (iR - 0.5)*delta, dlog(expmupos), dlog(expmuneg), dlog(expmuHplus), dlog(expmuOHmin) !array(dimR-iR+1)
        enddo
!        do i=0,dimR+1
!            write(45,*)(i-0.5)*delta,psi(i)
!        enddo
    close(45)
!******************************************************************
    format_string = '(A5,A1,I2.2,F0.2,A1,I3.3,A4)'
    ! Construye el array filename
    title='p_mus' ! en el poro
    !write(filename,format_string) title,'_', int(pHs(ipH)) , pHs(ipH)-int(pHs(ipH)), '_', ipH, '.dat'
    write(filename,format_string) title,'_', int(pH) , pH-int(pH), '_', ipH, '.dat'
    ! Abre el archivo y guarda la informacion 
    open(unit=45, file=filename)
        do iR=1, dimR
            write(45,*) (iR - 0.5)*delta, dlog( xpos(iR)/ vsalt / (xh(iR)**vsalt) ) + psi(iR), &!array(dimR-iR+1)
                                          dlog( xneg(iR)/ vsalt / (xh(iR)**vsalt) ) - psi(iR), &
                                          dlog( xHplus(iR) / xh(iR) ) + psi(iR), &
                                          dlog( xOHmin(iR) / xh(iR) ) - psi(iR)
        enddo
!        do i=0,dimR+1
!            write(45,*)(i-0.5)*delta,psi(i)
!        enddo
!******************************************************************
    close(45)
#endif 
end subroutine save_data
