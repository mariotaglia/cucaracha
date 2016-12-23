!
!
!

subroutine savetodisk(array, title, sigma, pHbulk, counter2)
    use globales
!    use csys, only: pHs
    implicit none
     
    real(KIND=8), dimension(dimR), intent(in) :: array
    character(len=5), intent(in) :: title 
    real(kind=8), intent(in) :: pHbulk, sigma
    integer, intent(in) :: counter2
    
    character(len=36) :: format_string, filename
    character(len=6) :: titlez
!    real(kind=8) :: posx, posy
    integer :: iR !, maxT, iT
    
! Un character(len=6) + character(len=1)+integer(len=3) with zeros on the left + character(len=1) ...
!    format_string='(A6, A1, I3.3, A1, I3.3, A4)'
    format_string='(A6,A1,F5.3,A1,I2.2,F0.2,A1,I3.3,A4)'

    titlez = title // 'z'
    ! Construye el array filename
    write(filename,format_string) titlez,'_', sigma ,'_', int(pHbulk) , pHbulk-int(pHbulk), '_',counter2, '.dat'
    ! Abre el archivo y guarda la informacion 
    open(unit=45, file=filename)
!        do iR=1, dimR
!            write(45,*)(-dimR + iR - 0.5)*delta,array(dimR-iR+1)
!        enddo

        do iR=1,dimR
            write(45,*)(iR-0.5+radio/delta)*delta,array(iR)
        enddo
    close(45)

!!     titlez = title // 't'
!! !   Construye el array filename
!! !    write(filename,format_string) titlez,'.', counter,'.', counter2, '.dat'
!!     write(filename,format_string) titlez,'_', int(pHbulk) , pHbulk-int(pHbulk), '_', counter2, '.dat'
!! !    write(filename,format_string) titlez,'_', pHs(counter2), '_', counter2, '.dat'
!!     open(unit=45, file=filename)
!!         maxT = 36
!!         do iT = 1, maxT
!!             do iR=1, dimR
!!                 posx = sin(dfloat(iT)/maxT*2.0*pi)*(iR)*delta
!!                 posy = cos(dfloat(iT)/maxT*2.0*pi)*(iR)*delta
!!                 write(45, *) posx, posy, array(iR)
!!             enddo
!!         enddo
!!     close(45)
!!     return
end subroutine savetodisk
