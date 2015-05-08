!#####################################################################
! Este programa calcula los kai (parametro de lennard-jones)
! para poor-solvent en 1D-coordenadas polares usando un metodo de MC
!
!#####################################################################

subroutine kai
#   include "control_run.h"
    use globales
    use csys, only: Xu
    implicit none

    real(kind=8), dimension(:,:), allocatable :: suma
    real(kind=8) :: x_1,x_2,y_1, y_2, z1, z2, vect
    real(kind=8) ::  R,  z,  theta
    integer :: ii,j,iR, iz, itheta
    integer :: MCstepsR, MCstepsz, MCstepstheta ! numero de steps de MC

    print*,'Kai calculation'
    allocate( suma(dimR,dimR) )
    open(unit=111, file='kais.dat')
!    open(unit=811, file='aux.dat')

    suma(:,:) = 0.0
    Xu = 0.0 ! vector Xu
    MCstepsz = 100
    MCstepsR = dimR*100
      ! NOTA: Dar vuelta los indices ii,j en Xu y suma, porque fortran guarda en memoria las cosas al reves! 
      !*********************************
    do ii = 1, dimR ! loop sobre cada posicion del segmento
!    ii=dimR
        x_1 = (dfloat(ii) - 0.5)*delta ! asume theta segmento = 0, z segmento = 0 y segmento en el centro de la layer
        y_1 = 0.0
        z1 = 0.0

!        write(811,*) '      ii       iz      z       itheta      theta: '
        MCstepstheta = ii*100
        do itheta = 1, MCstepstheta 
!            itheta = 1
            do iz = 1, MCstepsz + 1
!                iz = 1
!                do iR = 1, MCstepsR + 1 ! Esta es la linea que genera error PREGUNTAR A MARIO!!
                do iR = 1, MCstepsR  ! Con esta linea funciona! El tema es que j vale dimR+1 cuando iR=MCstepsR
                   ! iR = 1
                    theta = 2*pi*dfloat(itheta - 1)/dfloat(MCstepstheta)    ! numero 0 y 2pi
                    z = 3.0*(dfloat(iz - 1)/dfloat(MCstepsz)-0.5)*delta ! numero entre -1.5*delta y 1.5*delta
                    R = radio*(dfloat(iR-1)/dfloat(MCstepsR))         ! numero 0 y radio
! coordenadas del segmento (x_1,y_1,z1) y del punto a integrar (x_2,y_2,z2)
      !             x_1 = (dfloat(ii) - 0.5)*delta ! asume theta segmento = 0, z segmento = 0 y segmento en el centro de la layer
      !             y_1 = 0.0
      !             z1 = 0.0
!            write(811,*) 'inicio', ii, iz, z, itheta, theta 

                    x_2 = R*cos(theta)
                    y_2 = R*sin(theta)
                    z2 = z
                    vect = sqrt((x_1-x_2)**2 + (y_1-y_2)**2 + (z1-z2)**2) ! vector diferencia
!********************************
      !OJO! R=radio=dimR*delta --> j=21. por que +1?
!********************************
                    j = int(R/delta)+1 ! j tiene la celda donde cae el punto a integrar
!print*, j, iR, R, iz, z, itheta, theta
                    suma(ii, j) = suma(ii, j) + R
                    !write(811,*) 'iz, z, itheta, theta: ', iz, z, itheta, theta, ' coordenada iR, R:', iR, R!, "suma: ", suma(ii,j) , "Xu: ", Xu(ii,j)
! This conditional exclude interactions out of range
                    if( (vect.gt.(1.5*delta)).or.(vect.lt.lseg)) cycle

!!!!!!!! Old form of the conditional
!                    if(vect.gt.(1.5*delta)) then
!                        print*, "vect > 1.5*delta",  ii, iz
!                        cycle
!                    endif
!                    if(vect.lt.lseg) then
!                        print*, "vect < lseg",  ii, iz
!                        cycle
!                    endif
        !                write(811,*) 'ii, iR, R', ii, iR, R
                    Xu(ii, j) = Xu(ii, j) + ((lseg/vect)**6)*R ! incluye el jacobiano R(segmento)

                   enddo ! iz
            enddo ! itheta
        enddo ! ii

        do j = 1, dimR
            Xu(ii, j) = Xu(ii, j)/(MCstepsR*MCstepsz*MCstepstheta)*(3.0*delta)*2*pi*radio
            suma(ii, j) = suma(ii, j)/(MCstepsR*MCstepsz*MCstepstheta)*(3.0*delta)*2*pi*(radio)
            write(111,*) ii,j,Xu(ii,j) ! residual size of iteration vector
       !     print*, 'suma,i,j', suma(ii,j), ii, j
        enddo
       !         write(811,*) 'print 9'

    end do ! ii

    close(111)
!    close(811)
    deallocate( suma )
    print*, 'kai, OK!'
end subroutine kai
