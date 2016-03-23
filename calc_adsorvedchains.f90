!--------------------------------------------------------------------
!        Calcula Energia libre
! 2015-04-20 F.M.G. : Energía de Pong super chequeada!
!--------------------------------------------------------------------

subroutine calc_adsorvedchains(pHbulk)
#   include "control_run.h"
    use globales, only: radio, cuantas, long, dimR
    use pore, only: 
    use csys, only: pro, in1 , pR
    implicit none
    real(kind=8), intent(in) :: pHbulk
    real(kind=8) :: vect, radio_2
    integer :: nro_conf,j, nro_mon=0 
    integer, dimension(5) :: nro_cadenas(0:4)
    real(kind=8), dimension(5) :: p(0:4)

    print*, "Calculating adsorved chains..." 

    radio_2 = radio**2 
    nro_cadenas(:)=0
    p(:)=0.0
    do nro_conf=1,cuantas
        nro_mon=0
            do j=1,long
               ! vect=in1(nro_conf,j,1)**2 + in1(nro_conf,j,2)**2
               ! print*, radio_2, vect
               ! if( ( radio_2 -vect ) <= 1 ) then
               !     nro_mon=nro_mon+1
               ! endif
                if (pR(nro_conf,j) == dimR) nro_mon=nro_mon +1
            enddo
    !    if ( nro_mon == 0 ) then
     !       nro_cadenas(0)=nro_cadenas(0)+1
      !      p(0)=p(0)+pro(nro_conf)
        if ( nro_mon == 1 ) then
            nro_cadenas(1)=nro_cadenas(1)+1
            p(1)=p(1)+pro(nro_conf)
        else if( nro_mon == 2) then
            nro_cadenas(2)=nro_cadenas(2)+1
            p(2)=p(2)+pro(nro_conf)
        else if( nro_mon == 3) then
            nro_cadenas(3)=nro_cadenas(3)+1
            p(3)=p(3)+pro(nro_conf)
        else if( nro_mon > 3) then
            nro_cadenas(4)=nro_cadenas(4)+1
            p(4)=p(4)+pro(nro_conf)
        endif
    enddo

! Guarda numero de cadenas
         write(203,*) "# ",pHbulk
           p(0)=p(4) + p(3) + p(2) + p(1)
           nro_cadenas(0)=nro_cadenas(4) + nro_cadenas(3) + nro_cadenas(2) + nro_cadenas(1)
           
         do j=0,4
            write(203,*), j, p(j)*nro_cadenas(j), p(j), nro_cadenas(j)
         enddo
         write(203,*) " "
!c--------------------------- FIN -----------------
end subroutine calc_adsorvedchains
