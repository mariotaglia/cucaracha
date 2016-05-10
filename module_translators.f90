!
! Esta funcion translada el centro de masa de una dada conformacion del polimero 
! hasta algun lugar dentro del nanocanal. (Traslacion Pura)
!
! 1) Si el primer argumento es un entero debe corresponder al numero de layer donde se desea colocar la cadena
! 2) Si el primer argumento es la configuracion entonces se sortea un layer al azar. (OVERLOADING)
!
!***********************************************************************
module translators

interface transla
  module procedure transla_nlayer, transla_rand
end interface

contains

subroutine transla_nlayer(n_layer,xend,xendt,test)
#   include "control_run.h"
    use globales, only: radio, delta, long, dimR
    use csys
      implicit none
      REAL(KIND=8), intent(in) :: xend(3,long)
      REAL(KIND=8), intent(out) :: xendt(3,long)
      character*1, intent(out) :: test
      integer, intent(in) :: n_layer
      real(kind=8) :: vect
      integer :: i
      
! El layer donde cae el centro de masa debe ser aleatorio
           xendt(1, :) = xend(1,:) + (n_layer -0.5)*delta
           xendt(2, :) = xend(2,:) 
           xendt(3, :) = xend(3,:) 

! Testing configuration:
      do i=1,long ! loop sobre todos los segmentos
           vect = (xendt(1, i))**2 + (xendt(2, i))**2 ! distancia desde el centro al segmento ^ 2
           if ( (vect.gt.radio**2) ) then
                  test='N'  ! hay un segmento fuera del poro
                  cycle
           end if
      enddo
    return
end subroutine transla_nlayer

!!!!!!!! OVERLOADING FUNCTION TRANSLA
! Esta es la funcion que no necesita el layer donde poner la cadena
subroutine transla_rand(xend,xendt,test)
#   include "control_run.h"
    use globales, only: radio, delta, long, dimR
    use csys
      implicit none
      REAL(KIND=8), intent(in) :: xend(3,long)
      REAL(KIND=8), intent(out) :: xendt(3,long)
      character*1, intent(out) :: test
      real(kind=8) :: fac,vect
      real(kind=8) :: rands
      integer :: i
      
      fac=rands(seed)
      
      i = int( (dimR)*fac + 1)
! El layer donde cae el centro de masa debe ser aleatorio
           xendt(1, :) = xend(1,:) + (i -0.5)*delta
           xendt(2, :) = xend(2,:) 
           xendt(3, :) = xend(3,:) 

! Testing configuration:
      do i=1,long ! loop sobre todos los segmentos
           vect = (xendt(1, i))**2 + (xendt(2, i))**2 ! distancia desde el centro al segmento ^ 2
           if ( (vect.gt.radio**2) ) then
                  test='N'  ! hay un segmento fuera del poro
                  cycle
           end if
      enddo
      
    return
end subroutine transla_rand
end module translators
