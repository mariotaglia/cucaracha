!
! Esta funcion acomoda la cadena dentro del nanocanal. 
!  Rotacion + Traslacion
!
!***********************************************************************

subroutine rota36(xend,xendr,n_in,test)
#   include "control_run.h"
! Se acomodan las cadenas cerca de la posicion
! (x,y) = (0,-radio)
    use globales, only: pi, radio, delta
    use csys
      implicit none
      integer, intent(in) :: n_in ! long = #segments
      REAL(KIND=8), intent(in) :: xend(3,n_in)
      REAL(KIND=8), intent(out) :: xendr(3,n_in)
      character*1, intent(out) :: test

      real*8 :: rands
      real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
      real*8 vect, dist_ymin, rmax, tmax, smax, cmax, xmax, ymax, rmaxaux, tmaxaux
      real*8 rmin, ymin, xmin
      real*8 alfa, gama, cga, a, b, c, xc, yc, zc
      integer :: i, int_ymin = 0
      
      fac=rands(seed)
      fac1=rands(seed)
      fac2=rands(seed)
      alfa=fac*2*pi
      cbe=2.0d0*fac1-1.0d0
      gama=fac2*2*pi

      sbe=(1-cbe**2)**0.5
      cal=cos(alfa)
      sal=sin(alfa)
      cga=cos(gama)
      sga=sin(gama)
      
! Preparo los valores iniciales para buscar el maximo
      !dist_ymin=xend(1,2)*(cbe*cal*sga+sal*cga)+ xend(2,2)*(cbe*cal*cga-sal*sga)-xend(3,2)*sbe*cal
      dist_ymin=xend(1,1)*(cbe*cal*sga+sal*cga)+ xend(2,1)*(cbe*cal*cga-sal*sga)-xend(3,1)*sbe*cal
# if CRITERIO ==1
      rmax=0.0
      xmax=0.0
      ymax=0.0
      int_ymin = 1
# elif CRITERIO ==2 
<<<<<<< HEAD
      !rmax = xend(1,1)
      rmax=xend(1,1)*(-cbe*sal*sga+cal*cga) -xend(2,1)*(cbe*sal*cga+cal*sga) +xend(3,1)*sbe*sal
#endif
      !int_ymin = 2
      int_ymin = 1
=======
      xmax = 0.0 !xend(1,1)
      int_ymin = 1
# elif CRITERIO ==3
! Calc Center of Mass
         xc=0.0
         yc=0.0
         zc=0.0
#endif
>>>>>>> release1.1

!      dist_ymin = 0.0 ! seguro existe al menos 1 segmento con y_segmento menor que este valor
      do 1 i=1,n_in
         a=xend(1,i)
         b=xend(2,i)
         c=xend(3,i)

         xendr(1,i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
         xendr(2,i)=a*(cbe*cal*sga+sal*cga)+ b*(cbe*cal*cga-sal*sga)-c*sbe*cal
         xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
# if CRITERIO == 1         
! Busco el rmax
         rmaxaux = sqrt( xendr(1,i)**2 + xendr(2,i)**2 ) 
         if ( rmaxaux > rmax ) then 
            ! here the polar coordinates of the farest monomer from the center of the nanochannel.
              int_ymin = i
              rmax = rmaxaux        
              tmax = atan(xendr(2,i)/xendr(1,i)) ! This is the angle
              xmax = xendr(1,i)
              ymax = xendr(2,i)
            !  cmax = xendr(1,i)/rmaxaux 
            !  smax = xendr(2,i)/rmaxaux 
            ! Este if es necesario por la definición de atan(x) is in ( -pi/2 , pi/2 )
!            if (xendr(1,i) < 0 ) tmax = tmax + pi
!              dist_ymin = xendr(2,i) 
         endif
#elif CRITERIO == 2
        rmaxaux = xendr(1,i) ! Maxima coordenada X
!        rmaxaux = sqrt( xendr(1,i)**2 + xendr(2,i)**2 ) 
        if ( rmaxaux > xmax ) then 
            xmax = rmaxaux     
            int_ymin = i   
            ymax = xendr(2,i) ! Coordenada Y correspondiente a la maxima coordenada X 
         endif
# elif CRITERIO == 3
! Calc Center of Mass
         xc=xc+xend(1,i)
         yc=yc+xend(2,i)
         zc=zc+xend(3,i)
# endif
 1    continue

# if CRITERIO == 3
! Pone el centro de masa en la posicion (0,0,0) de forma que
! la rutina rota36, solo hace una rotación (sin la traslación respecto del origen)
! Calc Center of Mass
         xc=xc/n_in
         yc=yc/n_in
         zc=zc/n_in
      do i=1,n_in
         xendr(1,i)=xendr(1,i)-( xc*(-cbe*sal*sga+cal*cga)-yc*(cbe*sal*cga+cal*sga)+zc*sbe*sal )
         xendr(2,i)=xendr(2,i)-( xc*(cbe*cal*sga+sal*cga)+ yc*(cbe*cal*cga-sal*sga)-zc*sbe*cal )
         xendr(3,i)=xendr(3,i)-( xc*sbe*sga+yc*sbe*cga+zc*cbe )
      enddo
# endif

        rmaxaux = 0.1
      do 2 i=1,n_in ! loop sobre todos los segmentos
# if CHAIN == 1 /*Here I choose the system: monolayers*/
! Si tengo monolayer entonces acerco el segmento mas cercano a la superfice.
! la idea es asegurar al menos un monomer adsorbido a a la superficie del canal. 
!*! Vieja forma:
!*!         xendr(2, i) = xendr(2,i) - radio + dist_ymin + 1e-5 ! cambia el origen del eje x, ahora el 0 esta en el centro del poro
#   if CRITERIO == 1
! Acerco el segmento con el rmax, manteniendo el angulo.
         xendr(1, i) = xendr(1,i) -xmax + (radio-0.01)*xmax /rmax  ! lleva la cordenada y del xmax a cero, ahora el 0 esta en el centro del poro
         xendr(2, i) = xendr(2,i) -ymax + (radio-0.01)*ymax /rmax  ! lleva la cordenada y del xmax a cero, ahora el 0 esta en el centro del poro
#   elif CRITERIO == 2 
         xendr(1, i) = xendr(1,i) - xmax + radio - 1e-5  ! lleva la cordenada y del xmax a cero, ahora el 0 esta en el centro del poro
         xendr(2, i) = xendr(2,i) - ymax ! lleva la cordenada y del xmax a cero, y=0 corresponde al eje que cruza al xmax
!         xendr(1, i) = xendr(1,i) + ( sqrt( radio**2 - dist_ymin **2) - rmax - 1e-5)  ! lleva la cordenada y del xmax a cero, ahora el 0 esta en el centro del poro
#   endif 

# elif CHAIN == 2 /*Here I choose the system: brushes*/
! Si brushes entonces acerco el primer segmento.
         xendr(1, i) = radio - xendr(1, i) + 1e-5 ! cambia el origen del eje y, ahora el 0 esta en el centro del poro
                                                  ! sumo 1e-5 para evitar que el primer segmento quede exactamento sobre la pared del poro
# endif /*Here I choose the system*/

         vect = (xendr(1, i))**2 + (xendr(2, i))**2 ! distancia desde el centro al segmento ^ 2
          
         !if ( (abs(dist_ymin) .gt. radio) .or. (vect.gt.radio**2) ) test='N'  ! hay un segmento fuera del poro
!         if ( (vect.ge.radio**2) ) test='N'  ! hay un segmento fuera del poro
         if ( (sqrt(vect) > radio) ) test='N'  ! hay un segmento fuera del poro
# ifdef fdebug_rota36
  print*,"rota36, i, (r^2 -vect), test: ", i, (radio**2 - vect), radio**2 - ( xendr(1,int_ymin)**2 + xendr(2,int_ymin)**2), test 
# endif
 2    continue

      return
end subroutine rota36
