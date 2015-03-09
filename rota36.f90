!
! Esta funcion acomoda la cadena dentro del nanocanal. 
!  Rotacion + Traslacion
!
!***********************************************************************

subroutine rota36(xend,xendr,n_in,test)
! Se acomodan las cadenas cerca de la posicion
! (x,y) = (0,-radio)
    use globales, only: pi, radio
    use csys
      implicit none
      integer, intent(in) :: n_in ! long = #segments
!      REAL(KIND=8), dimension(3,200), intent(in) :: xend 
!      REAL(KIND=8), dimension(3,200), intent(out) :: xendr
      REAL(KIND=8), intent(in) :: xend(3,n_in)
      REAL(KIND=8), intent(out) :: xendr(3,n_in)
      character*1, intent(out) :: test

      real*8 :: rands
      real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
      real*8 vect, dist_ymin
      real*8 alfa, gama, cga, a, b, c
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
      dist_ymin=xend(1,2)*(cbe*cal*sga+sal*cga)+ xend(2,2)*(cbe*cal*cga-sal*sga)-xend(3,2)*sbe*cal
      int_ymin = 2
!      dist_ymin = 0.0 ! seguro existe al menos 1 segmento con y_segmento menor que este valor
      do 1 i=1,n_in

         a=xend(1,i)
         b=xend(2,i)
         c=xend(3,i)

         xendr(1,i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
         xendr(2,i)=a*(cbe*cal*sga+sal*cga)+ b*(cbe*cal*cga-sal*sga)-c*sbe*cal
         xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
         
! Busco el xmax
         if ( dist_ymin > xendr(2,i) ) then 
              int_ymin = i 
              dist_ymin = xendr(2,i) 
         endif
 
 1    continue

! Chequea que la cadena rotada quepa dentro del poro
!
! (x seg - x radio)^2 + (y seg - y radio)^2 < radio^2
!
!      dist_ymax = xendr(2,int_xmax)

      do 2 i=1,n_in ! loop sobre todos los segmentos
! Si tengo monolayer entonces acerco el segmento mas cercano a la superfice.
! la idea es asegurar al menos un monomer adsorbido a a la superficie del canal. 
         xendr(2, i) = xendr(2,i) - radio + dist_ymin + 1e-5 ! cambia el origen del eje x, ahora el 0 esta en el centro del poro
!         xendr(2, i) = xendr(1,i) - dist_ymax ! lleva la cordenada y del xmax a cero, ahora el 0 esta en el centro del poro

! Si brushes entonces acerco el primer segmento.
!         xendr(1, i) = radio - xendr(1, i) + 1e-5 ! cambia el origen del eje y, ahora el 0 esta en el centro del poro
                                                  ! sumo 1e-5 para evitar que el primer segmento quede exactamento sobre la pared del poro
 
         vect = (xendr(1, i))**2 + (xendr(2, i))**2 ! distancia desde el centro al segmento ^ 2
         if (vect.gt.radio**2) test='N'  ! hay un segmento fuera del poro
                  
 2    continue

      return
end subroutine rota36
