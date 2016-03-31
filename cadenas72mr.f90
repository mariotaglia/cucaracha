!****************************************************************
! Esta función crea una cadena, utiliza el modelo 
! RIS - Rotational Isomeric Structure 
! luego rota36 la acomoda en el nanocanal.
!

subroutine cadenas72mr(chains,ncha)
#   include "control_run.h"
    use globales, only: lseg, long, pi,dimR, delta, radio
    use csys
    use translators
    implicit none
    integer :: ncha
!    REAL(KIND=8), intent(out) :: chains(3,200,100)
    REAL(KIND=8), intent(out) :: chains(3,200,130)
    
    real*8 rn,rands,state1,dista, vect
    real*8 sitheta,cotheta,siphip,cophip
    real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
    real*8 x(3),xend(3,200),xendr(3,200),xendt(3,200)

    integer i,state,ii,j,jj,ive,jve
    character*1 test

    sitheta=sin(68.0*pi/180.0)
    cotheta=cos(68.0*pi/180.0)
    siphip=sin(120.0*pi/180.0)
    cophip=cos(120.0*pi/180.0)
      
 223  x(1)=lseg
      x(2)=0.0
      x(3)=0.0
      
      xend(1,1)=lseg
      xend(2,1)=0.0
      xend(3,1)=0.0
      
      tt(1,1)=cotheta
      tt(1,2)=sitheta
      tt(1,3)=0.0
      tt(2,1)=sitheta
      tt(2,2)=-cotheta
      tt(2,3)=0.0
      tt(3,1)=0.0
      tt(3,2)=0.0
      tt(3,3)=-1.0
      
      tp(1,1)=cotheta
      tp(1,2)=sitheta
      tp(1,3)=0.0
      tp(2,1)=sitheta*cophip
      tp(2,2)=-cotheta*cophip
      tp(2,3)=siphip
      tp(3,1)=sitheta*siphip
      tp(3,2)=-cotheta*siphip
      tp(3,3)=-cophip
      
      tm(1,1)=cotheta
      tm(1,2)=sitheta
      tm(1,3)=0.0
      tm(2,1)=sitheta*cophip
      tm(2,2)=-cotheta*cophip
      tm(2,3)=-siphip
      tm(3,1)=-sitheta*siphip
      tm(3,2)=cotheta*siphip
      tm(3,3)=-cophip
      
 222  rn=rands(seed)
      
      state1=0.0
      
      m(1,1)=cotheta
      m(1,2)=sitheta
      m(1,3)=0.0
      
      m(2,1)=cos(state1)*sitheta
      m(2,2)=-cos(state1)*cotheta
      m(2,3)=sin(state1)
      m(3,1)=sin(state1)*sitheta
      m(3,2)=-sin(state1)*cotheta
      m(3,3)=-cos(state1)
      
      x(1)=m(1,1)*lseg
      x(2)=m(2,1)*lseg
      x(3)=m(3,1)*lseg
      
      xend(1,2)=lseg+x(1)
      xend(2,2)=x(2)
      xend(3,2)=x(3)
      
      do 10 i=3,long
         
! 123     rn=rands(seed)
         rn=rands(seed)
         state=int(rn*3)

         if (state.eq.3) then 
            state=2
         endif
         if (state.eq.0) then
            
            call mrrrr(m,tt,mm)
            do 30 ii=1,3
               do 40 j=1,3
                  m(ii,j)=mm(ii,j)
 40            continue
 30         continue
            
            
         elseif (state.eq.1) then
            
            call mrrrr(m,tp,mm)
            do 31 ii=1,3
               do 41 j=1,3
                  m(ii,j)=mm(ii,j)
 41            continue
 31         continue

         elseif (state.eq.2) then
            
            call mrrrr(m,tm,mm)
            do 32 ii=1,3
               do 42 j=1,3
                  m(ii,j)=mm(ii,j)
 42            continue
 32         continue
            
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
         
 10   continue
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chequea cadenas
! Selfavoiding entre segmentos
 
      dista=0.0
      do 300 ive=4,long
         do 310 jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
!cccccccccwrite(2,*)'noself-a.',k
               goto 222
            endif
 310     continue
 300  continue

      ncha=0
      do 400 i=1,300
         test='S'
         call rota36(xend,xendr,long,test)
         if (test.eq.'N') goto 400
# if CRITERIO == 3 
!Criterio 3 Le agrega un grado de libertad a las configuraciones de polimero
!         do 403 jj=1,3000 ! Mueve por layers al azar.
!            test='S'
!            call transla(jj,xendr,test)
!            if (test.eq.'N') go to 403 ! Mueve por layers tambien al azar.

          do jj=1,dimR ! Mueve por layers
            test='S'
            call transla(jj,xendr, xendt,test)
            if (test.eq.'N') cycle ! Mueve por layers
# endif
            ncha=ncha+1

# if CRITERIO == 3 
            do j=1,long
               chains(1,j,ncha)=xendt(1,j) ! y
               chains(2,j,ncha)=xendt(2,j) ! x
               chains(3,j,ncha)=xendt(3,j) ! z
            enddo
         !if (ncha.eq.25) goto 402
         if (ncha.eq.125) goto 402
# ifdef fdebug
         print*, " WARNING! (L180, cadenas72) Esto supone un sesgo en las configuraciones: if (ncha.eq.25) goto 402 "
         print*, " WARNING! (L181, cadenas72) Para radios grandes hay que pensar otra cosa."
# endif
!403     continue
        enddo ! 403     continue
# else
            do j=1,long
               chains(1,j,ncha)=xendr(1,j) ! y
               chains(2,j,ncha)=xendr(2,j) ! x
               chains(3,j,ncha)=xendr(3,j) ! z
            enddo
         if (ncha.eq.25) goto 402
# endif
 400  continue
 402  if (ncha.eq.0) goto 223

      return
end subroutine cadenas72mr
