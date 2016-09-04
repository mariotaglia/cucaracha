!****************************************************************
! Generates a RIS chain inside a nanpore
! (RIS - Rotational Isomeric Structure)
! using Rosenbluth sampling to favor conformations next to the surface
! 
!
!
!
!
!

subroutine cadenas72mr_rosen(chains,ncha, wchains)
    use globales, only: lseg, long, pi,dimR, delta, radio
    use csys
    use translators
    implicit none
    integer :: ncha
!    REAL(KIND=8), intent(out) :: chains(3,200,100)
    REAL(KIND=8), intent(out) :: chains(3,200,130)
    REAL(KIND=8), intent(out) :: wchains(200)
    
    real*8 rn,rands,state1,dista!, vect
    real*8 sumexpu
    real*8 sitheta,cotheta,siphip,cophip
    real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
    real*8 x(3),xend(3,200),xendt(3,200)
    integer i,state,ii,j,ive,jve
    real*8 expu(3), wrosen, utot
    real*8 tmp
    real*8 expenergy

    sitheta=sin(68.0*pi/180.0)
    cotheta=cos(68.0*pi/180.0)
    siphip=sin(120.0*pi/180.0)
    cophip=cos(120.0*pi/180.0)
 
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
      
! m is the initial rotation matrix


! first segment 
 222   xend(1,1)=radio-delta/2 ! first segment is on the wall
       xend(2,1)=0.0
       xend(3,1)=0.0

      utot = -log(expenergy(xend,1)) ! exp(-energy) of a segment a position xpot
      wrosen = expenergy(xend,1) ! rosenbluth weight up to first segment    

      call randommatrix(m) ! random direction

      do i = 2,long ! loop over segments
      do state = 1,3 ! three RIS positions
      select case (state)
        case(1)
         call mrrrr(m,tt,mm)
        case(2)
         call mrrrr(m,tp,mm)
        case(3)
         call mrrrr(m,tp,mm)
       endselect ! mm contains the rotated matrix

         x(1)=mm(1,1)*lseg
         x(2)=mm(2,1)*lseg
         x(3)=mm(3,1)*lseg
 
         xendt = xend
         xendt(:,i) = xendt(:,i-1) + x(:) ! trial position 

         expu(state) = expenergy(xendt, i) ! boltzmann factor of segment(state), energy must consider repulsions with other segments and wall
       enddo ! state
     
       sumexpu = sum(expu) !   
       if(sumexpu.eq.0.0) then  ! dead end... we should differentiate here if all possible choinces are drop due to collisions with the wall, because in that case the chain should be counted in bulk
       goto 222
       endif
    

       expu = expu/sumexpu ! probability of selecting segment state from 0 to 1
       rn = rands(seed)
       tmp = 0.0
       do state = 1, 3
         tmp = tmp + expu(state)
         if(tmp.gt.rn)exit ! pick state
       enddo
        
        select case (state)
        case(1)
         call mrrrr(m,tt,mm)
        case(2)
         call mrrrr(m,tp,mm)
        case(3)
         call mrrrr(m,tp,mm)
       endselect ! mm contains the rotated matrix

       m=mm

       x(1)=m(1,1)*lseg
       x(2)=m(2,1)*lseg
       x(3)=m(3,1)*lseg
  
       xend(1,i)=xend(1,i-1)+x(1)
       xend(2,i)=xend(2,i-1)+x(2) 
       xend(3,i)=xend(3,i-1)+x(3)

      utot =  utot -log(expenergy(xend,i)) ! exp(-energy) of a segment a position xpot
      wrosen = wrosen * sumexpu/3.0 ! rosenbluth weight up to first segment   
     
 enddo ! i

      ncha=1
            do j=1,long
               chains(1,j,ncha)=xend(1,j) ! y
               chains(2,j,ncha)=xend(2,j) ! x
               chains(3,j,ncha)=xend(3,j) ! z
            enddo

      wchains(1) = 3.0**float(long)*exp(-utot)/wrosen
return
end subroutine 


function expenergy(xend,pos)
use globales, only: delta, radio, lseg, eps_rosen
implicit none
real*8 expenergy
real*8 dista
real*8 xend(3,200)
integer pos, ive, jve

expenergy = 1.0 
dista = sqrt(xend(1,pos)**2 + xend(2,pos)**2)

!if(dista.gt.radio)expenergy = 0.0 ! out-of-pore
if((dista.lt.radio).and.(dista.gt.(radio-delta)))expenergy = exp(-eps_rosen)

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chequea cadenas
! Selfavoiding entre segmentos
 
      dista=0.0
      do ive=4,pos
         do jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
               expenergy = 0.0
            endif
        enddo
        enddo

end function




subroutine randommatrix(m)
      use csys
      use globales
      implicit none
      real*8 m(3,3)
      real*8 :: rands
      real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
      real*8 vect, dist_ymin, rmax, tmax, smax, cmax, xmax, ymax, rmaxaux, tmaxaux
      real*8 rmin, ymin, xmin
      real*8 alfa, gama, cga, a, b, c, xc, yc, zc
      
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
!      dist_ymin=xend(1,2)*(cbe*cal*sga+sal*cga)+ xend(2,2)*(cbe*cal*cga-sal*sga)-xend(3,2)*sbe*cal

      m(1,1) = cbe*cal*sga
      m(2,1) = sal*cga
      m(3,1) = 0.0
      m(1,2) = cbe*cal*cga
      m(2,2) = -sal*sga
      m(2,3) = 0.0
      m(1,3) = 0.0
      m(2,3) = 0.0
      m(3,3) = -sbe*cal

      return
end subroutine

