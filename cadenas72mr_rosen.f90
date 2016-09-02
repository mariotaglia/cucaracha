!****************************************************************
! Esta función crea una cadena, utiliza el modelo 
! RIS - Rotational Isomeric Structure 
! luego rota36 la acomoda en el nanocanal.
!

subroutine cadenas72mr_rosen(chains,ncha)
    use globales, only: lseg, long, pi,dimR, delta!, radio
    use csys
    use translators
    implicit none
    integer :: ncha
!    REAL(KIND=8), intent(out) :: chains(3,200,100)
    REAL(KIND=8), intent(out) :: chains(3,200,130)
    
    real*8 rn,rands,state1,dista!, vect
    real*8 sitheta,cotheta,siphip,cophip
    real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
    real*8 x(3),xend(3,200),xendr(3,200)
    integer i,state,ii,j,ive,jve
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
      
      do i=3,long
         
! 123     rn=rands(seed)
         rn=rands(seed)
         state=int(rn*3)

         if (state.eq.3) then 
            state=2
         endif
         if (state.eq.0) then
            
            call mrrrr(m,tt,mm)
            m=mm
            
         elseif (state.eq.1) then
            
            call mrrrr(m,tp,mm)
            m=mm

         elseif (state.eq.2) then
            
            call mrrrr(m,tm,mm)
            m=mm
            
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
 enddo        
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chequea cadenas
! Selfavoiding entre segmentos
 
      dista=0.0
      do ive=4,long
         do jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
!cccccccccwrite(2,*)'noself-a.',k
               goto 222
            endif
        enddo
        enddo

      ncha=0
      do i=1,300
         test='S'
         call rota36(xend,xendr,long,test)
            ncha=ncha+1
            
            do j=1,long
               chains(1,j,ncha)=xendr(1,j) ! y
               chains(2,j,ncha)=xendr(2,j) ! x
               chains(3,j,ncha)=xendr(3,j) ! z
            enddo

         if (ncha.eq.25)exit
       enddo
      if (ncha.eq.0) goto 223

      return
end subroutine cadenas72mr_rosen
