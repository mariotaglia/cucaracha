!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculo de conductancia 
! 2015-08-12 FMG: Conductance output is in S units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_conductance(pHbulk,cpolbulk)
#   include "control_run.h"
!
! Calculo de conductancia !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
use globales
use pore
use csys
real(kind=8), intent(in) :: pHbulk
real(kind=8) :: Gvacio, Grel, Gporo, Gporopos, Gporoneg, GporoOHmin, GporoHplus
real(kind=8) :: Grel_coefD, Gporo_coefD, Gporo_coefDpos, Gporo_coefDneg, Gporo_coefDOHmin, Gporo_coefDHplus
integer :: iR
print*, "Calculating conductance"
Gvacio = 0.0
Grel = 0.0
Gporo = 0.0
Gporopos = 0.0
Gporoneg = 0.0
GporoHplus = 0.0
GporoOHmin = 0.0

Grel_coefD = 0.0
Gporo_coefD = 0.0
Gporo_coefDpos = 0.0
Gporo_coefDneg = 0.0
Gporo_coefDHplus = 0.0
Gporo_coefDOHmin = 0.0

Gvacio = Gvacio + xposbulk/vsalt/vsol*movpos
Gvacio = Gvacio + xnegbulk/vsalt/vsol*movneg
Gvacio = Gvacio + xHplusbulk/vsol*movHplus
Gvacio = Gvacio + xOHminbulk/vsol*movOHmin
Gvacio = Gvacio * 1e24/Na ! Corrige unidades concentracion
Gvacio = Gvacio * 1.0d-18/1d-6 ! Corrige unidades pasa nm2 -> m2 y micrometros -> m
Gvacio = Gvacio / longporo * pi * radio**2 ! Corrige unidades concentracion
!# Gvacio queda expresado en siemmens

do iR = 1, dimR
         Gporo = Gporo     +xpos(iR)/vsalt/vsol*movpos  *(dfloat(iR)-0.5)
      Gporopos = Gporopos  +xpos(iR)/vsalt/vsol*movpos  *(dfloat(iR)-0.5)
         Gporo = Gporo     +xneg(iR)/vsalt/vsol*movneg  *(dfloat(iR)-0.5)
      Gporoneg = Gporoneg  +xneg(iR)/vsalt/vsol*movneg  *(dfloat(iR)-0.5)
         Gporo = Gporo     + xHplus(iR) / vsol *movHplus*(dfloat(iR)-0.5)
    GporoHplus = GporoHplus+ xHplus(iR) / vsol *movHplus*(dfloat(iR)-0.5)
         Gporo = Gporo     + xOHmin(iR) / vsol *movOHmin*(dfloat(iR)-0.5)
    GporoOHmin = GporoOHmin+ xOHmin(iR) / vsol *movOHmin*(dfloat(iR)-0.5)
enddo
     Gporo = Gporo * 1e24/Na ! Corrige unidades concentracion
     Gporo = Gporo * 1.0d-18/1d-6 ! Corrige unidades
     Gporo = Gporo / longporo * 2*pi * delta**2
  
  Gporopos = Gporopos * 1e24/Na ! Corrige unidades concentracion
  Gporopos = Gporopos * 1.0d-18/1d-6 ! Corrige unidades
  Gporopos = Gporopos / longporo * 2*pi * delta**2
  
  Gporoneg = Gporoneg * 1e24/Na ! Corrige unidades concentracion
  Gporoneg = Gporoneg * 1.0d-18/1d-6 ! Corrige unidades
  Gporoneg = Gporoneg / longporo * 2*pi * delta**2

GporoHplus = GporoHplus * 1e24/Na ! Corrige unidades concentracion
GporoHplus = GporoHplus * 1.0d-18/1d-6 ! Corrige unidades
GporoHplus = GporoHplus / longporo * 2*pi * delta**2

GporoOHmin = GporoOHmin * 1e24/Na ! Corrige unidades concentracion
GporoOHmin = GporoOHmin * 1.0d-18/1d-6 ! Corrige unidades
GporoOHmin = GporoOHmin / longporo * 2*pi * delta**2

Grel = Gporo/Gvacio

! Conductance calculus with mobilities dependent on polymer volume fraction 
! ************************************************************************
do iR = 1, dimR
         Gporo_coefD = Gporo_coefD     +xpos(iR)/vsalt/vsol*( ((1-avpol(iR))/(1+avpol(iR)))**2)*movpos  *(dfloat(iR)-0.5)
      Gporo_coefDpos = Gporo_coefDpos  +xpos(iR)/vsalt/vsol*( ((1-avpol(iR))/(1+avpol(iR)))**2)*movpos  *(dfloat(iR)-0.5)
         Gporo_coefD = Gporo_coefD     +xneg(iR)/vsalt/vsol*( ((1-avpol(iR))/(1+avpol(iR)))**2)*movneg  *(dfloat(iR)-0.5)
      Gporo_coefDneg = Gporo_coefDneg  +xneg(iR)/vsalt/vsol*( ((1-avpol(iR))/(1+avpol(iR)))**2)*movneg  *(dfloat(iR)-0.5)
         Gporo_coefD = Gporo_coefD     + xHplus(iR) / vsol *( ((1-avpol(iR))/(1+avpol(iR)))**2)*movHplus*(dfloat(iR)-0.5)
    Gporo_coefDHplus = Gporo_coefDHplus+ xHplus(iR) / vsol *( ((1-avpol(iR))/(1+avpol(iR)))**2)*movHplus*(dfloat(iR)-0.5)
         Gporo_coefD = Gporo_coefD     + xOHmin(iR) / vsol *( ((1-avpol(iR))/(1+avpol(iR)))**2)*movOHmin*(dfloat(iR)-0.5)
    Gporo_coefDOHmin = Gporo_coefDOHmin+ xOHmin(iR) / vsol *( ((1-avpol(iR))/(1+avpol(iR)))**2)*movOHmin*(dfloat(iR)-0.5)
enddo

     Gporo_coefD = Gporo_coefD * 1e24/Na ! Corrige unidades concentracion
     Gporo_coefD = Gporo_coefD * 1.0d-18/1d-6 ! Corrige unidades
     Gporo_coefD = Gporo_coefD / longporo * 2*pi * delta**2
  
  Gporo_coefDpos = Gporo_coefDpos * 1e24/Na ! Corrige unidades concentracion
  Gporo_coefDpos = Gporo_coefDpos * 1.0d-18/1d-6 ! Corrige unidades
  Gporo_coefDpos = Gporo_coefDpos / longporo * 2*pi * delta**2
  
  Gporo_coefDneg = Gporo_coefDneg * 1e24/Na ! Corrige unidades concentracion
  Gporo_coefDneg = Gporo_coefDneg * 1.0d-18/1d-6 ! Corrige unidades
  Gporo_coefDneg = Gporo_coefDneg / longporo * 2*pi * delta**2

Gporo_coefDHplus = Gporo_coefDHplus * 1e24/Na ! Corrige unidades concentracion
Gporo_coefDHplus = Gporo_coefDHplus * 1.0d-18/1d-6 ! Corrige unidades
Gporo_coefDHplus = Gporo_coefDHplus / longporo * 2*pi * delta**2

Gporo_coefDOHmin = Gporo_coefDOHmin * 1e24/Na ! Corrige unidades concentracion
Gporo_coefDOHmin = Gporo_coefDOHmin * 1.0d-18/1d-6 ! Corrige unidades
Gporo_coefDOHmin = Gporo_coefDOHmin / longporo * 2*pi * delta**2

Grel_coefD = Gporo_coefD/Gvacio

write(315,*) cpolbulk, Gporo, Gporo_coefD
write(316,*) cpolbulk, Gvacio
write(317,*) cpolbulk, Grel, Grel_coefD

write(320,*) cpolbulk, Gporopos, Gporo_coefDpos
write(321,*) cpolbulk, Gporoneg, Gporo_coefDneg
write(322,*) cpolbulk, GporoHplus, Gporo_coefDplus
write(323,*) cpolbulk, GporoOHmin, Gporo_coefDOHmin
! output.aux
    write(324,*) cpolbulk, (Gporo - Gvacio), 7.352*(sigmaq*delta/vsol)*(3.141592*dimR/12000)/(6.02*1.0d7)
end subroutine calc_conductance
