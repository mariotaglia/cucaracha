!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculo de conductancia ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2015-08-12 FMG: Conductance output is in *** units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_conductance(pHbulk)
#   include "control_run.h"
!
! Calculo de conductancia !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    use globales
    use pore
    use csys
    real(kind=8), intent(in) :: pHbulk
    real(kind=8) :: Gvacio, Grel, Gporo,Gporopos, Gporoneg, GporoOHmin, gporoHplus
    integer :: iR
    print*, "Calculating conductance"
    Gvacio = 0.0
    Grel = 0.0
    Gporo = 0.0
    Gporopos = 0.0
    Gporoneg = 0.0
    GporoHplus = 0.0
    GporoOHmin = 0.0

    Gvacio = Gvacio + xposbulk/vsalt/vsol*movpos
    Gvacio = Gvacio + xnegbulk/vsalt/vsol*movneg
    Gvacio = Gvacio + xHplusbulk/vsol*movHplus
    Gvacio = Gvacio + xOHminbulk/vsol*movOHmin
    Gvacio = Gvacio * 1e24/Na ! Corrige unidades concentracion
    Gvacio = Gvacio * 1.0d-18/1d-6 ! Corrige unidades pasa nm2 -> m2 y micrometros -> m
    Gvacio = Gvacio / longporo * pi * radio**2 ! Corrige unidades concentracion
    Gvacio = Gvacio * 10 ! Deja el output en S/cm

    do iR = 1, dimR
        Gporo=Gporo+xpos(iR)/vsalt/vsol*movpos*(dfloat(iR)-0.5)
        Gporopos=Gporopos+xpos(iR)/vsalt/vsol*movpos*(dfloat(iR)-0.5)
        Gporo=Gporo+xneg(iR)/vsalt/vsol*movneg*(dfloat(iR)-0.5)
        Gporoneg=Gporoneg+xneg(iR)/vsalt/vsol*movneg*(dfloat(iR)-0.5)
        Gporo=Gporo+xHplus(iR)/vsol*movHplus*(dfloat(iR)-0.5)
        GporoHplus=GporoHplus+xHplus(iR)/vsol*movHplus*(dfloat(iR)-0.5)
        Gporo=Gporo+xOHmin(iR)/vsol*movOHmin*(dfloat(iR)-0.5)
        GporoOHmin=GporoOHmin+xOHmin(iR)/vsol*movOHmin*(dfloat(iR)-0.5)
    enddo

    Gporo = Gporo * 1e24/Na ! Corrige unidades concentracion
    Gporo = Gporo * 1.0d-18/1d-6 ! Corrige unidades
    Gporo = Gporo / longporo * 2*pi * delta**2
    Gporo = Gporo * 10 ! Deja el output en S/cm

    Gporopos = Gporopos * 1e24/Na ! Corrige unidades concentracion
    Gporopos = Gporopos * 1.0d-18/1d-6 ! Corrige unidades
    Gporopos = Gporopos / longporo * 2*pi * delta**2
    Gporopos = Gporopos *10 ! Deja el output en S/cm

    Gporoneg = Gporoneg * 1e24/Na ! Corrige unidades concentracion
    Gporoneg = Gporoneg * 1.0d-18/1d-6 ! Corrige unidades
    Gporoneg = Gporoneg / longporo * 2*pi * delta**2
    Gporoneg = Gporoneg * 10 ! Deja el output en S/cm

    GporoHplus = GporoHplus * 1e24/Na ! Corrige unidades concentracion
    GporoHplus = GporoHplus * 1.0d-18/1d-6 ! Corrige unidades
    GporoHplus = GporoHplus / longporo * 2*pi * delta**2
    GporoHplus = GporoHplus * 10 ! Deja el output en S/cm

    GporoOHmin = GporoOHmin * 1e24/Na ! Corrige unidades concentracion
    GporoOHmin = GporoOHmin * 1.0d-18/1d-6 ! Corrige unidades
    GporoOHmin = GporoOHmin / longporo * 2*pi * delta**2
    GporoOHmin = GporoOHmin * 10 ! Deja el output en S/cm 

    Grel = Gporo/Gvacio

    write(315,*) pHbulk, Gporo
    write(316,*) pHbulk, Gvacio
    write(317,*) pHbulk, Grel
 
    write(320,*) pHbulk, Gporopos
    write(321,*) pHbulk, Gporoneg
    write(322,*) pHbulk, GporoHplus
    write(323,*) pHbulk, GporoOHmin
    write(324,*) pHbulk, (Gporo - Gvacio), 7.352*(sigmaq*delta/vsol)*(3.141592*dimR/12000)/(6.02*1.0d7)
end subroutine calc_conductance
