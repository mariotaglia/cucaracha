! Here all variables related with bulk

subroutine set_bulk_properties(pHbulk)
#   include "control_run.h"
    use globales
    use csys
    implicit none
    real(KIND=8), intent(in) :: pHbulk

        cHplus = 10**(-pHbulk) ! concentration H+ in bulk
        xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
! Relacion que sale del Atkins íntrinseco a las definiciones de equilibrio con pH
        pOHbulk= pKw -pHbulk 
        cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
!------- Volume fractions
        cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
        xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
        if(pHbulk.le.7) then  ! pH<= 7
            xposbulk= xsalt/zpos
            !xnegbulk=-xsalt/zneg -xsalt2/zneg +(xHplusbulk -xOHminbulk)*vsalt ! NaCl+ HCl
            xnegbulk=-xsalt/zneg +(xHplusbulk -xOHminbulk)*vsalt ! NaCl+ HCl
        else                  ! pH >7
            xposbulk= xsalt/zpos +(xOHminbulk -xHplusbulk)*vsalt ! NaCl+ NaOH
            !xnegbulk=-xsalt/zneg -xsalt2/zneg
            xnegbulk=-xsalt/zneg
        endif
! Totalmente add-hoc y util solo para calcular mupol 11-6-2015
        xpolbulk = (cpol*Na/(1.0d24))*vpol ! mupol 11-6-2015
!!#if MUPOL == 0
! Polymer Chemical potential in bulk
           xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk -xpolbulk ! mupol 11-6-2015
          expmupol = ( cpol*0.602*vsol/( xsolbulk**(vpol*long) ) )            ! mupol 11-6-2015

! En bulk no hay polimero.
         xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk !- xposbulk2
!! NO hay polimero en solucion!
!! #elif MUPOL == 1 /* Mal! una vez que se funcionaliza el nanocanal ya no hay poliemro en la solucion de medicion */
!!          xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk -xpolbulk !- xposbulk2
! Polymer chemical potential!
!!           expmupol = xpolbulk/vpol / xsolbulk**(vpol*long)
!!        print*, "expmupol ", pHbulk, expmupol
!!#endif

! Ojo Kw y Ka estan en mol/l mientras que vsol y xsolbulk estan en nm3
! 1.0d24 factor de conversion de volumen de litro a nm3  1nm3 = 10^-24 l = 10^-24 dm3

! remember that vsalt, vHplus and vOHmin are in vsol unit!
           expmupos = xposbulk/vsalt / xsolbulk**vsalt 
           expmuneg = xnegbulk/vsalt / xsolbulk**vsalt 
!ojo! expmuHplus ya tiene el signo menos!!
         expmuHplus = xHplusbulk / xsolbulk**vHplus   ! vsol = vHplus
         expmuOHmin = xOHminbulk / xsolbulk**vOHmin  ! vsol = vOHmin
! intrinstic equilibruim constants in mol/nm3
            Ka0 = (Ka*vsol/xsolbulk)*(Na/1.0d24)
            Kb0 = (Kb*vsol/xsolbulk)*(Na/1.0d24)
         Kwall0 = (Kwall*vsol/xsolbulk)*(Na/1.0d24)

end subroutine set_bulk_properties
