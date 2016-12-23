function fmix()
    use globales, only: dimR, Radio, delta, vsol, vsalt
    use csys, only: xsolbulk, expmuHplus, xHplusbulk, expmuOHmin, xOHminbulk &
                            , expmupos, xposbulk, expmuneg, xnegbulk
    use pore, only: xh, xHplus, xOHmin, xpos, xneg
    use FreeEnergy, only: F_Mix_s, F_Mix_pos,  F_Mix_neg,  F_Mix_Hplus,  F_Mix_OHmin
 
    implicit none
    real(kind=8) :: fmix, fmixs, fmixHplus, fmixOHmin, fmixpos, fmixneg
    real(kind=8) :: rdrR
    integer :: iR
! El output de la funcion es la suma de todo esto
    fmix=0
    fmixs=0
    fmixHplus=0
    fmixOHmin=0
    fmixpos=0
    fmixneg=0
    

! Siempre se calcula la energia respecto de la de bulk
    do iR = 1, dimR
        rdrR = delta*(dfloat(iR)-0.5+radio/delta)*delta/Radio
       
        fmixs = fmixs + rdrR* ( xh(iR) /vsol) *( dlog(xh(iR)) -1.0) 
        fmixs = fmixs - rdrR* (xsolbulk/vsol) *(dlog(xsolbulk)-1.0)
        ! H+
        fmixHplus = fmixHplus + (xHplus(iR)/vsol)*(dlog(xHplus(iR))-1.0 -dlog(expmuHplus)) *rdrR
        fmixHplus = fmixHplus - (xHplusbulk/vsol)*(dlog(xHplusbulk)-1.0 -dlog(expmuHplus)) *rdrR
        ! OH-
        fmixOHmin = fmixOHmin + (xOHmin(iR)/vsol)*(dlog(xOHmin(iR))-1.0 -dlog(expmuOHmin)) *rdrR
        fmixOHmin = fmixOHmin - (xOHminbulk/vsol)*(dlog(xOHminbulk)-1.0 -dlog(expmuOHmin)) *rdrR
! Ver NOTAR sobre mix entropy recordar concentracion referencia 1/vsol 
        ! K+
        fmixpos = fmixpos + (xpos(iR)/(vsalt*vsol)) *(dlog(xpos(iR)/vsalt) -1.0 -dlog(expmupos) ) *rdrR
        fmixpos = fmixpos - (xposbulk/(vsalt*vsol)) *(dlog(xposbulk/vsalt) -1.0 -dlog(expmupos) ) *rdrR
        ! Cl-
        fmixneg = fmixneg + (xneg(iR)/(vsalt*vsol)) *(dlog(xneg(iR)/vsalt) -1.0 -dlog(expmuneg) ) *rdrR
        fmixneg = fmixneg - (xnegbulk/(vsalt*vsol)) *(dlog(xnegbulk/vsalt) -1.0 -dlog(expmuneg) ) *rdrR
    enddo
        fmix = fmixs + fmixHplus + fmixOHmin + fmixpos + fmixneg
!    print*, "fmixs llama checknumber: fmixs", fmixs
!    call checknumber(fmixs,'fmixs')
!    call checknumber(fmixHplus, 'Energia fmixHplus')
!    call checknumber(fmixOHmin, 'Energia fmixOHmin')
!    call checknumber(fmixpos, 'fmixpos')
!    call checknumber(fmixneg, 'fmixneg')
    F_Mix_s     = fmixs 
    F_Mix_pos   = fmixHplus
    F_Mix_neg   = fmixOHmin
    F_Mix_Hplus = fmixpos
    F_Mix_OHmin = fmixneg
    return
    contains 
        subroutine checknumber(var, arg)
            implicit none
            real(kind=8), intent(in) :: var
            character(len=*), intent(in) :: arg
            ! Check if Not a Number
            if ( var /= var ) then
                print*, arg, " real number is NaN"
                call printstate('checknumber real NaN')
                stop
            endif
        
            if ( var-1 == var ) then
                print*, arg, " real number is infinity"
                call printstate('checknumber real infinity')
                stop
            endif
        endsubroutine checknumber
end function fmix
