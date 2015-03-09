
module FreeEnergy
    real(kind=8) :: F_Mix_pos, F_Mix_neg, F_Mix_Hplus, F_Mix_OHmin, &
                    F_Conf, F_Eq, F_Eq_wall, F_vdW, F_electro, F_eps
! En contains va el codigo de cada funcion!
    interface
        real(kind=8) function fmixs()
        end function
        real(kind=8) function fmixpos()
        end function
        real(kind=8) function fmixneg()
        end function
        real(kind=8) function fmixHplus()
        end function
        real(kind=8) function fmixOHmin()
        end function
        real(kind=8) function fconf_pol()
        end function
        real(kind=8) function fchem_eq()
        end function
        real(kind=8) function fvdW()
        end function
        real(kind=8) function fpol_sup()
        end function
        real(kind=8) function fchem_eq_wall()
        end function
        real(kind=8) function pong_energy()
        end function
    end interface

end module FreeEnergy
