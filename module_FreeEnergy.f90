
module FreeEnergy
    real(kind=8) :: F_Mix, F_Mix_s, F_Mix_pos, F_Mix_neg, F_Mix_Hplus, F_Mix_OHmin, F_Mix_pol,&
                    F_Conf, F_Eq, F_Eq_wall, F_vdW, F_electro, F_ospi, F_eps
    real(kind=8) :: Free_Energy, Free_Energy2
! En contains va el codigo de cada funcion!
    interface
        real(kind=8) function fmix()
        end function
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
        real(kind=8) function fospi()
        end function
        real(kind=8) function pong_energy()
        end function
    end interface

end module FreeEnergy
