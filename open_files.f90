subroutine open_files(m)
#   include "control_run.h"
!    use csys, only: no_chains
    integer, intent(in) :: m
    select case ( m)
        case (1)
        ! Energy 
            open(unit=301, file='F_tot.dat')
            open(unit=302, file='F_mixs.dat')
            open(unit=303, file='F_mixpos.dat')
            open(unit=304, file='F_mixneg.dat')
            open(unit=305, file='F_mixH.dat')
            open(unit=306, file='F_mixOH.dat')
            open(unit=307, file='F_conf.dat')
            open(unit=308, file='F_eq.dat')
# ifdef VDW
            open(unit=309, file='F_vdW.dat')
# endif
            open(unit=310, file='F_eps.dat')
            open(unit=311, file='F_electro.dat')
            open(unit=312, file='F_tot2.dat')
        ! Mean values
            open(unit=313, file='fmedio.dat')
            write(313,*) '# pHbulk, fmedio, fmedio2, fdiswall'
            open(unit=314, file='pKas')
            open(unit=318, file='Rmedio')
            write(318,*) '# pHbulk, Rmedio, sumcharge'
        ! Conductance
            open(unit=315, file='Gporo')
            open(unit=316, file='Gvacio')
            open(unit=317, file='Grel')
            open(unit=320, file='Gpos')
            open(unit=321, file='Gneg')
            open(unit=322, file='GHplus')
            open(unit=323, file='GOHmin')
        case ( 0 )
        ! Energy 
            close(301)
            close(302)
            close(303)
            close(304)
            close(305)
            close(306)
            close(307)
            close(308)
#       ifdef VDW
            close(309)
#       endif
            close(310)
            close(311) 
            close(312)

        ! Conductance and mean values
            close(313) 
            close(314)
            close(315)
            close(316)
            close(317)
            close(318)
            close(320)
            close(321)
            close(322)
            close(323)
    end select 
     
end subroutine
