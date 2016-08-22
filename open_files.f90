subroutine open_files(m,ipH)
#   include "control_run.h"
    use csys, only: pHs
    integer, intent(in) :: m, ipH
    real :: pHbulk
    character(len=14) :: format_string
    character(len=9) :: pH_dat_str
    character(len=15) :: name_301
    character(len=16) :: name_302
    character(len=18) :: name_303
    character(len=18) :: name_304
    character(len=16) :: name_305
    character(len=17) :: name_306
    character(len=18) :: name_206
    character(len=16) :: name_307
    character(len=14) :: name_308
    character(len=18) :: name_319
    character(len=19) :: name_202
    character(len=25) :: name_203
    character(len=15) :: name_309
    character(len=16) :: name_201
    character(len=15) :: name_310
    character(len=19) :: name_311
    character(len=16) :: name_312
    character(len=16) :: name_313
    character(len=14) :: name_314
    character(len=16) :: name_318
    character(len=15) :: name_315
    character(len=16) :: name_316
    character(len=14) :: name_317
    character(len=14) :: name_320
    character(len=14) :: name_321
    character(len=16) :: name_322
    character(len=16) :: name_323
    character(len=20) :: name_324
    character(len=18) :: name_325

    
    select case ( m)
        case (1)
            pHbulk = pHs(ipH)
            format_string='(I2.2,F0.2,A4)'
            write(pH_dat_str,format_string) int(pHbulk) , pHbulk-int(pHbulk), '.dat'
        ! Energy 
            write(name_301,'(A6, A9)') 'F_tot_', pH_dat_str
            open(unit=301, file=name_301)
           !open(unit=301, file='F_tot.dat')
            write(name_302,'(A7, A9)') 'F_mixs_', pH_dat_str
            open(unit=302, file=name_302)
            !open(unit=302, file='F_mixs.dat')

            write(name_303,'(A9, A9)') 'F_mixpos_', pH_dat_str
            open(unit=303, file=name_303)
            !open(unit=303, file='F_mixpos.dat')
            write(name_304,'(A9, A9)') 'F_mixneg_', pH_dat_str
            open(unit=304, file=name_304)
            !open(unit=304, file='F_mixneg.dat')
            write(name_305,'(A7, A9)') 'F_mixH_', pH_dat_str
            open(unit=305, file=name_305)
            !open(unit=305, file='F_mixH.dat')
            write(name_306,'(A8, A9)') 'F_mixOH_', pH_dat_str
            open(unit=306, file=name_306)
            !open(unit=306, file='F_mixOH.dat')
            write(name_206,'(A9, A9)') 'F_mixpol_', pH_dat_str
            open(unit=206, file=name_206)
            !open(unit=206, file='F_mixOH.dat')
            write(name_307,'(A7, A9)') 'F_conf_', pH_dat_str
            open(unit=307, file=name_307)
            !open(unit=307, file='F_conf.dat')
            write(name_308,'(A5, A9)') 'F_eq_', pH_dat_str
            open(unit=308, file=name_308)
            !open(unit=308, file='F_eq.dat')
            write(name_319,'(A9, A9)') 'F_eqwall_', pH_dat_str
            open(unit=319, file=name_319)
            !open(unit=319, file='F_eqwall.dat')
            write(name_202,'(A10, A9)') 'sys_mupol_', pH_dat_str
            open(unit=202, file=name_202)
            !open(unit=202, file='sys_mupol.dat')

!            write(name_203,'(A16, A9)') 'adsorved_chains_', pH_dat_str
!            open(unit=203, file=name_203)
            !open(unit=203, file='adsorved_chains.dat')
# ifdef VDW
            write(name_309,'(A6, A9)') 'F_vdW_', pH_dat_str
            open(unit=309, file=name_309)
!            open(unit=309, file='F_vdW.dat')
# endif
            write(name_201,'(A7, A9)') 'F_ospi_', pH_dat_str
            open(unit=201, file=name_201)
            !open(unit=201, file='F_ospi.dat')
            write(name_310,'(A6, A9)') 'F_eps_', pH_dat_str
            open(unit=310, file=name_310)
            !open(unit=310, file='F_eps.dat')
            write(name_311,'(A10, A9)') 'F_electro_', pH_dat_str
            open(unit=311, file=name_311)
            !open(unit=311, file='F_electro.dat')
            write(name_312,'(A7, A9)') 'F_tot2_', pH_dat_str
            open(unit=312, file=name_312)
            !open(unit=312, file='F_tot2.dat')
        ! Mean values
            write(name_313,'(A7, A9)') 'fmedio_', pH_dat_str
            open(unit=313, file=name_313)
            !open(unit=313, file='fmedio.dat')
#if CHAIN != 0
            write(313,*) '# pHbulk, fmedio, fdiswall, qwall'
#else            
            write(313,*) '# pHbulk, fdiswall, qwall'
#endif
!            write(name_314,'(A5, A9)') 'pKas_', pH_dat_str
!            open(unit=314, file=name_314)
            !open(unit=314, file='pKas')
            write(name_318,'(A7, A9)') 'Rmedio_', pH_dat_str
            open(unit=318, file=name_318)
            !open(unit=318, file='Rmedio')
            write(318,*) '# pHbulk, Rmedio, sumcharge'
        ! Conductance
            write(name_315,'(A6, A9)') 'Gporo_', pH_dat_str
            open(unit=315, file=name_315)
!            open(unit=315, file='Gporo')
            write(name_316,'(A7, A9)') 'Gvacio_', pH_dat_str
            open(unit=316, file=name_316)
!            open(unit=316, file='Gvacio')
            write(name_317,'(A5, A9)') 'Grel_', pH_dat_str
            open(unit=317, file=name_317)
!            open(unit=317, file='Grel')
            write(name_320,'(A5, A9)') 'Gpos_', pH_dat_str
            open(unit=320, file=name_320)
!            open(unit=320, file='Gpos')
            write(name_321,'(A5, A9)') 'Gneg_', pH_dat_str
            open(unit=321, file=name_321)
!            open(unit=321, file='Gneg')
            write(name_322,'(A7, A9)') 'GHplus_', pH_dat_str
            open(unit=322, file=name_322)
!            open(unit=322, file='GHplus')
            write(name_323,'(A7, A9)') 'GOHmin_', pH_dat_str
            open(unit=323, file=name_323)
!            open(unit=323, file='GOHmin')
! AUXILIAR FILES
 !           write(name_324,'(A11, A9)') 'output_aux_', pH_dat_str
 !           open(unit=324, file=name_324)
!            open(unit=324, file='output.aux')
            write(name_325,'(A9, A9)') 'sigmapol_', pH_dat_str
            open(unit=325, file=name_325)
            write(325,*) "#cpol, sigmapol # Mupol_cte" 
        case ( 0 )
        ! Energy
            close(201) ! F_ospi(?) 
            close(202) ! sys_mupol
           ! close(203) ! adsorved_chains
            close(301) ! F_tot 
            close(302) ! F_mixs
            close(303) ! F_mixpos
            close(304) ! F_mixneg
            close(305) ! F_mixH
            close(306) ! F_mixOH
            close(206) ! F_mixpol
            close(307) ! F_conf
            close(308) ! F_eq
#       ifdef VDW
            close(309) ! F_vdw
#       endif
            close(310) ! F_eps 
            close(311) ! F_electro
            close(312) ! F_tot2 

        ! Conductance and mean values (313 - 323)
            close(313) ! fmedio 
            !close(314) ! pKas
            close(315) ! Gporo
            close(316) ! Gvacio 
            close(317) ! Grel
            close(318) ! Rmedio
            close(319) ! F_eqwall
            close(320) ! Gpos
            close(321) ! Gneg
            close(322) ! GHplus
            close(323) ! GOHmin
! AUXILIAR FILES
  !          close(324) ! output.aux
            close(325) ! output.aux
    end select 
     
end subroutine
