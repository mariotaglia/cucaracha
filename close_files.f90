subroutine close_files

      close(unit=313)! , file='fmedio.dat')
      close(unit=314)! , file='pKas')
      close(unit=315)! , file='Gporo')
      close(unit=316)! , file='Gvacio')
      close(unit=317)! , file='Grel')
      close(unit=318)! , file='Rmedio')
      close(unit=320)! , file='Gpos')
      close(unit=321)! , file='Gneg')
      close(unit=322)! , file='GHplus')
      close(unit=323)! , file='GOHmin')

end subroutine
