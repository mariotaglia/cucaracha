subroutine imprimir_cadenas()
    use globales, only: cuantas, long
    use csys, only: in1
    
    integer  :: il,j
    il = 1  
    print*, 'Se imprimen las conformaciones en conformations.xyx, use: vmd file.xyz '

    open(unit=314, file='conformations.xyz')
! Los comentarios en las primeras lineas no estan soportados en VMD
!    write(314,'(A73)') '# intenta formato xyz como en http://en.wikipedia.org/wiki/XYZ_file_format'
!    write(314,'(A32,I7)') '# Total number of Conformations: ', cuantas

    do while (il.lt.cuantas)
        write(314,'(I2)') long
        write(314,'(A5,I7,A16)') '#Conf', il, ', position in nm'

        do j=1,long         
            write(314,'(A1,3F7.3)') 'C', in1(il,j,:)
        end do
        il=il+1
!        write(314, *) " "
    end do
    close(314)
end subroutine
