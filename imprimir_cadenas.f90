subroutine imprimir_cadenas()
    use globales, only: cuantas, long
    use csys, only: in1
    
    integer  :: il,j
    il = 1  
    print*, 'Se imprimen las conformaciones en conformations.xyx, use: vmd file.xyz '

    open(unit=678, file='conformations.xyz')
! Los comentarios en las primeras lineas no estan soportados en VMD
!    write(678,'(A73)') '# intenta formato xyz como en http://en.wikipedia.org/wiki/XYZ_file_format'
!    write(678,'(A32,I7)') '# Total number of Conformations: ', cuantas

    do while (il.le.cuantas)
        write(678,'(I2)') long
        write(678,'(A5,I7,A16)') '#Conf', il, ', position in nm'

        do j=1,long         
            write(678,'(A1,3F7.3)') 'C', in1(il,j,:)
        end do
        il=il+1
!        write(678, *) " "
    end do
    close(678)
end subroutine
