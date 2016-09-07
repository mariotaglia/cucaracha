subroutine imprimir_cadenas(in1)
    use globales, only: cuantas, long
    real*8 in1(long,3) 
    integer  :: j
    il = 1  
    print*, 'Se imprimen las conformaciones en conformations.xyz, use: vmd conformations.xyz '

    open(unit=678, file='conformations.xyz')
! Los comentarios en las primeras lineas no estan soportados en VMD
!    write(678,'(A73)') '# intenta formato xyz como en http://en.wikipedia.org/wiki/XYZ_file_format'
!    write(678,'(A32,I7)') '# Total number of Conformations: ', cuantas

        write(678,'(I2)') long
        write(678,'(A5,I7,A16)') '#Conf position in nm'

        do j=1,long         
            write(678,'(A1,3F7.3)') 'C', in1(j,:)
        end do
!        write(678, *) " "
    close(678)
end subroutine
