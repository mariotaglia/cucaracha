subroutine print_ent(xend)
use globales
implicit none

real*8 xend(3,200)
integer i
character*21 filename
integer indexncha
indexncha = 1


! Imprime cadenas en formato ENT

write(filename,'(A6,A1, I3.3, A4)') 'cadena','.', indexncha, '.ent'

open(unit=4400, file=filename)

do i=1, long ! Imprime todo
WRITE(4400,'(A6,I5,A3,I12,A4,F8.3,F8.3,F8.3)') &
"HETATM",i,"  C",i,"    ",xend(1, i)*10,  &
xend(2, i)*10,xend(3, i)*10
end do

do i = 1, long ! Une segmentos
WRITE((4400),'(A6,I5,I5)')"CONECT", i, i+1
end do

WRITE(4400,*)"END"

close(4400)
indexncha = indexncha + 1
if(indexncha.eq.1000)stop
end subroutine

