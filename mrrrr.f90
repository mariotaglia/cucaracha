!      
!C****************************************************************
!
!

subroutine mrrrr(a,b,c)
    implicit none
      
    REAL(KIND=8), dimension(3,3), intent(in) :: a,b
    REAL(KIND=8), dimension(3,3), intent(out) :: c
    integer :: i, j, k

    do 1 i=1,3
        do 1 j=1,3
            c(i,j)=0
 1    continue

      do 2 i=1,3
         do 2 j=1,3
            do 2 k=1,3
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
 2    continue

      return
end subroutine mrrrr
