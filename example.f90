!This is file : example.f90
program example
    implicit None
    integer, parameter :: sp = selected_real_kind(1),     &
                          dp = selected_real_kind(2*sp),    &
                          qp = selected_real_kind(2*dp),    &
                          Bp = selected_real_kind(15,308)  
    
    real(kind=8) ::  ereal
    real(kind=sp) :: sreal
    real(kind=dp) :: dreal
    real(kind=qp) :: qreal
    real(kind=Bp) :: breal
    
    print*, "Precision de un kind=sp: ", precision(sreal) , "Rango de un kind=sp: ",     range(sreal)
    print*, "Precision de un kind=8 : ", precision(ereal) , "Rango de un kind=8 : ",     range(ereal)
    print*, "Precision de un kind=dp: ", precision(dreal) , "Rango de un kind=dp: ",     range(dreal)
    print*, "Precision de un kind=qp: ", precision(qreal) , "Rango de un kind=qp: ",     range(qreal)
    print*, "Precision de un kind=bp: ", precision(breal) , "Rango de un kind=bp: ",     range(breal)

    contains
        subroutine imprimir(a)
            implicit none
            real, intent(in) :: a
            print*, "Precision de un kind=8: ", precision(a)
            print*, "Rango de un kind=8: ", range(a)
        end subroutine imprimir    
end program example
