!This is file : checknum.f90
subroutine checknum(var, arg)
    implicit none
    real(kind=8), intent(in) :: var
    character(len=*), intent(in) :: arg
    ! Check if Not a Number
    
    if ( var /= var ) then
        print*, arg, " real number is NaN"
        !call printstate('checknumber real NaN')
        call printstate(arg)
        stop
    endif

    if ( var-1 == var ) then
        print*, arg, " real number is infinity"
        !call printstate('checknumber real infinity')
        call printstate(arg)
        stop
    endif
endsubroutine checknum
