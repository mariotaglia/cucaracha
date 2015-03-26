
subroutine set_initial_guess(m)
    use globales, only: dimR
    use csys, only: x1, xsolbulk, sigmaq
    use pore, only: xh, psi
    implicit none
    integer, intent(in) :: m
    integer :: i

    select case (m)
        case ( 0 )  ! Bulk solution
             xh(:dimR) = xsolbulk
!            xg1(:dimR) = xsolbulk
             x1(:dimR) = xsolbulk
   ! Psi Initial values
!          xg1(dimR+1:) = 0.0d0
             x1(dimR+1:) = 0.0d0
             psi(0:dimR+1) = 0.0d0
!            do i = 1, dimR
!                 xh(i) = xsolbulk
!                xg1(i) = xsolbulk
!                 x1(i) = xsolbulk
!            end do
!            if (sigmaq /= 0.0 ) then
!                do i = 1, dimR
!                    xg1(dimR+i) = -0.001
!                     x1(dimR+i) = -0.001
!                end do
!            else 
!                do i = 1, dimR
!                    xg1(dimR+i) = 0.0d0
!                     x1(dimR+i) = 0.0d0
!                end do
!            end if 
        case ( 1 ) 
        open(unit=92,file='initial_guess')
        read(92,*);read(92,*) ! ignore first and second lines
       
            do i=1,2*dimR
                read(92,*) x1(i)
                !read(92,*) xg1(i)
            end do
       
        close(92)
   ! xh and Psi Initial values are in xg1
         xh(:dimR) = x1(:dimR) 
         psi(0) = 0.
         psi(1:dimR)= x1(dimR+1:)
         psi(dimR+1) = 0.
         !xh(:dimR) = xg1(:dimR) 
         !psi(1:dimR)= xg1(dimR+1:)
         !x1(:) = xg1(:)
    end select
! Por unica vez. Luego es llamada desde fkfun
    call set_pore_distrib
end subroutine set_initial_guess
