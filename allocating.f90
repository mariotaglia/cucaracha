!
! Easy way to remember to  clean the memory
!
!
subroutine allocating(m)
    use globales, only: dimR, long, cuantas
    use csys
    use pore
    implicit none
    integer, intent(in) :: m
    integer :: var1=0

    select case ( m ) 
        case ( 0 ) ! Deallocating 
            ! Declared in module csys
            deallocate( pHs, pR, in1, Xu,pro,stat=var1 )
            if ( var1 /= 0 ) print*, "There was an error in  memory deallocation. See last lines in main program. "
            var1=0
            deallocate( x1, stat=var1)
            !deallocate( x1, xg1, xflag, stat=var1)
            if ( var1 /= 0 ) print*, "There was an error in  memory deallocation. See last lines in main program. "
            deallocate(xh,    &  
                     avpol, &        
                     xpot,  &
                     qtot,  &
                     psi, &
                     xpos,  &
                     xpos2, &
                     xneg,  &
                     xHplus,&
                     xOHmin,&
                     xtotal,&
                     fdis, fdis2, pp, stat=var1)
            if ( var1 /= 0 ) print*, "allocating.f90: There is an erro while deallocating memory for variables: xh - fdis" 
        case ( 1 ) ! Allocating
            allocate( x1(2*dimR), stat=var1)
            !allocate( x1(2*dimR), xg1(2*dimR), xflag(2*dimR) , stat=var1)
            if ( var1 /= 0 ) print*, "There is no sufficient memory for variables: x1, xg1, xflag" 
            var1=0
            allocate( pR(cuantas,long), in1(cuantas,long,3), Xu(dimR,dimR),pro(cuantas) )
            if ( var1 /= 0 ) print*, "There is no sufficient memory for variables: x1, xg1, xflag" 
            var1=0
            allocate(xh(dimR),    &  
                     avpol(dimR), &        
                     xpot(dimR),  &
                     qtot(dimR),  &
                     psi(0:dimR+1), &
                     xpos(dimR),  &
                     xpos2(dimR), &
                     xneg(dimR),  &
                     xHplus(dimR),&
                     xOHmin(dimR),&
                     xtotal(dimR),&
                     fdis(dimR), fdis2(dimR), pp(2*dimR), stat=var1)
            if ( var1 /= 0 ) print*, "allocating.f90: There is no sufficient memory for variables: xh - fdis" 

    end select
! End allocating

end subroutine allocating
