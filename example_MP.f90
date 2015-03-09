!This is file : example_MP.f90
program main
  use mpmodule
  implicit none
  integer, parameter :: sp=kind(1.0), dp = selected_real_kind(2*sp), qp = selected_real_kind(2*dp)
  real(kind=dp) :: creal, dreal, creal2, dreal2
  type (mp_real) a, b, c, d
  call mpinit (500)
!  a = 1.d0
  c = '7213.423d307'
!  creal = 7213.423d307
  creal2= c
  call mpwrite(8,c)
!  b = cos(a)**2 + sin(a)**2 - 1.d0
  d = log(c)
!  dreal = dlog(creal)
  dreal2 = d
  call mpwrite(8,c,d)
!  call mpwrite(6, b)
!   write(8,*) creal, dreal
   write(8,*) creal2, dreal2
  stop
end program
