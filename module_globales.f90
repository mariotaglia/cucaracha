! El objetivo es poner aca todas las variables globales del sistema
! constante y propiedades fisicas

module globales
    implicit none
! System Parameters, constant, etc
    real(KIND=8), parameter :: &
        Na=6.02d23, & ! Avogadro Number
        pi=3.14159265359, &
!        infinity = HUGE(dbl_prec_var), & ! computer infinity
! Charges        
        zpos=1.0, zpos2=3.0, zneg=-1.0, &
        zwall=-1.0, zpol=1.0, &
        zH=1.0, zOH=-1.0, &
! Lengths
        lb=0.714, & !Bjerrum Length (nm)
        delta=0.5, & ! Layer width (nm)
        lseg=0.5, & ! Largo del segmento (nm)

! Volumes        
        vsol=1.0, &
        vpol=1.0, &
!        vsol=0.030, &
!        vpol= 0.095/vsol, &  ! Polymer volume segment in units of vsol: ((4.0/3.0)*pi*(0.283)**3)/vsol
!        vHplus=0.030/vsol, vOHmin=0.030/vsol, & ! vHplus=1 ! vOHplus=1 
        vHplus=1.0, vOHmin=1.0, &
        vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol, &  ! Salt volume in units of vsol 0.2nm=radius salt
!        vsalt2=((4.0/3.0)*pi*(0.27)**3)/vsol, & ! vpol > vsalt > vsol
        vsalt2 = 0.0, &
       
        pKw = 14, Kw = 1.0e-14, &
       
        error = 1e-6, &  ! para comparar con la norma...
        betae = 38.94, & ! beta * e
        errel=1d-6
!                                                        vsol,   vH+,  vOH-,    vK+,   vCl-,  vPol
!    real(KIND=8), dimension(6), parameter :: v_all = (/ 0.030, 0.030, 0.030, 0.0824, 0.0824, 0.095 /)
        
!    real(KIND=8), dimension(5), parameter :: v_all = (/ 0.030, 0.030, 0.030, ( (4.0/3.0)*pi*(0.27**3)) , 0.095 /)
!    real(KIND=8) :: radio=dimR*delta
    real(KIND=8) :: radio

    integer, parameter :: itmax=200
! System parameter (INPUTS)
    integer :: long, cuantas ! numero de layers del radio del poro  (tiene que ser par?)
!    integer, parameter :: long=14, dimR=20, cuantas=100 ! numero de layers del radio del poro  (tiene que ser par?)
! Kinsol Parameters
!    integer :: neq = 2*dimR
    integer :: dimR, neq

end module globales
