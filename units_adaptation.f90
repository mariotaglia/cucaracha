! Here the main idea is to adapt the units of the input parameters

subroutine units_adaptation
    use globales
    use csys
    implicit none

    sigma = sigma *vsol/ delta
    sigmaq = sigmaq ! chequear esto?!
    constq=delta*delta*4.0*pi*lb/vsol ! multiplicative factor in poisson eq units? without units

! Chemical Equilibrium
    Ka=10**(-pKa)
    Kb=10**(-pKb)
    Kwall=10**(-pKawall)

! Volume fraction second salt in mol/l, csalt2=0! Check!
      xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   
!      xsalt2=(csalt2*Na/(1.0d24))*(vsalt2*vsol)   

    print*, 'Units adaptation OK!'
end subroutine
