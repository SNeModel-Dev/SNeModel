program lightcurve
    !
    implicit none
    !******* Comments through out the code

    double precision m(2000), energy(2000), v(2000), r(0:2000), &
            rho(2000), temp(2000), edot(2000), kappa(2000), vshr(2000)
    double precision rhosh0(2000), vsh0(2000), tsh0(2000)

    double precision rhosh(2000), vsh(2000), tsh(2000)
    double precision rhoshw, rhow, brad, dr(2000)
    double precision vshw, tshw, accel, msh, mw
    double precision a, pi43, pi4, sigma, &
            difffac, efac, tauni, rho0, alpha
    double precision rstar, v0, mexp, eexp, rold
    double precision mni, mcore, etest, gam, vwind, mdot
    character filename*128 
    integer jedge

    double precision mtot, esh(2000), rw, alpha2
    !
    !-- evolution variables
    !
    double precision dt, time, rhowi
    double precision tau, lum, diff, rphot, tphot, lumt, &
            den, lums, p, tsh1, rho1, vel
    double precision trecom, precom, factor(2000)
    !
    integer i, j, jcore, jphot
    !
    read(*,*) alpha,rstar,mexp,eexp,filename
    ! rstar = 8.d12
    !mexp = 6.d34
    !eexp = 5.d51
    trecom = 1.d0 * 1.16d4
    precom = 0.25d0
    mni = 2.d32

    mcore = mexp/10.0d0
    do i = 1, 2000
        kappa(i) = 0.4d0
    end do
    pi4 = 3.14159258d0 * 4.d0
    pi43 = 3.14159258d0 * 4.d0 / 3.d0
    a = 7.5656d-15
    difffac = a * 4.d0 * 3.14159258d0 * 3.d10&
            / 3.d0
    sigma = 5.6704d-5 * 4.0d0 * 3.14159258d0
    efac = 4.78d10
    tauni = 7.605d5
    gam = 4.d0 / 3.d0
    vwind = 4.d8
    mdot = 2.26d26
    !alpha = 2.0d0
    
    open(69, file=trim(filename)//'BB.dat',status='new')
    open(79, file=trim(filename)//'prop.dat',status='new')
    open(89, file=trim(filename)//'initial.dat',status='new')
    102  format(5(1pe12.4), I5)
    103  format(6(1pe12.4))
    104  format(I5,5(1pe12.4))

    rho0 = mexp / pi4 / (rstar**(3.0d0 - alpha) - 1.d9**(3.0d0 - alpha))
    do i = 1, 2000
        r(i) = rstar / 2.d3 * dble(i)
        if (i==1) then
            rho(i) = rho0 * (0.5d0 * (r(i) + 1.d9))**(-alpha)
        else
            rho(i) = rho0 * (0.5d0 * (r(i) + r(i - 1)))**(-alpha)
        end if
    end do
    r(0) = 0.d0
    mtot = 0.d0
    do i = 1, 2000
        m(i) = pi43 * rho(i) * (r(i)**3 - r(i - 1)**3)
        mtot = mtot + m(i)
    end do
    print *, mtot, mexp
    do i = 1, 2000
        m(i) = m(i) * mexp / mtot
    end do
    !
    !--  The nickel is placed in the center.
    !
    mtot = 0.d0
    do i = 1, 2000
        factor(i) = 1.0
        mtot = mtot + m(i)
        if (mtot<mcore) then
            kappa(i) = 0.2
            jcore = i
        end if

        if (m(i)<mni) then
            edot(i) = m(i) * efac
            mni = mni - m(i)
        else if (mni>0) then
            edot(i) = mni * efac
            mni = 0.d0
        else
            edot(i) = 0.d0
        end if
        v0 = v0 + m(i) * r(i)**2
    end do
    print *, mtot
    print *, 'core zone', jcore
    print *, v0, m(2000) / 1.9d33, r(2000) - r(1999), trecom

   
    v0 = dsqrt(2.0 * eexp / v0)
    etest = 0.d0
    do i = 1, 2000
        v(i) = v0 * r(i)
        rhosh0(i) = ((gam + 1) / (gam - 1)) * rho(i)
        vsh0(i) = v(i) * (rho(i) / rhosh0(i))**(-0.19)
        vshr(i) = v(i) * (rho(i) / rhosh0(i))
        tsh0(i) = (((3 / 2) * (gam + 1) / a) * rhosh0(i) * vshr(i)**2)**0.25
        energy(i) = (a * tsh0(i)**4 * pi43 * (r(i)**3 - r(i - 1)**3))
        if (mod(i, 100)==0) then
            print *, i, r(i), rho(i), tsh0(i), vsh0(i), m(i)
	    write(89, 104) i, r(i), rho(i), tsh0(i), vsh0(i), m(i)
        end if
        etest = etest + 0.5 * m(i) * v(i)**2
    end do
    close(89)
    print *, etest
    print *, 'shock temp ', tsh0(2000)


    

    time = 0.d0
    dt = 1.0
    lumt = 0.d0
    jedge = 2000
    jphot = 2000
    lums = 0.d0
    rhosh = 0.d0
    rhow = 0.d0
    vshw = 0.d0
    accel = 1.0
    rho1 = rho(2000)
    msh = 0.d0
    !**** This accel variable is defined later; from FAST RADIATION MEDIATED SHOCKS AND SUPERNOVA SHOCK BREAKOUTS
    !cccc  Boaz Katz, Ran Budnik, and Eli Waxman

    do i = 1, 10000000
        do j = 1, 2000
            r(j) = r(j) + vsh0(j) * accel * dt
            rho(j) = m(j) / pi43 / (r(j)**3 - r(j - 1)**3)
        end do
        do j = 1, 2000
            if (energy(j)<1.d0) then
                temp(j) = 0.0
            else
                temp(j) = (energy(j) / pi43 / (r(j)**3 - r(j - 1)**3) / a)**0.25
            end if
            if (j>jcore) then
                if (temp(j)<trecom) then
                    factor(j) = (temp(j) / trecom)**precom
                    factor(j) = max(1e-4, factor(j))
                else
                    factor(j) = 1.0
                end if
            end if
            if (j<=jphot - 1) then
                diff = r(j)**2 * difffac / 0.5 / (kappa(j + 1) * factor(j + 1) + &
                        kappa(j) * factor(j)) / &
                        (rho(j + 1) + rho(j)) * &
                        (temp(j + 1)**4 - temp(j)**4) / (r(j + 1) - r(j))
                energy(j + 1) = energy(j + 1) - diff * dt
                energy(j) = energy(j) + diff * dt
            end if
        end do
        do j = 1, 2000
            if (edot(j)>0) then
                energy(j) = energy(j) + edot(j) * dexp(-time / tauni) * dt
            end if
        end do
        tau = 0.d0
        lum = 0.d0

        do j = 1, 2000
            rhosh(j) = ((gam + 1) / (gam - 1)) * rho(j)
            vsh(j) = (rho(j) * vsh0(j) * accel) / rhosh(j)
            tsh(j) = (((3 / 2) * (gam + 1) / a) * rhosh(j) * vsh(j)**2)**0.25
            esh(j) = ((a * tsh(j)**4) * pi43 * (r(j)**3 - r(j - 1)**3))
        end do
        do j = 2000, 1, -1
            tau = tau + rho(j) * kappa(j) * factor(j) * (r(j) - r(j - 1))
            rhow = mdot / (pi4 * r(j)**2 * vwind)
            mw = pi4 * r(j)**2 * tau / kappa(j)

            if (tau>0.4) then
                if (mod(i, 10000)==0) then
                    print *, 'Val', tau, i, jedge, r(j) / 8.d12, temp(j), &
                            energy(j), accel, (rho(j) / rhow)**(0.19)
                end if
                if (tau>10.) then
                    lum = 0
                    msh = rho(j) * pi43 * (rstar**3)
                    accel = ((mw * kappa(j) * tau) / pi4 / &
                            rstar**2)**(0.19)
                    rold = r(j)
                else
                    if (jedge == jphot) then
                        msh = (rho(j) * pi43 * (r(j)**3 - rold**3))
                        rhosh(j) = ((gam + 1) / (gam - 1)) * rho(j)
                        accel = ((msh * kappa(j) * tau) / pi4 / &
                                r(j)**2)**(0.19)
                    end if
                    lum = (sigma * r(j)**2 * temp(j)**4 * dexp(-tau))
                end if
                lum = min(lum, energy(j) * 0.5 / dt)
                lum = max(lum, 0.0)
                lums = lums + lum
                energy(j) = energy(j) - lum * dt
                rphot = r(j)
                tphot = temp(j)
                den = rho(j)
                jphot = j
                vel = vsh0(j) * accel
                goto 20
            else
                jedge = j
                msh = rhow * pi43 * (r(jedge))**3
                lum = sigma * r(j)**2 * temp(j)**4
                lum = min(lum, energy(j) * 0.99 / dt)
                lum = max(lum, 0.d0)
                lums = lums + lum
                energy(j) = energy(j) - lum * dt
                if (energy(j)<1.d10) then
                    energy(j) = 0.d0
                end if
            end if
        end do
        20      continue
        lumt = lumt + lum
        lums = 0.d0
        time = time + dt
        if (mod(i, 10000)==0) then
            print *, 'Lum', tau, time / 3600. / 24., lumt / 10000.d0, esh(2000)
            write(69, 102) time / 3600. / 24., dlog10(lumt / 10000.d0 + 1.d-5), &
                    rphot, tphot, tsh(j), jedge
            write(79, 103) (time / 3600. / 24.), den, tau, vel, msh, rhow
            lumt = 0.d0
        end if
        if (time>200. * 3600. * 24.) then
            stop
        end if
    end do

    close(69)
    close(79)
    

end



! This program is open source under the BSD-3 License.
! Redistribution and use in source and binary forms, with or without modification, are permitted
! provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of c conditions andthe following disclaimer.
! 
! 2.Redistributions in binary form must reproduce the above copyright notice, this list of !c conditions
! and the following disclaimer in the documentation and/or other materials provided with the
! distribution.
! 
! 3.Neither the name of the copyright holder nor the names of its contributors may be used to !c endorse
! or promote products derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
! OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! © (or copyright) 2020. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
! National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
! Department of Energy/National Nuclear Security Administration. All rights in the program are
! reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others acting on its behalf a
! nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
! derivative works, distribute copies to the public, perform publicly and display publicly, and to c permit
! others to do so.
