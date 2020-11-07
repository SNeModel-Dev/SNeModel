program spec
    !
    implicit none
    !
    double precision temp, tempdc, h, lam
    double precision lumr, lumtot, &
            lumop, lumv, lumu, lumb, &
            lumuvw1, lumuvw2, lumuvm2, ene
    double precision time, djnk, llum, den, c
    double precision rat, taux, jedge, tau, tphot
    !      double precision taux,trecom
    double precision lambda(8000), ev(8000)
    integer i, j, ijnk, ir, nzones, n
    character*32 filename
    print *, 'about to read filename and zones'
    read(*,*) filename,nzones
    print *, "read in filename"
    open(42, file=filename//'BB.dat')
    open(43, file=filename//'spect.dat')
    h=4.13d-15
    c=2.998d18
    ! lam=100.25
    ! Changing lam to match previous verions
    lam=10.25
    print *,"About to build lambdas"
    do n=1, 8000
        lambda(n)=lam + 0.75 * dble(n)
        ev(n)=(h * c) / lambda(n)
    end do
    print *,"Finished building lambdas"
    print *, lambda(8000), lambda(1), ev(1), ev(8000)
    print *,"Going to iterate over the zones"
    do j=1, nzones
        read(42, *) time, llum, djnk, tphot, tempdc, &
                jedge

        !         if (time.gt.0.9.and.llum.gt.48.) read *, ijnk
        !            tempdc=(7.5646d-15*temp**4)/(1.38d-16*6.022d23*den)
        !            temp=max(temp,tempdc)*3.0*8.6d-5
        !			if file name tshmsh then use shocking heating temperature during SBO
        !           if (jedge .eq. 2000.0) then
        !               temp=tempdc/1.16d4
        !			else
        temp=tphot / 1.16d4
        !            endif
        ! 1.16d4 is temp K to eV
        lumtot=0.d0
        lumop=0.d0
        lumr=0.d0
        lumv=0.d0
        lumu=0.d0
        lumb=0.d0
        lumuvw1=0.d0
        lumuvw2=0.d0
        lumuvm2=0.d0
        taux=0.0
        do i=1, 8000
            ene=ev(i)
            if (ene / temp.lt.300) then
                rat=ene**3 / (dexp(ene / temp) - 1.d0)
            else
                rat=0.
            end if
            lumtot=lumtot + rat
            if (ene.gt.2.60.and.ene.lt.7.75) then
                lumuvw1=lumuvw1 + rat
                !*dexp(-taux)
            end if
            if (ene.gt.4.11.and.ene.lt.7.51) then
                lumuvm2=lumuvm2 + rat
                !*dexp(-taux)
            end if
            if (ene.gt.3.50.and.ene.lt.7.75) then
                lumuvw2=lumuvw2 + rat

            end if
            if (ene.gt.2.06.and.ene.lt.2.53) then
                lumv=lumv + rat
            end if
            if (ene.gt.2.47.and.ene.lt.3.41) then
                lumb=lumb + rat
            end if
            if (ene.gt.2.60.and.ene.lt.4.11) then
                lumu=lumu + rat
            end if
            if (ene.gt.1.77.and.ene.lt.2.60) then
                lumop=lumop + rat
            endif
            if (ene.gt.1.77.and.ene.lt.2.06) then
                lumr=lumr + rat
                !*(1.d0-dexp(-tau))
            end if
            !               ene=ene+0.1
        end do
        !            lumv=lumv+0.25*sub
        !            lumb=lumb+0.25*sub
        if (lumr.gt.1.d-10) then
            write(43, 102) time, llum, &
                    llum+ dlog10(lumr / lumtot), &
                    llum + dlog10(lumop / lumtot), &
                    llum + dlog10(lumv / lumtot), &
                    llum + dlog10(lumb / lumtot), &
                    llum + dlog10(lumu / lumtot), &
                    llum + dlog10(lumuvw1 / lumtot), &
                    llum + dlog10(lumuvm2 / lumtot), &
                    llum + dlog10(lumuvw2 / lumtot)
        else
            write(43, 102) time, llum , &
                    llum+ dlog10(lumr / lumtot),&
                    llum + dlog10(lumop / lumtot), &
                    llum + dlog10(lumv / lumtot), &
                    llum + dlog10(lumb / lumtot), &
                    llum + dlog10(lumu / lumtot), &
                    llum + dlog10(lumuvw1 / lumtot), &
                    llum + dlog10(lumuvm2 / lumtot), &
                    llum + dlog10(lumuvw2 / lumtot)
        end if
    end do
    close(42)
    close(43)
    102  format(10(1pe12.4))
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
