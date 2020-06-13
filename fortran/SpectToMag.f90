program SpectToMag
implicit none
    !
    double precision temp, avetemp,h,r
    double precision lumr, lumtot, &
            lumop, lumv, lumu, lumb, &
            lumuvw1, lumuvw2, lumuvm2, ene
    double precision rabsmag, vabsmag, &
            babsmag, uabsmag, uvw1absmag,  &
            uvw2absmag, uvm2absmag, absmag
    double precision time(1000), rad(1000), &
            llum(1000),tphot(1000),jedge(1000),&
            shtemp(1000)
    double precision rat, tau, dmbc,lam, c,d
    doubleprecision uvw1mag,uvm2mag,uvw2mag,&
            umag,bmag,vmag,rmag
    doubleprecision lambda(100000), ev(100000),lum
    integer i, j, ijnk, ir, nzones, n
    character*13 filename
    h=4.13d-15 ! planck's constant in ev*s
    c=2.998d18 ! speed of light in angstrom/second
    d=10*(3.08d18) ! Distance of 10 parsecs in cm ( 1 parsec = 3.08d18 cm)
    !...wavelength range from LSST R band to X-Ray...!
    do n=1, 100000
        lambda(n)=0.5+ 0.10 * dble(n)
        ev(n)=(h * c) / lambda(n)
    end do
    !...Read in BB filename and number of lines...!
    read(*,*) filename,nzones
    open(42, file=filename//'BB.dat')
      do j=1, nzones
        read(42, *) time(j), llum(j), rad(j), tphot(j), shtemp(j), &
                jedge(j)
      end do
    close(42)
    !...Create output file for Absolute Magntidues...!
    open(43, file=filename//'Absmag.dat')

      do j=1, nzones
    !...Generates spectra for each modeled epoch...!
        if (jedge(j) .eq. 2000.0) then
    !...Blackbody properties during shock breakout...!
    !...For additonal shock heating and xray spect set temp = shtemp(j)/1.16d4...!
               temp=tphot(j)/1.16d4
               r=rad(j)
               lum=llum(j)
        	else
    !...Data Smoothing: Use average BB photospheric temperature and radius beyond shock breakout...!
                avetemp=sum(tphot(j-4:j+4))/9.0
                r = sum(rad(j-4:j+4))/9.0
                lum = llum(j)
                temp=avetemp / 1.16d4 !1.16d4 is temp K to eV
            end if
        lumtot=0.d0
        lumop=0.d0
        lumr=0.d0
        lumv=0.d0
        lumu=0.d0
        lumb=0.d0
        lumuvw1=0.d0
        lumuvw2=0.d0
        lumuvm2=0.d0
    !...Spectral Brightness of BB radiation as a function of wavelength for each temperature defined above...!
            do i=1, 100000
                ene=ev(i)! hc/wavenlength in eV
                lam=lambda(i)! wavelength in Angstrom
             if (ene / temp<300) then
                rat=ene**3 / (dexp(ene / temp) - 1.d0)! Derivation of Planck's Law
                dmbc = -5.0*(dlog10(d**2./r**2.))+lum !Distance Mod and intrinstic bolometric luminosity
                else
                    rat=0.
                    dmbc=0.
            end if
            lumtot=lumtot + rat

              if (ene>2.60.and.ene<7.75) then
               !FWHM energy range for UVW1
                lumuvw1=lumuvw1 + rat !spectral distribution for given wavelength and temperature
            end if

             if (ene>4.11.and.ene.lt.7.75) then
                !FWHM energy range for UVM2
                 lumuvm2=lumuvm2 + rat!spectral distribution for given wavelength and temperature
            end if

            if (ene>3.50.and.ene<7.75) then
                !FWHM energy range for UVW2
                lumuvw2=lumuvw2 + rat!spectral distribution for given wavelength and temperature
            end if

            if (ene>2.06.and.ene<2.53) then
                !FWHM energy range for V
                lumv=lumv + rat!spectral distribution for given wavelength and temperature
            end if

            if (ene>2.55.and.ene<3.00) then
                !FWHM energy range for B
                lumb=lumb + rat!spectral distribution for given wavelength and temperature
            end if

            if (ene>3.50.and.ene<4.11) then
                !FWHM energy range for U
                lumu=lumu + rat!spectral distribution for given wavelength and temperature
            end if

            if (ene>1.80.and.ene<2.23) then
                !FWHM energy range for R
                lumr=lumr + rat!spectral distribution for given wavelength and temperature
            end if

            end do
            rabsmag=-2.5*(log10(lumr))- dmbc !Modeled Absolute R Magnitude
            call readUVOTresponse(lumuvw1, lumuvm2,lumuvw2,&
                lumu,lumb,lumv,dmbc, uvw1mag,uvm2mag,uvw2mag,&
                umag,bmag,vmag,lam)
            call readLSSTresponse(lumr,dmbc,lam,rmag)
              write(43, 102) time(j), dmbc, uvw1mag, &
                    uvm2mag,uvw2mag,umag,bmag,&
                      vmag,rmag,rabsmag
            end do
102  format(10(1pe12.4))
    close(43)

end program SpectToMag
!...Read in UVOT Reponse Curve...!
subroutine readUVOTresponse(lumuvw1, lumuvm2,lumuvw2,&
        lumu,lumb,lumv,dmbc, uvw1mag,uvm2mag,uvw2mag,&
        umag,bmag,vmag,lam)
    doubleprecision, intent(in):: lumuvw1,lumuvm2,lumuvw2
    doubleprecision, intent(in):: lumu,lumb,lumv
    doubleprecision, intent(in):: dmbc,lam
    doubleprecision, intent(out):: uvw1mag,uvm2mag,uvw2mag
    doubleprecision, intent(out):: umag,bmag,vmag
    doubleprecision wavelength, w1r,m2r,w2r,ur,br,vr
    doubleprecision uvw1,uvm2,uvw2,u,b,v
    integer i
    open(45, file='C:\Users\Janie\PycharmProjects\SNeModels\UVOTResponseCurve')
    uvw1 = 0.0
    uvm2 = 0.0
    uvw2 = 0.0
    u = 0.0
    b = 0.0
    v = 0.0
    do i=1, 640
        ! Read in Swift effective area transmission per bandpass filter
        read(45,*) wavelength, w1r,m2r,w2r,ur,br,vr
        !filter spectral distrubution to UVOT response
        uvw1 = uvw1+ lumuvw1*(w1r/30.0)
        uvm2 = uvm2+ lumuvm2*(m2r/30.0)
        uvw2 = uvw2+ lumuvw2*(w2r/30.0)
        u = u+ lumu*(ur/60.0)
        b = b+ lumb*(br/60.0)
        v = v+ lumv*(vr/60.0)
    end do
    close (45)
    !convert filtered spectral distribution to Swift UVOT magnitudes
    uvw1mag=-2.5*(log10(uvw1))- dmbc
    uvm2mag=-2.5*(log10(uvm2))- dmbc
    uvw2mag=-2.5*(log10(uvw2))- dmbc
    umag=-2.5*(log10(u))- dmbc
    bmag=-2.5*(log10(b))- dmbc
    vmag=-2.5*(log10(v))- dmbc
    return
end subroutine readUVOTresponse
!...Read in LSST Throughput Response...!
subroutine readLSSTresponse(lumr,dmbc,lam,rmag)
    doubleprecision, intent(in):: lumr
    doubleprecision, intent(in):: dmbc,lam
    doubleprecision, intent(out):: rmag
    doubleprecision wavelength, rr,r,response
    integer i
     ! Read in LSST R band throughput
    open(45, file='C:\Users\Janie\PycharmProjects\SNeModels\LSST_Rbandpass.txt')
    r = 0.0
    do i=1, 1688
        read(45,*) wavelength, rr,response
        !filter spectral distrubution to LSST R band response
        r = r + (lumr*response)
    end do
    close (45)
    !convert filtered spectral distribution to LSST R magnitudes
    rmag=-2.5*(log10(r))-dmbc+8.0
    return
end subroutine readLSSTresponse

c This program is open source under the BSD-3 License.
c Redistribution and use in source and binary forms, with or without modification, are permitted
c provided that the following conditions are met:
c 1. Redistributions of source code must retain the above copyright notice, this list of c conditions andthe following disclaimer.
c 
c 2.Redistributions in binary form must reproduce the above copyright notice, this list of !c conditions
c and the following disclaimer in the documentation and/or other materials provided with the
c distribution.
c 
c 3.Neither the name of the copyright holder nor the names of its contributors may be used to !c endorse
c or promote products derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
c IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
c PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
c CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

c © (or copyright) 2020. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
c National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
c Department of Energy/National Nuclear Security Administration. All rights in the program are
c reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
c Security Administration. The Government is granted for itself and others acting on its behalf a
c nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
c derivative works, distribute copies to the public, perform publicly and display publicly, and to c permit
c others to do so.