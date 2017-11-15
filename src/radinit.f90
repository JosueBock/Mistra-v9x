!
! Copyright 1996-2017 the Authors
!
! Licensed under the EUPL, Version 1.1 only (the "Licence");
!
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!   https://joinup.ec.europa.eu/software/page/eupl
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an
! "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
! either express or implied.
!
! See the Licence for the specific language governing permissions
! and limitations under the Licence.




! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! radinit.f: radiation initialisation and start (SR str)
!
! contains the following subroutines:
!     - intrad
!     - initstr
!     - ipdata
!     - initr
!     - load1
!     - str
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine intrad

! Description :
! -----------
!    Linear interpolation of radiation parameters qabs0, qext0, asym0
!    which are tabulated for droplet radii xw0 (in microns)
!    and percentage part of pure water xa0 on the refractive index m
!      with: m = ma + (mw -ma) *xa0,
!      mw: refractive index pure water, ma: refractive index pure aerosol
!    on the radiation parameters qabs, qext and asym
!    for actual droplet radii xw1 (in microns) and percentage part of
!    pure water xa1 on the refractive index.
!
!    More information about the calculations of the tabulated values qabs0, qext0, asym0
!    using Mie theorie can be found in the following paper:
!    Bott, Sievers & Zdunkowski, J. Atmos. Sci., 47 (18), 2153-2166, 1990
!    Note that in this paper, alpha = 1-xa


! Interface :
! ---------
!    SR intrad is called by main program during initialisation

! Input :
! -----
!    - data files (3 aerosol types, short (visible) and long (IR) bands)
!    - rn (dry particle radius) and rq (total particle radius) from /cb50/

! Output :
! ------
!    - interpolated values qabs, qext, asym in /cb49/

! Externals :
! ---------
!    none


! Method :
! ------
!    Find index iw0 (2 <= iw0 <= nw0) such that xw0(iw0-1) < xw1 <= xw0(iw0)
!    Then define an interpolation factor dx (proportionality factor, 0 <= dx <= 1)
!      such that the interpolated value is dx*f(iw0) + (1-dx)*f(iw0-1) where f is one of the
!      radiation parameter qabs, qext, asym.

!    If xw1 < xw0(1) or xw1 > xw0(nw0), then the boundary value (xw0(1) or xw0(nw0)) will be
!       used without interpolation nor extrapolation (dx=0. or dx=1., respectively). Warning
!       messages are displayed in these specific cases.  

!    Similarly, find index ia0 (2 <= ia0 <= na0) such that xa0(ia0-1) < xw1 <= xa0(ia0), and define an
!       interpolation factor dy.
!    Note the different treatment of special cases:
!       - if xa1 < 0, an error message arise and the program is stopped
!       - if xa1 > 1, another error message arise and the program is stopped


! Author :
! ------
!    Andreas Bott


! Modifications :
! -------------
!     ?         Roland       Introduction of special cases to avoid out-of-bounds errors, and deal
!               von Glasow      with particle radius smaller / larger than tabulated values
!
!    Jul-2016   Josue Bock   Header with comments
!                            Use module for parameters
!                            Explicit declaration of all variables, implicit none
!                            Reindexed qabs/qext/asym for efficiency (innermost index is leftmost)
!                            Removed labeled do loops
!                            Rewritten, with more generic cases, error/warning messages, comments and cleaned code structure
!
!    Oct-2016   Josue Bock   Further commenting after personnal communication with A. Bott
!
! 10-Mar-2017   Josue Bock   Improvement of j case (now accept xa1=0, error if xa<0)
!                            Cosmetic reorganisation, new header organisation
!                            Set nkaer as a new global parameter (later replaced by jptaerrad)
!
! 04-Apr-2017   Josue Bock   Fortran90 conversion
!
! 05-Nov-2017   Josue Bock   Removed module directories, replaced by config
!  
! 09-Nov-2017   Josue Bock   Final cleaning after second F90 conversion (GitHub version),
!                              rewritten interpolation, coding standards for variable names

! == End of header =============================================================


! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       cinpdir

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfunerr,      &
       jpfunout,      &
       jpfunaerrad

  USE global_params, ONLY : &
! Imported Parameters:
       jptaerrad,           &
       nka,                 &
       nkt,                 &
       mb,                  &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local parameters:
  integer, parameter :: na0 = 11  ! Number of tabulated values, volume fraction of water
  integer, parameter :: nw0 = 40  ! Number of tabulated values, total aerosol radius

  real (kind=dp), parameter :: xa0(na0) = &  ! volume fraction of pure water (tabulated values) []
       (/0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.7_dp, 0.8_dp, 0.85_dp, 0.9_dp, 0.95_dp, 0.975_dp, 1.0_dp/)
  real (kind=dp), parameter :: xw0(nw0) = &  ! total radius of the scattering spheres (tabulated values) [micro m]
       (/0.01_dp, 0.0125_dp, 0.015_dp, 0.02_dp, 0.025_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp, 0.08_dp, &
         0.1_dp,  0.125_dp,  0.15_dp,  0.2_dp,  0.25_dp,  0.3_dp,  0.4_dp,  0.5_dp,  0.6_dp,  0.8_dp,  &
         1.0_dp,  1.25_dp,   1.5_dp,   2.0_dp,  2.5_dp,   3.0_dp,  4.0_dp,  5.0_dp,  6.0_dp,  8.0_dp,  &
        10.0_dp, 12.5_dp,   15.0_dp,  20.0_dp, 25.0_dp,  30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp, 80.0_dp/)
  

! Local scalars:
  character (len=len_trim(cinpdir)+11) :: fname
  real (kind=dp) :: dx,dy,dxdy,dxdy1,dx1dy,dx1dy1 ! Interpolation coefficients
  real (kind=dp) :: xa1                           ! Actual percentage of water computed for each particle class (2D spectrum)
  real (kind=dp) :: xw1                           ! Actual total aerosol radius [micro m] for each particle class (2D spectrum)
  integer :: ifun                                 ! File unit number
  integer :: ia0,iw0                              ! upper indexes of tabulated values for interpolation
  integer :: ja0,jw0                              ! loop indexes, tabulated values
  integer :: jaer                                 ! loop index for aerosol type
  integer :: jb                                   ! loop index for spectral band number
  integer :: jka,jkt                              ! loop indexes, 2D microphysical grid

! Local arrays:
  real (kind=dp) :: qabs0(mb,nw0,na0) ! tabulated values of absorption coefficient
  real (kind=dp) :: qext0(mb,nw0,na0) ! tabulated values of extinction coefficient
  real (kind=dp) :: asym0(mb,nw0,na0) ! tabulated values of asymmetry factor


! Common blocks:
  common /cb49/ qabs(mb,nkt,nka,jptaerrad),qext(mb,nkt,nka,jptaerrad),asym(mb,nkt,nka,jptaerrad)
  real (kind=dp) :: qabs,qext,asym

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)                      ! only rn (dry radius), rq (total radius) are used here
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

! == End of declarations =======================================================


! Open input data files
! ---------------------

! first unit number to read
  ifun=jpfunaerrad

! input qabs, qext and asym for radiation parameters of particles
!   *kw.dat: short wave data
!   *lw.dat: long wave data
  fname=TRIM(cinpdir)//"urbankw.dat"
  open (unit=ifun, file=fname, status='old')
  fname=TRIM(cinpdir)//"urbanlw.dat"
  open (unit=ifun+1, file=fname, status='old')
  fname=TRIM(cinpdir)//"ruralkw.dat"
  open (unit=ifun+2, file=fname, status='old')
  fname=TRIM(cinpdir)//"rurallw.dat"
  open (unit=ifun+3, file=fname, status='old')
  fname=TRIM(cinpdir)//"ozeankw.dat"
  open (unit=ifun+4, file=fname, status='old')
  fname=TRIM(cinpdir)//"ozeanlw.dat"
  open (unit=ifun+5, file=fname, status='old')

! input format to read these files
5000 format (3e16.8)

  do jaer=1,jptaerrad  ! jaer=1 urban, jaer=2 rural, jaer=3 ocean

! Read input data file: tabulated values
! --------------------------------------
     do ja0=1,na0
        do jw0=1,nw0
           do jb=1,mbs
              read(ifun,5000) qabs0(jb,jw0,ja0),qext0(jb,jw0,ja0),asym0(jb,jw0,ja0)
           enddo
        enddo
     enddo
     close (ifun)
     ifun=ifun+1
     do ja0=1,na0
        do jw0=1,nw0
           do jb=mbs+1,mb
              read(ifun,5000) qabs0(jb,jw0,ja0),qext0(jb,jw0,ja0),asym0(jb,jw0,ja0)
           enddo
        enddo
     enddo
     close (ifun)
     ifun=ifun+1

! Loop over the aerosol class
! ---------------------------
     do jka=1,nka
        do jkt=1,nkt
           ! Values of total radius (xw1) and volume mixing ratio of water (xa1) for each aerosol class
           ! ------------------------------------------------------------------------------------------
           xw1=rq(jkt,jka)
           xa1=1._dp-(rn(jka)/rq(jkt,jka))**3

           ! Locate xw1 in the tabulated xw0 grid, and calculate interpolation factor dx
           ! ---------------------------------------------------------------------------
           if (xw1 < xw0(1)) then
              ! case rq < xw0(1), warn user, set iw0 and dx values as if xw1=xw0(1)
              iw0 = 2
              dx  = 0._dp
              write(jpfunout,*)'Warning: in SR intrad, rq < xw0(1)',jka,jkt,xw1
           else if (xw1 > xw0(nw0)) then
              ! case rq > xw0(nw0), warn user, set iw0 and dx values as if xw1=xw0(nw0)
              iw0 = nw0
              dx  = 1._dp
              write(jpfunout,*)'Warning: in SR intrad, rq > xw0(nw0)',jka,jkt,xw1
           else
              ! general case: xw0(iw-1) < rq <= xw0(iw)
              iw0 = 2
              do while (xw1 > xw0(iw0))
                 iw0 = iw0 + 1
              end do
              dx=(xw1-xw0(iw0-1))/(xw0(iw0)-xw0(iw0-1))
           end if

           ! Locate xa1 in the tabulated xa0 grid, and calculate interpolation factor dy
           ! ---------------------------------------------------------------------------
           ! First test to detect potential model inconsistency
           !  (this case should never happend: if the particle grid has been properly defined,
           !   rq(jt,ia) >= rn(ia) thus xa1 >= 0.)
           if (xa1 < xa0(1)) then
              write(jpfunerr,*)'Error in SR intrad: xa1 < xa0(1)=0.',jka,jkt,xa1
              stop ' Stopped by SR intrad'
           else if (xa1 > xa0(na0)) then
              write(jpfunerr,*)'Error in SR intrad: xa1 > xa0(na0)=1.',jka,jkt,xa1
              stop ' Stopped by SR intrad'
           else
              ! general case: xa0(ia-1) < xa1 <= xa0(ia)
              ia0 = 2
              do while (xa1 > xa0(ia0))
                 ia0 = ia0 + 1
              end do
              dy=(xa1-xa0(ia0-1))/(xa0(ia0)-xa0(ia0-1))
           end if


           ! Interpolation factors
           ! ---------------------
           dxdy   = dx*dy
           dxdy1  = dx*(1._dp-dy)
           dx1dy  = (1._dp-dx)*dy
           dx1dy1 = (1._dp-dx)*(1._dp-dy)

           ! Interpolated values of qabs, qext and asym
           ! ------------------------------------------
           do jb=1,mb
              qabs(jb,jkt,jka,jaer) = dxdy   * qabs0(jb,iw0,ia0)   &
                                    + dxdy1  * qabs0(jb,iw0,ia0-1) &
                                    + dx1dy  * qabs0(jb,iw0-1,ia0)   &
                                    + dx1dy1 * qabs0(jb,iw0-1,ia0-1)

              qext(jb,jkt,jka,jaer) = dxdy   * qext0(jb,iw0,ia0)   &
                                    + dxdy1  * qext0(jb,iw0,ia0-1) &
                                    + dx1dy  * qext0(jb,iw0-1,ia0)   &
                                    + dx1dy1 * qext0(jb,iw0-1,ia0-1)

              asym(jb,jkt,jka,jaer) = dxdy   * asym0(jb,iw0,ia0)   &
                                    + dxdy1  * asym0(jb,iw0,ia0-1) &
                                    + dx1dy  * asym0(jb,iw0-1,ia0)   &
                                    + dx1dy1 * asym0(jb,iw0-1,ia0-1)
           enddo

        end do ! jkt loop
     end do ! jka loop

  end do ! jaer loop

end subroutine intrad
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine initstr
!
! Description :
! -----------
!    first calculation of radiative fluxes and heating rates


! Interface :
! ---------
!    SR initstr is called by main program during initialisation

! Input :
! -----

! Output :
! ------
!    - /cb15/ fnseb,flgeg
!
!
! History:
! Version   Date     Comment
! -------   ----     -------
  ! 10-Nov-2017   Josue Bock   Removal of cb55, which is useless
  !                            Calling nstrahl twice during initialisation has absolutely no effect,
  !                              since the only apparent reason was to fill cb55 data, which were not
  !                              actually used. This simplifies much this subroutine.
  !
! 1.1      07/2016   Removal of labeled do loops.                 <Josue Bock>
!                    Use module for parameters
!                    All explicit declarations and implicit none
!
! 1.0       ?        Original code.                               <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

  implicit none

! Local scalars:
  logical :: linit

!- End of header ---------------------------------------------------------------

! initialisation of radiation code
! --------------------------------
  call ipdata
  call initr

  linit = .true.
  call rotate_in(linit)

  call nstrahl
  call profr

  call rotate_out

end subroutine initstr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine ipdata
!
! Description :
! -----------
!    input data for radiation code:
!    reads from file 'pifm2_yymmdd.dat'
!    writes in several common blocks
!    descriptions of variables can be found at the end of file 'pifm2_yymmdd.dat'
!    (part of these description written in German)


! Modifications :
! -------------
  !    Jul-2016  Josue Bock  Header, comments
  !                          Use modules for parameters
  !
  !    Oct-2016  Josue Bock  More cleaning, after checking the variables actually used
  !                            (several old variables from PIFM1, apparently)
  !
  ! 28-Oct-2017  Josue Bock  Reindexed pibtab (thus changed input file) so that it is read
  !                            in optimum Fortran order (leftmost is innermost in loops)
  !                          Replaced dimension 12 by mbir
  !
  ! 15-Nov-2017  Josue Bock  Merged input data from pifm2*.dat and initr*.dat
  !                          Cleaned some unused data from PIFM1 (mostly in cb56)
  !                          Adopted the same data file structure as WS did in his Mistra-F90 version
  !                            ==> all comment lines are from his version
  !                          Each data block is now separated by a comment line.
  !                          Also added a former common block /sol/ to remove hardcoded values from
  !                            the radiative code.
  !                          Final cleaning after Fortran90 conversion

! == End of header =============================================================

! Declarations :
! ------------
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       cinpdir,      &
! Imported Routines:
       abortM

  USE file_unit, ONLY : &
! Imported Parameters:
       jpfundatarad

  USE global_params, ONLY : &
! Imported Parameters:
       mb,                  &
       mbs,                 &
       mbir,                &
       ncw,                 &
       nrlay

  USE precision, ONLY : &
       dp

  implicit none

! Local scalars:
  character (len=len_trim(cinpdir)+16) :: clfname
  integer :: istat

! Common blocks:
  common /aeag/ feux(8),seanew(8,mb,4),saanew(8,mb,4),ganew(8,mb,4)
  real (kind=dp) :: feux, seanew, saanew, ganew
  common /band1/ hk1(10),fk1o3(10)
  real (kind=dp) :: hk1, fk1o3
  common /band2/ hk2(8),c2h2o(3,11,8)
  real (kind=dp) :: hk2, c2h2o
  common /band3/ hk3(12),c3h2o(3,11,12)
  real (kind=dp) :: hk3, c3h2o
  common /band4/ hk4(7),c4h2o(3,11,7)
  real (kind=dp) :: hk4, c4h2o
  common /band5/ hk5(12),c5h2o(3,11,12)
  real (kind=dp) :: hk5, c5h2o
  common /band6/ hk6(5),c6h2o(3,11,5)
  real (kind=dp) :: hk6, c6h2o
  common /band7/ hk7(2),c7h2o(3,19,2)
  real (kind=dp) :: hk7, c7h2o
  common /band8/ hk8(3),c8h2o(3,19,3)
  real (kind=dp) :: hk8, c8h2o
  common /band9/ hk9(4),c9h2o(3,19,4)
  real (kind=dp) :: hk9, c9h2o
  common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
  real (kind=dp) :: hk10, c10h2o, c10ch4, c10n2o
  common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
  real (kind=dp) :: hk11, c11h2o, c11ch4, c11n2o
  common /band12/ hk12(5),c12o3(3,19,5),c12h2o(3,19)
  real (kind=dp) :: hk12, c12o3, c12h2o
  common /band13/ hk13(2),c13h2o(3,19,2)
  real (kind=dp) :: hk13, c13h2o
  common /band14/ hk14(10),c14hca(3,19,10),c14hcb(3,19,10)
  real (kind=dp) :: hk14, c14hca, c14hcb
  common /band15/ hk15(12),c15hca(3,19,12),c15hcb(3,19,12)
  real (kind=dp) :: hk15, c15hca, c15hcb
  common /band16/ hk16(7),c16h2o(3,19,7)
  real (kind=dp) :: hk16, c16h2o
  common /band17/ hk17(7),c17h2o(3,19,7)
  real (kind=dp) :: hk17, c17h2o
  common /band18/ hk18(8),c18h2o(3,19,8)
  real (kind=dp) :: hk18, c18h2o

  common /cb19/ berayl(mbs),bea(mb,nrlay),baa(mb,nrlay),ga(mb,nrlay)
  real (kind=dp) :: berayl, bea, baa, ga

  common /cb56/ o3un(52)
  real (kind=dp) :: o3un

  common /sol/ s0b(mbs), s0tot
  real (kind=dp) :: s0b, s0tot

  common /plancd/ ttab(35),pibtab(35,mbir) ! jjb was used in FN fst4 (called by SR plancktab), both unused now
  real (kind=dp) :: ttab, pibtab           !     left for backward compatibility, might be deleted in the future

  common /was1/ ret(ncw),r2wt(ncw),b2wt(ncw,mb),w2wt(ncw,mb),g2wt(ncw,mb)
  real (kind=dp) :: ret, r2wt, b2wt, w2wt, g2wt

! == End of declarations =======================================================

  clfname=TRIM(cinpdir)//'pifm2_171115.dat'

  open(unit=jpfundatarad, file=clfname, status='old', iostat=istat)
  if (istat /= 0) call abortM ('Error in SR ipdata: cannot open radiation data file: '//clfname)

! input formats to read this file
10 format(a)
11 format(8e16.8)

! planck table                 ! could be deleted if fst4 and plancktab are deleted
  read(jpfundatarad,10)
  read(jpfundatarad,11) ttab
  read(jpfundatarad,10)
  read(jpfundatarad,11) pibtab
! effective droplet radii and radiation coefficients for water
  read(jpfundatarad,10)
  read(jpfundatarad,11) ret
  read(jpfundatarad,10)
  read(jpfundatarad,11) r2wt
  read(jpfundatarad,10)
  read(jpfundatarad,11) b2wt
  read(jpfundatarad,10)
  read(jpfundatarad,11) w2wt
  read(jpfundatarad,10)
  read(jpfundatarad,11) g2wt
! table of reference relative humidities and tabulated radiation coefficients for aerosols
! last index of seanew, saanew, ganew = aerosol type: 1=rural, 2=urban, 3=ocean, 4=background
  read(jpfundatarad,10)
  read(jpfundatarad,11) feux
  read(jpfundatarad,10)
  read(jpfundatarad,11) seanew
  read(jpfundatarad,10)
  read(jpfundatarad,11) saanew
  read(jpfundatarad,10)
  read(jpfundatarad,11) ganew
! solar energy in the 6 SW spectral bands
  read(jpfundatarad,10)
  read(jpfundatarad,11) s0b
! radiation coefficients for gas absorption for every spectral band
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk1
  read(jpfundatarad,10)
  read(jpfundatarad,11) fk1o3
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk2
  read(jpfundatarad,10)
  read(jpfundatarad,11) c2h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk3
  read(jpfundatarad,10)
  read(jpfundatarad,11) c3h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk4
  read(jpfundatarad,10)
  read(jpfundatarad,11) c4h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk5
  read(jpfundatarad,10)
  read(jpfundatarad,11) c5h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk6
  read(jpfundatarad,10)
  read(jpfundatarad,11) c6h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk7
  read(jpfundatarad,10)
  read(jpfundatarad,11) c7h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk8
  read(jpfundatarad,10)
  read(jpfundatarad,11) c8h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk9
  read(jpfundatarad,10)
  read(jpfundatarad,11) c9h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk10
  read(jpfundatarad,10)
  read(jpfundatarad,11) c10h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) c10ch4
  read(jpfundatarad,10)
  read(jpfundatarad,11) c10n2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk11
  read(jpfundatarad,10)
  read(jpfundatarad,11) c11h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) c11ch4
  read(jpfundatarad,10)
  read(jpfundatarad,11) c11n2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk12
  read(jpfundatarad,10)
  read(jpfundatarad,11) c12o3
  read(jpfundatarad,10)
  read(jpfundatarad,11) c12h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk13
  read(jpfundatarad,10)
  read(jpfundatarad,11) c13h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk14
  read(jpfundatarad,10)
  read(jpfundatarad,11) c14hca
  read(jpfundatarad,10)
  read(jpfundatarad,11) c14hcb
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk15
  read(jpfundatarad,10)
  read(jpfundatarad,11) c15hca
  read(jpfundatarad,10)
  read(jpfundatarad,11) c15hcb
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk16
  read(jpfundatarad,10)
  read(jpfundatarad,11) c16h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk17
  read(jpfundatarad,10)
  read(jpfundatarad,11) c17h2o
  read(jpfundatarad,10)
  read(jpfundatarad,11) hk18
  read(jpfundatarad,10)
  read(jpfundatarad,11) c18h2o
! unreduced ozone amount from Craig table
  read(jpfundatarad,10)
  read(jpfundatarad,11) o3un
! coefficients for rayleigh scattering
  read(jpfundatarad,10)
  read(jpfundatarad,11) berayl

  close(jpfundatarad)

end subroutine ipdata
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine initr
!
! Description:
!    standard atmosphere between etw(n) and 50000 meters.
!    input of constant radiation parameters
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
!        14/11/2017 <jjb> removed ntypdx(nrlay) (droplet type), initialised =4 but unused
!                         removed ntypax(nrlay) (aerosol type) from cb01, local declaration
!        10/03/2017 <jjb> corrected rf=0.08 for level 30000 (was 0.8)
!
! 1.2       10/2016  corrected several array size
!                    initialisation of qmo3(nrlev) was missing
!                    initialisation of rnaer(nrlay) was missing
!                    removed rnaer from /cb02/, turned into a local array
!
!           07/2016  Removal of labeled do loops.                        <Josue Bock>
!                    Use modules for parameters
!                    data file updated with berayl(6)
!                      (previously, an old version berayl(4) was read here, then
!                       berayl was redefined with 'data' in strahl). But for some
!                       unknown reason, only the last two values were updated,
!                       the first four were NOT overwritten, while they were expected to be)
!
! 1.1       ?        Calculation of qmo3 for layers 1-80 (bug fix)       <Roland von Glasow>
!
! 1.0       ?        Original code.                                      <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

  USE constants, ONLY : &
! Imported Parameters:
       r0               ! Specific gas constant of dry air, in J/(kg.K)

  USE global_params, ONLY : &
! Imported Parameters:
       n,                   &
       nrlay,               &
       nrlev,               &
       nka,                 &
       mb,                  &
       mbs,                 &
       mbir

  USE precision, ONLY : &
! Imported parameters:
       dp

  implicit none

! Local scalars:
  integer :: i, j, jj, jp, k, na, naf, nafm, nw, jb
  real (kind=dp) :: gamma                      ! Vertical temperature gradients [K/m]
  real (kind=dp) :: rf, zj, zjp
  real (kind=dp) :: dd, emdd, xdd, xemdd, xnaer
  real (kind=dp) :: dz                         ! Height increment of extra layers up to the tropopause [m]

! Local arrays:
  integer :: ntypax(nrlay)
  real (kind=dp) :: etax(nrlev), usav(nrlay)
  real (kind=dp) :: rnaer(nrlay)               ! total number concentration of aerosol particles [cm**-3]

! Internal function:
  real (kind=dp) :: p21, tt

! Common blocks:
  common /aeag/ feux(8),seanew(8,mb,4),saanew(8,mb,4),ganew(8,mb,4)
  real (kind=dp) :: feux, seanew, saanew, ganew

  common /cb01/ tx(nrlev),px(nrlev),rhox(nrlev),xm1x(nrlev),      &
                rho2wx(nrlay),fracx(nrlay),zx(nrlev),thkx(nrlay), &
                beax(mb,nrlay),baax(mb,nrlay),gax(mb,nrlay),      &
                qmo3x(nrlev),rewx(nrlay), tsx
  real (kind=dp) :: tx,px,rhox,xm1x,rho2wx,fracx,zx,thkx, &
                    beax, baax, gax, qmo3x, rewx, tsx

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb19/ berayl(6),bea(mb,nrlay),baa(mb,nrlay),ga(mb,nrlay)
  real (kind=dp) :: berayl, bea, baa, ga

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks, &
                bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
  real (kind=dp) :: g,a0m,b0m,ug,vg,z0,ebs,psis,aks, &
                    bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

  common /cb56/ o3un(52)
  real (kind=dp) :: o3un

  common /sol/ s0b(mbs), s0tot
  real (kind=dp) :: s0b, s0tot

  common /tmp2/ as(mbs),ee(mbir)
  real (kind=dp) :: as, ee

!- End of header ---------------------------------------------------------------

  p21(tt)=610.7*exp(17.15*(tt-273.15)/(tt-38.33))


! albedo of the ground for the six solar wavelength regions of the radiation code
  albedo(:)=0.05
!  albedo(:)=0.8       ! snow
  as = albedo
! emissivity of the ground
  ee(:)=1.0

! total solar energy = sum of solar energies for each spectral bands
  s0tot = s0b(1)
  do jb = 2, mbs
     s0tot = s0tot + s0b(jb)
  end do


! actual meteorological profiles in the model domain z(1)-z(n-1)
  call load1
  do k=1,n-1
     zx(k)=etw(k)
  enddo

! Even if etw(n-1) should be much lower than 11 km, check anyway
  if(zx(n-1) >= 11000.d0) then
     write(0,*) 'Error: etw(n-1) >= 11 km, check the vertical grid'
     stop 'Stopped by SR initr'
  end if
      ! dz: height increment up to the tropopause [m]
  dz=(11000.-zx(n-1))/7.

  gamma=0.0065
  do k=n,n+6
     zx(k)=zx(k-1)+dz
     tx(k)=tx(k-1)-gamma*dz
     px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
     rf=.3
     xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
     rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
     rnaer(k)=100.
  enddo

  k=n+7
  zx(k)=20000.
  tx(k)=tx(k-1)
  px(k)=px(k-1)*dexp(-g*(zx(k)-zx(k-1))/(r0*tx(k)))
  rf=.02
  xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
  rnaer(k)=0.

  k=k+1
  zx(k)=30000.
  gamma=-0.001
  tx(k)=tx(k-1)-gamma*10000.
  px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
  rf=.005
  xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
  rnaer(k)=0.

  k=k+1
  zx(k)=40000.
  gamma=-0.0026
  tx(k)=tx(k-1)-gamma*10000.
  px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
  rf=.00005
  xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
  rnaer(k)=0.

  k=k+1
  zx(k)=50000.
  gamma=-0.0018
  tx(k)=tx(k-1)-gamma*10000.
  px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
  rf=.000002
  xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
  rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
  rnaer(k)=0.

! fictitious level at infinity
  zx(nrlev)=zx(nrlay)+50000._dp
  tx(nrlev)=210.
  px(nrlev)=0._dp
  xm1x(nrlev)=0._dp
  rhox(nrlev)=0._dp

! Layer thicknesses calculated during initialisation
  do i=1,nrlay
     thkx(i)=zx(i+1)-zx(i)
  end do

! aerosol type: 1 rural 2 urban 3 maritme 4 tropospheric
  do k=1,nrlay
     ntypax(k)=2
  enddo

! ozone concentration profile
!----------------------------
! interpolate unreduced and reduced ozone amounts from craig
! table, save preliminarily the total amounts in arrays eta
! and xi which will be used for h2o and co2 later.
  do i=1,nrlay
     do j=1,51
        jj=j
        zj=(j-1.)*1000.
        jp=jj+1
        zjp=zj+1000.
        if (j.eq.51) zjp=0.
        if (zx(i).le.zjp) go to 2000
     enddo
     stop 'error in subroutine initr'
2000 continue
     dd=(zx(i)-zj)/(zjp-zj)
     etax(i)=o3un(jj)+(o3un(jp)-o3un(jj))*dd
  enddo
  etax(nrlev)=0.
! path lengths for every layer
  do i=1,nrlay
     j=nrlev-i
     jp=j+1
     usav(i)=(etax(j)-etax(jp))*0.01
  enddo

! start calculate qmo3
  do i=1,nrlay
     j=nrlev-i
     jp=j+1
     qmo3x(i)=usav(j)/(2.3808*(px(j)-px(jp)))
  end do
  qmo3x(nrlev)=0.

  do i=n,nrlay
!    ip=i+1         ! jjb potential problem here, ip is defined but not used
     na=ntypax(i)
     if (na.gt.0.and.rnaer(i).gt.0.) then
        rf=xm1x(i)*px(i)/(p21(tx(i))*(.62198+.37802*xm1x(i)))
        xnaer=rnaer(i)*1.d+06
        do nw=1,8
           naf=nw
           if (rf.le.feux(nw)) go to 2030
        enddo
2030    nafm=naf-1                                ! jjb potential problem here, nafm can be equal to 0 and lead to out of bounds index below
        dd=(rf-feux(nafm))/(feux(naf)-feux(nafm))
        emdd=1.-dd
        xdd=xnaer*dd
        xemdd=xnaer*emdd
        do k=1,mb
           beax(k,i)=xemdd*seanew(nafm,k,na)+xdd*seanew(naf,k,na)
           baax(k,i)=xemdd*saanew(nafm,k,na)+xdd*saanew(naf,k,na)
           gax(k,i)=emdd*ganew(nafm,k,na)+dd*ganew(naf,k,na)
        enddo
     else
        do k=1,mb
           beax(k,i)=0.
           baax(k,i)=0.
           gax(k,i)=0.
        end do
     endif
  enddo

! Initialise frac, rho2wx, and rew: no clouds
  do i=1,nrlay
     fracx(i)=0.
     rho2wx(i)=0.
     rewx(i) = 0.d0
  enddo

end subroutine initr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine load1
!
! Description:
!    loading actual variables and actual optical parameters
!

!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1       10/2016  Removed /nox/ tauno2(nrlay) which was set to 0    <Josue Bock>
!                      (see SR tau in nrad.f90 for explanations)

!           07/2016  Removal of labeled do loops.                   <Josue Bock>
!                    Header
!                    Use module for parameters
!                    All explicit declarations and implicit none
!
! 1.0       ?        Original code.                                 <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

  USE config, ONLY : &
! Imported Parameters:
       mic

  USE constants, ONLY : &
! Imported Parameters:
       pi,              &
       r0                 ! Specific gas constant of dry air, in J/(kg.K)

  USE global_params, ONLY : &
! Imported Parameters:
       jptaerrad,           &
       n,                   &
       nf,                  &
       nrlay,               &
       nrlev,               &
       nka,                 &
       nkt,                 &
       mb,                  &
       mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: ia, jt, k, ka, ka0, l1
  integer ::lu0, luf
  real (kind=dp) :: zeit, horang, rlat, rdec, u00, ru0, x0
  real (kind=dp) :: znum,zdenom,zfix   ! calculation of effective drop radius

! Common blocks:
  common /cb01/ tx(nrlev),px(nrlev),rhox(nrlev),xm1x(nrlev),      &
                rho2wx(nrlay),fracx(nrlay),zx(nrlev),thkx(nrlay), &
                beax(mb,nrlay),baax(mb,nrlay),gax(mb,nrlay),      &
                qmo3x(nrlev),rewx(nrlay), tsx
  real (kind=dp) :: tx,px,rhox,xm1x,rho2wx,fracx,zx,thkx, &
                    beax, baax, gax, qmo3x, rewx, tsx

  common /cb08/ re1(nkt), re2(nkt), re3(nkt)
  real (kind=dp) :: re1, re2, re3

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb18/ alat,declin                ! for the SZA calculation
  real (kind=dp) :: alat,declin

  common /cb40/ time,lday,lst,lmin,it,lcl,lct
  real (kind=dp) :: time
  integer :: lday, lst, lmin, it, lcl, lct

  common /cb41/ detw(n),deta(n),eta(n),etw(n)
  real (kind=dp) :: detw, deta, eta, etw

  common /cb49/ qabs(mb,nkt,nka,jptaerrad),qext(mb,nkt,nka,jptaerrad),asym(mb,nkt,nka,jptaerrad)
  real (kind=dp) :: qabs,qext,asym

  common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), &
                e(nkt),dew(nkt),rq(nkt,nka)
  real (kind=dp) :: enw,ew,rn,rw,en,e,dew,rq

  common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
  real (kind=dp) :: ff,fsum
  integer :: nar

  common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
  real (kind=dp) :: theta, thetl, t, talt, p, rho

  common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
  real (kind=dp) :: xm1, xm2, feu, dfddt, xm1a, xm2a

      !common /neu/ icld ! jjb 13/10/2016: added, see comment at the end of this SR
      !integer icld

!- End of header ---------------------------------------------------------------

      tsx=t(1)
      tx(1)=t(2)                               ! jjb CHECK this! For me, this should be t(1) (same as pressure below)
      px(1)=p(1)
      xm1x(1)=xm1(2)
      rhox(1)=px(1)/(r0*tx(1)*(1.+.608*xm1x(1)))
      rho2wx(1)=xm2(2)
      do k=2,n-1
         x0=0.5*detw(k)/deta(k)
         tx(k)=t(k)+(t(k+1)-t(k))*x0
         px(k)=p(k)+(p(k+1)-p(k))*x0
         xm1x(k)=xm1(k)+(xm1(k+1)-xm1(k))*x0
         rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
         rho2wx(k)=xm2(k+1)
      enddo

      if (mic) then
      
! calculate u0 from geogr. latitude, declination and hourangle
! make correction because of spherical surface of the earth
      zeit=lst*3600.+lmin*60.
      horang=7.272205e-05*zeit-pi
! pi/180=1.745329e-02
      rlat=alat*1.745329e-02
      rdec=declin*1.745329e-02
      u00=cos(rdec)*cos(rlat)*cos(horang)+sin(rdec)*sin(rlat)
      ru0=6371.*u00
      u0=8./(dsqrt(ru0**2+102000.)-ru0)

      ! Wavelenght bands: 1-6 = shortwave (sun), 7-18 = longwave (earth)
      ! 24h per day: earth longwave radiation taken into account
      lu0=7
      luf=mb
      ! Depending on solar zenith angle: sun shortwave radiation also accounted for
      if (u0.gt.1.d-02) lu0=1

      ! Initialisation
      do k=1,n-1
         do l1=lu0,luf
            baax(l1,k)=0.
            beax(l1,k)=0.
            gax(l1,k)=0.
         end do
      end do

      do k=1,n-1
         ka0=nar(k+1)
! the quantities needed in the radiation code are given by
! the sum over ia and jt of: qabs(ia,jt,l)*pi*r**2*ff(jt,ia,k)
! (same formula for qext) asym has weighting factor qsca*f(k)
     do ia=1,nka
        ka=ka0
        if (rn(ia).lt.0.5.and.ka0.eq.3) ka=2 ! ka=1 urban, ka=2 rural, ka=3 ocean
        do jt=1,nkt
           x0=pi*1.d-6*rq(jt,ia)**2*ff(jt,ia,k+1)
           do l1=lu0,luf
              baax(l1,k)=baax(l1,k)+qabs(l1,jt,ia,ka)*x0
              beax(l1,k)=beax(l1,k)+qext(l1,jt,ia,ka)*x0
              gax(l1,k)=gax(l1,k)+asym(l1,jt,ia,ka)*x0*(qext(l1,jt,ia,ka)-qabs(l1,jt,ia,ka))
           end do
        end do
     end do

     do l1=lu0,luf
        gax(l1,k)=gax(l1,k)/(beax(l1,k)-baax(l1,k)+.1e-15)
     end do

  end do ! k loop over vertical layers

! define frac, for use in SR frr,
!     and rew, for use in SR water
! start at index 2 (the case lcl = lct = 1 means no clouds)
!     and go up to index nf+1 (lct should be < nf anyway)
!     the remaining indexes nf+2:nrlay have been initialised to 0. in SR initr
      do k=2,nf+1
         if (k.ge.lcl .and. k.le.lct) then
            fracx(k-1) = 1.d0

            znum=0.d0
            zdenom=0.d0
            do jt=1,nkt
               zfix=0.d0
               do ia=1,nka
                  zfix = zfix + ff(jt,ia,k)
               end do
               znum   = znum + re3(jt)*zfix
               zdenom = zdenom + re2(jt)*zfix
            end do
            rewx(k-1) = znum/zdenom

         else
            fracx(k-1) = 0.d0
            rewx(k-1) = 0.d0
         end if
      end do

      end if

! < jjb 13/10/2016. Quick fix for consistency:
!
!     the variable icld was still needed in SR strahl (see file nrad.f), but it
!     was not currently defined in SR load1, only in the former SR load0.
!
!     According to a comment in file nrad.f, SR water should return zero optical
!     depths, unless SR load0 is used. This is achieved by defining icld = 0.
!     See comment after SR strahl, and the SR water in file nrad.f
!
!     icld is also used in SR frr (in file nrad.f), and has to be equal to 0 so that
!     frac is set to 1 between lcl and lct.
!
!     In order to keep the compatibility with SR load0, the variable icld has thus
!     been kept, and defined to zero in this subroutine.
!     However, it is likely that SR load0 will never be used again: if so, some
!     cleaning and little tuning will be required to completely remove it.

      !icld = 0

! jjb 13/10/2016 >

end subroutine load1
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine str
!
! Description:
!    radiative fluxes and heating rates at each minute

!

!
! History:
! Version   Date     Comment
! -------   ----     -------
  ! 16-Nov-2017  Josue Bock  The mic=false case is now directly handled in SR load1.
  !                            This avoids to duplicate the part of this SR which is
  !                            common to both cases.
! 1.2      10/2016   'clouds' removed after discussion with Andeas Bott   <Josue Bock>
!
!          07/2016   Common block /cb20/ was missing for 'clouds'         <Josue Bock>
!                    Comments / header
!                    Use module for parameters
!                    All explicit declarations and implicit none
!
! 1.1       ?        Change to add 'mic' case, and call load1             <Roland von Glasow>
!                       every 2 minutes
!
! 1.0       ?        Original code.                                       <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

  implicit none

!- End of header ---------------------------------------------------------------

!    if (lmin/2*2.eq.lmin) call load1
     call load1
     ! difference is nearly unplottable (cloudless case) but saves 20 % CPU time !

  call nstrahl

  call rotate_out

end subroutine str
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine rotate_in(linit)

! Description :
! -----------
!    Rotate arrays filled from main program, that are indexed from ground to top (bottom-up: _bu),
!    towards radiative code input arrays, indexed from top to ground (top-down: _td)


! Interface :
! ---------
!    SR rotate_in is called by SR initstr during initialisation, and by SR str during the run

! Input :
! -----
!    - /cb01/ contains all "_bu" arrays that need to be rotated to be used in the radiative code

! Output :
! ------
!    - /cb02/
!    - /cb16/    thk(nrlay)
!    - /height/  zx (nrlev)

! Externals :
! ---------
!    none


! Method :
! ------
!    All extra layer variables (T, P, ...) are constant.
!    Thus, the arrays are rotated over all layers/levels during initialisation only.
!    Else, only the meteorological layers will be rotated.



! Author :
! ------
!    Josué Bock


! Modifications :
! -------------
  !    02-Apr-2017  Josue Bock   First version of this routine
  !    13-Nov-2017  Josue Bock   removed ntypa and ntypd from /cb02/, unused
  !    14-Nov-2017  Josue Bock   further removed ntypd from /cb01/
!
! End modifications
!-----------------------------------------------------------------------------------------------------------------------



! Declarations:
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       n, &
       nrlay, &
       nrlev, &
       mb, mbs

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

  logical, intent(in) :: linit          ! called during initialisation, or not.

  integer :: jlbu, jltd                 ! running indexes (bu = bottom-up, td = top-down)
  integer :: nmin

  common /cb01/ tx(nrlev),px(nrlev),rhox(nrlev),xm1x(nrlev),         &
                rho2wx(nrlay),fracx(nrlay),zx_bu(nrlev),thkx(nrlay), &
                beax(mb,nrlay),baax(mb,nrlay),gax(mb,nrlay),         &
                qmo3x(nrlev),rewx(nrlay), tsx
  real (kind=dp) :: tx,px,rhox,xm1x,rho2wx,fracx,zx_bu,thkx, &
                    beax, baax, gax, qmo3x, rewx, tsx

  common /cb02/ t(nrlev),p(nrlev),rho(nrlev),xm1(nrlev),rho2w(nrlay), &
                frac(nrlay),ts
  real (kind=dp) :: t, p, rho, xm1, rho2w, frac, ts

  common /cb09/ rew(nrlay)
  real (kind=dp) :: rew

  common /cb16/ u0,albedo(mbs),thk(nrlay)
  real (kind=dp) :: u0, albedo, thk

  common /cb19/ berayl(6),bea(mb,nrlay),baa(mb,nrlay),ga(mb,nrlay)
  real (kind=dp) :: berayl, bea, baa, ga

  common /height/ zx(nrlev)        ! Level height [m] indexed from ground to top
  real (kind=dp) :: zx

  common /ozon/ qmo3(nrlev)
  real (kind=dp) :: qmo3


! Vertical grid arrays: rotate only during initialisation
  if(linit) then
     ! thk(nrlay)
     do jltd = 1,nrlay
        jlbu = nrlay - jltd + 1

        thk(jltd) = thkx(jlbu)
     end do
     ! zx(nrlev)
     do jltd = 1,nrlev
        jlbu = nrlev - jltd + 1

        zx(jltd) = zx_bu(jlbu)
     end do
  end if


! Define max index that will be rotated
  if(linit) then
     nmin = 1
  else
     nmin = nrlev - n + 1 ! fits nrlev size arrays. For nrlev size arrays, use nmin+1 instead
     !nmin = 1 ! fits nrlev size arrays. For nrlev size arrays, use nmin+1 instead
  end if

  do jltd = nmin, nrlev
     jlbu = nrlev - jltd + 1

     t(jltd) = tx(jlbu)
     p(jltd) = px(jlbu)
     rho(jltd) = rhox(jlbu)
     xm1(jltd) = xm1x(jlbu)
     qmo3(jltd) = qmo3x(jlbu)
     !(jltd) = x(jlbu)
  end do

  do jltd = MAX(nmin-1,1), nrlay
     jlbu = nrlay - jltd + 1

         rho2w(jltd) = rho2wx(jlbu)
         rew(jltd) = rewx(jlbu)
         frac(jltd) = fracx(jlbu)
         bea(:,jltd) = beax(:,jlbu)
         baa(:,jltd) = baax(:,jlbu)
         ga(:,jltd) = gax(:,jlbu)

  end do

  ts = tsx

end subroutine rotate_in
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine rotate_out

! Description :
! -----------
  !    Rotate arrays outputed by the radiative code, that are indexed from top to the ground
  !    (top-down: _td) to the corresponding variables used in the main program, that are indexed
  !    from the ground to the top (bottom-up: _bu).


! Interface :
! ---------
  !    SR rotate_out is called after the radiative code

! Input :
! -----
  !    - /cb10/ totrad(mb,nrlay)
  !    - /cb15/ fnseb,flgeg,hr(nrlay)

! Output :
! ------
  !    - /cb11/ totrad(mb,n)
  !    - /cb48/ sk,sl,dtrad(n)

! Externals :
! ---------
  !    none


! Method :
! ------
  !    Mistra layer #1 is surface
  !    Thus, fill dtrad and totrad_bu from indexes #2 to #n
  !    Index #2 in Mistra must match index #nrlay in the radiative code


! Author :
! ------
  !    Josue Bock


! Modifications :
! -------------
  !    13-Nov-2017  Josue Bock   First version of this routine

! == End of header =============================================================


! Declarations:
! ------------
! Modules used:

  USE global_params, ONLY : &
! Imported Parameters:
       n,                   &
       nrlay,               &
       mb

  USE precision, ONLY : &
! Imported Parameters:
       dp

  implicit none

! Local scalars:
  integer :: jzbu, jztd                 ! running indexes (bu = bottom-up, td = top-down)

! Common blocks:
  common /cb10/ totrad_td (mb,nrlay)
  real (kind=dp) :: totrad_td

  common /cb11/ totrad_bu (mb,n)
  real (kind=dp) :: totrad_bu

  common /cb15/ fnseb,flgeg,hr(nrlay)
  real (kind=dp) :: fnseb, flgeg, hr

  common /cb48/ sk,sl,dtrad(n),dtcon(n)
  real (kind=dp) :: sk, sl, dtrad, dtcon

! == End of declarations =======================================================


  sk = fnseb
  sl = flgeg

  dtrad(1)       = 0._dp
  totrad_bu(:,1) = 0._dp

  do jzbu = 2,n
     jztd = nrlay-jzbu+2

     dtrad(jzbu)=hr(jztd)
     totrad_bu(:,jzbu)=totrad_td(:,jztd)
  end do

end subroutine rotate_out
