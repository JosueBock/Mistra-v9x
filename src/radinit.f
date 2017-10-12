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
! radinit.f : radiation initialisation and start (SR str)
!
! contains the following subroutines:
!     - intrad
!     - initstr
!     - ipdata
!     - initr
!     - load0 (commented out, unused)
!     - load1
!     - str
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine intrad
!
! Description:
!    linear interpolation of radiation parameters qabs0, qext0, asym0
!    which are tabulated for droplet radii xw0 (microns)
!    and the percentage part xa0 of water on the refractive index m
!    with: m = ma + (mw -ma) *xa0,
!    mw: refractive index pure water, ma: refractive index pure aerosol
!    on the radiation parameters qabs, qext and asym
!    for droplet radii xw1 (in microns) and percentage parts of
!    pure water xa1 on the refractive index.
!
!    More information about the calculations of the tabulated values qabs0, qext0, asym0
!    using Mie theorie can be found in the following paper:
!    Bott, Sievers & Zdunkowski, J. Atmos. Sci., 47 (18), 2153-2166, 1990
!
!    Note that in this paper, alpha = 1-xa
! 
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      10/2016   Comments after personnal communication with A. Bott  <Josué Bock>
!
!          07/2016   Removal of labeled do loops
!                    Comments / header
!                    Use module for parameters
!                    All explicit declarations and implicit none
!                    Reindexing qabs/qext/asym for efficiency
!                    Partly rewritten, generic cases
!
! 1.1       ?        Improvements to avoid out of bounds indexes.         <Roland von Glasow>
!
! 1.0       ?        Original code.                                       <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE directories, ONLY :
! Imported Parameters:
     &     inpdir

      USE global_params, ONLY :
! Imported Parameters:
     &     nka,
     &     nkt,
     &     mb,
     &     mbs

      implicit none

! Local parameters:
      integer na0,nw0            ! Number of tabulated values for: percentage of water (na0) and total aerosol radius (nw0)
      parameter (na0=11,nw0=40)

! Local scalars:
      character (len=len_trim(inpdir)+11) fname
      double precision dx,dy,dxdy,dxdy1,dx1dy,dx1dy1 ! Interpolation coefficients
      double precision xa1,xw1                       ! Actual percentage of water (xa1, []) and total aerosol radius (xw1, [µm]) for each aerosol class
      integer ipc,i,j,k,ka,ia,jt

! Local arrays:
      double precision qabs0(mb,nw0,na0) ! tabulated values of absorption coefficient
      double precision qext0(mb,nw0,na0) ! tabulated values of extinction coefficient
      double precision asym0(mb,nw0,na0) ! tabulated values of asymmetry factor

      double precision xa0(na0) ! percentage part of pure water, in volume (tabulated values) [%]
      data xa0 /0.0,0.2,0.4,0.6,0.7,0.8,0.85,0.9,0.95,0.975,1.0/

      double precision xw0(nw0) ! total radius of the scattering spheres (tabulated values) [µm]
      data xw0 /0.01,0.0125,0.015,0.02,0.025,0.03,0.04,0.05,0.06,0.08
     & ,0.1,0.125,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.8
     & ,1.0,1.25,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0
     & ,10.0,12.5,15.0,20.0,25.0,30.0,40.0,50.0,60.0,80.0/

! Common blocks:
      common /cb49/ qabs(mb,nkt,nka,3),qext(mb,nkt,nka,3),
     &              asym(mb,nkt,nka,3)
      double precision qabs,qext,asym

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), ! only rn (dry radius), rq (total radius) are used here
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

!- End of header ---------------------------------------------------------------


! Open input data files
! ---------------------

! input qabs, qext and asym for radiation parameters of particles
!   odd unit numbers: short wave data
!   even unit numbers: long wave data
      fname=TRIM(inpdir)//"urbankw.dat"
      open (unit=51, file=fname, status='old')
      fname=TRIM(inpdir)//"urbanlw.dat"
      open (unit=52, file=fname, status='old')
      fname=TRIM(inpdir)//"ruralkw.dat"
      open (unit=53, file=fname, status='old')
      fname=TRIM(inpdir)//"rurallw.dat"
      open (unit=54, file=fname, status='old')
      fname=TRIM(inpdir)//"ozeankw.dat"
      open (unit=55, file=fname, status='old')
      fname=TRIM(inpdir)//"ozeanlw.dat"
      open (unit=56, file=fname, status='old')
! input format to read these files
 5000 format (3e16.8)
      ipc=51

      do ka=1,3  ! ka=1 urban, ka=2 rural, ka=3 ocean

! Read input data file: tabulated values
! --------------------------------------
      do j=1,na0
         do i=1,nw0
            do k=1,mbs
               read(ipc,5000) qabs0(k,i,j),qext0(k,i,j),asym0(k,i,j)
            enddo
         enddo
      enddo
      close (ipc)
      ipc=ipc+1
      do j=1,na0
         do i=1,nw0
            do k=mbs+1,mb
               read(ipc,5000) qabs0(k,i,j),qext0(k,i,j),asym0(k,i,j)
            enddo
         enddo
      enddo
      close (ipc)
      ipc=ipc+1

! Loop over the aerosol class
! ---------------------------
      do ia=1,nka
         do jt=1,nkt
            ! Values of total radius (xw1) and volume mixing ratio of water (xa1) for each aerosol class
            ! ------------------------------------------------------------------------------------------
            xw1=rq(jt,ia)
            xa1=1.-(rn(ia)/rq(jt,ia))**3
            xa1=max(xa1,1.d-2)

            ! Locate xw1 in the tabulated xw0 grid, and calculate interpolation factor dx
            ! ---------------------------------------------------------------------------
            do i=1,nw0
               if (xw1.le.xw0(i)) then
                  if (i.eq.1) then ! if rq <= xw0(1)
                     dx=1.d0
                     print*,'Warning: in SR intrad, rq <= xw0(1)',ia,jt,
     &                      xw1
                  else ! general case: xw0(i-1) < rq <= xw0(i)
                     dx=(xw1-xw0(i-1))/(xw0(i)-xw0(i-1))
                  end if
                  go to 2000 ! exit do loop
               end if
            end do
            ! Last case: do loop ended without going to 2000, means xw0(nw0) < rq
            print*,'Warning: in SR intrad, rq > xw0(nw0)',ia,jt,xw1
            dx=0
            i=nw0+1

 2000       continue

            ! Locate xa1 in the tabulated xa0 grid, and calculate interpolation factor dy
            ! ---------------------------------------------------------------------------
            do j=1,na0
               if (xa1.le.xa0(j)) then
                  if (j.eq.1) then ! if xa1 <= xa0(1), this case shouldn't happen
                     dy=1.d0
                     print*,'Error in SR intrad: xa1 < xa0(1)',ia,jt,xa1
                     stop
                  else ! general case xa0(i-1) < xa1 <= xa0(i)
                     dy=(xa1-xa0(j-1))/(xa0(j)-xa0(j-1))
                  end if
                  go to 2010 ! exit do loop
               end if
            end do
            ! Last case: do loop ended without going to 2010, means xa0(na0) < xa1
            print*,'Error in SR intrad: xa1 > xa0(na0)',ia,jt,xa1
            stop

 2010       continue

            ! Interpolation factors
            ! ---------------------
            dxdy=dx*dy
            dxdy1=dx*(1.d0-dy)
            dx1dy=dy*(1.d0-dx)
            dx1dy1=(1.d0-dy)*(1.d0-dx)

            ! Interpolated values of qabs, qext and asym
            ! ------------------------------------------
            do k=1,mb
               qabs(k,jt,ia,ka)=dxdy*qabs0(k,MIN(i,nw0),MIN(j,na0))
     &              +dxdy1*qabs0(k,MIN(i,nw0),MAX(1,j-1))
     &              +dx1dy*qabs0(k,MAX(1,i-1),MIN(j,na0))
     &              +dx1dy1*qabs0(k,MAX(1,i-1),MAX(1,j-1))
               qext(k,jt,ia,ka)=dxdy*qext0(k,MIN(i,nw0),MIN(j,na0))
     &              +dxdy1*qext0(k,MIN(i,nw0),MAX(1,j-1))
     &              +dx1dy*qext0(k,MAX(1,i-1),MIN(j,na0))
     &              +dx1dy1*qext0(k,MAX(1,i-1),MAX(1,j-1))
               asym(k,jt,ia,ka)=dxdy*asym0(k,MIN(i,nw0),MIN(j,na0))
     &              +dxdy1*asym0(k,MIN(i,nw0),MAX(1,j-1))
     &              +dx1dy*asym0(k,MAX(1,i-1),MIN(j,na0))
     &              +dx1dy1*asym0(k,MAX(1,i-1),MAX(1,j-1))
            enddo

            end do ! jt loop
         end do ! ia loop
      end do ! ka loop

      end subroutine intrad
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine initstr
!
! Description:
!    first calculation of radiative fluxes and heating rates
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
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
      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     n1

      implicit none

! Local scalars:
      integer istr, k, lmin0, lst0
! Common blocks:
      common /cb15/ fnseb,flgeg,hr(n1)
      double precision fnseb, flgeg, hr

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb55/ dtrad0(n),dtrad1(n),sk0,sl0,sk1,sl1,time0,time2
      double precision dtrad0, dtrad1, sk0, sl0, sk1, sl1, time0, time2
!- End of header ---------------------------------------------------------------

! initialisation of radiation code
! --------------------------------
      call ipdata
      call initr
      call nstrahl
      call profr

      sk0=fnseb
      sl0=flgeg
      sk=sk0
      sl=sl0
      dtrad(1)=0.d0
      do k=2,n
         dtrad(k)=hr(k-1)
         dtrad0(k)=hr(k-1)
      end do
      istr=30

! addition of istr to the actual time
! -----------------------------------
      lmin0=lmin
      lst0=lst
      lmin=lmin+istr
      if (lmin.ge.60) then
         lmin=lmin-60
         lst=lst+1
         if (lst.eq.24) lst=0
      end if
      time2=float(istr)*60.d0
      time0=time2

      call load1
      call nstrahl
      call profr

      sk1=fnseb
      sl1=flgeg
      do k=2,n
         dtrad1(k)=hr(k-1)
      end do

! recorrection of time
! --------------------
      lmin=lmin0
      lst=lst0

      end subroutine initstr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine ipdata
!
! Description:
!    input data for radiation code:
!    reads from file 'pifm2.dat'
!    writes in several common blocks
!    descriptions of variables can be found at the end of file 'pifm2*.dat'
!    (part of these description written in German)
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1       10/2016  More cleaning, after checking the variables actually used   <Josue Bock>
!                      (several old variables from PIFM1, apparently)
!
!           07/2016  Header, comments                                            <Josue Bock>
!                    Use modules for parameters
!
! 1.0       ?        Original code.                                              <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE directories, ONLY :
! Imported Parameters:
     &     inpdir

      USE global_params, ONLY :
! Imported Parameters:
     &     mb,
     &     ncw

      implicit none

! Common blocks:
      common /aeag/ seanew(8,mb,4),saanew(8,mb,4),ganew(8,mb,4),ff2(8)
      double precision seanew, saanew, ganew, ff2
      common /band1/ hk1(10),fk1o3(10)
      double precision hk1, fk1o3
      common /band2/ hk2(8),c2h2o(3,11,8)
      double precision hk2, c2h2o
      common /band3/ hk3(12),c3h2o(3,11,12)
      double precision hk3, c3h2o
      common /band4/ hk4(7),c4h2o(3,11,7)
      double precision hk4, c4h2o
      common /band5/ hk5(12),c5h2o(3,11,12)
      double precision hk5, c5h2o
      common /band6/ hk6(5),c6h2o(3,11,5)
      double precision hk6, c6h2o
      common /band7/ hk7(2),c7h2o(3,19,2)
      double precision hk7, c7h2o
      common /band8/ hk8(3),c8h2o(3,19,3)
      double precision hk8, c8h2o
      common /band9/ hk9(4),c9h2o(3,19,4)
      double precision hk9, c9h2o
      common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
      double precision hk10, c10h2o, c10ch4, c10n2o
      common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
      double precision hk11, c11h2o, c11ch4, c11n2o
      common /band12/ hk12(5),c12o3(3,19,5),c12h2o(3,19)
      double precision hk12, c12o3, c12h2o
      common /band13/ hk13(2),c13h2o(3,19,2)
      double precision hk13, c13h2o
      common /band14/ hk14(10),c14hca(3,19,10),c14hcb(3,19,10)
      double precision hk14, c14hca, c14hcb
      common /band15/ hk15(12),c15hca(3,19,12),c15hcb(3,19,12)
      double precision hk15, c15hca, c15hcb
      common /band16/ hk16(7),c16h2o(3,19,7)
      double precision hk16, c16h2o
      common /band17/ hk17(7),c17h2o(3,19,7)
      double precision hk17, c17h2o
      common /band18/ hk18(8),c18h2o(3,19,8)
      double precision hk18, c18h2o

!      common /plancd/ ttab(35),pibtab(12,35) ! jjb was used in FN fst4 (called by SR plancktab), both unused now
!      double precision ttab, pibtab          !     left for backward compatibility
      double precision ttab(35),pibtab(12,35)

      common /was1/ ret(ncw),r2wt(ncw),b2wt(ncw,mb),w2wt(ncw,mb),
     &              g2wt(ncw,mb)
      double precision ret, r2wt, b2wt, w2wt, g2wt

! Local scalars:
      character (len=len_trim(inpdir)+16) fname

!- End of header ---------------------------------------------------------------

      fname=TRIM(inpdir)//'pifm2_161019.dat'
      open(unit=57, file=fname, status='old')

! input format to read this file
 10   format(5e16.8)

      read(57,10) ttab,pibtab               ! could be deleted if fst4 and plancktab are deleted
      read(57,10) ret,r2wt,b2wt,w2wt,g2wt
      read(57,10) seanew,saanew,ganew
      read(57,10) hk1,fk1o3
      read(57,10) hk2,c2h2o
      read(57,10) hk3,c3h2o
      read(57,10) hk4,c4h2o
      read(57,10) hk5,c5h2o
      read(57,10) hk6,c6h2o
      read(57,10) hk7,c7h2o
      read(57,10) hk8,c8h2o
      read(57,10) hk9,c9h2o
      read(57,10) hk10,c10h2o,c10ch4,c10n2o
      read(57,10) hk11,c11h2o,c11ch4,c11n2o
      read(57,10) hk12,c12o3,c12h2o
      read(57,10) hk13,c13h2o
      read(57,10) hk14,c14hca,c14hcb
      read(57,10) hk15,c15hca,c15hcb
      read(57,10) hk16,c16h2o
      read(57,10) hk17,c17h2o
      read(57,10) hk18,c18h2o

      close(57)

      end subroutine ipdata
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine initr
!
! Description:
!    standard atmosphere between etw(n) and 50000 meters.
!    input of constant radiation parameters
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2       10/2016  corrected several array size
!                    initialisation of qmo3(np) was missing
!                    initialisation of rnaer(nrfl=n1) was missing
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

      USE directories, ONLY :
! Imported Parameters:
     &     inpdir

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     n1,
     &     n4,
     &     nka,
     &     mb,
     &     mbs,
     &     mbir,
     &     nrfl,
     &     np

      implicit none

! Local scalars:
      character  (len=len_trim(inpdir)+12) fname
      integer i, j, jj, jp, k, na, naf, nafm, nw
      double precision gamma, rf, zj, zjp
      double precision dd, emdd, xdd, xemdd, xnaer

! Local arrays:
      double precision etax(n4), usav(n1)
      double precision rnaer(n1)          ! total number concentration of aerosol particles [cm**-3]

! Internal function:
      double precision p21, tt

! Common blocks:
      common /aeag/ seanew(8,mb,4),saanew(8,mb,4),ganew(8,mb,4),ff2(8)
      double precision seanew, saanew, ganew, ff2

      common /cb02/ tx(n4),px(n4),rhox(n4),xm1x(n4),rho2x(n1),frac(n1),
     & ts,ntypa(n1),ntypd(n1)
      double precision tx,px,rhox,xm1x,rho2x,frac,ts
      integer ntypa,ntypd

      common /cb16/ u0,albedo(mbs),thk(n1)
      double precision u0, albedo, thk

      common /cb19/ berayl(6),bea(mb,n1),baa(mb,n1),ga(mb,n1)
      double precision berayl, bea, baa, ga

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ r0,r1,g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision r0,r1,g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw

!      common /cb56/ o3un(52),o3r(52),vis(n2),sigea(8,6,4),    ! /cb56/ was only used with SR load0
!     &              sigaa(8,6,4),gaa(8,6,4),sean(8,4),feux(8)
!      double precision o3un,o3r,vis,sigea,sigaa,gaa,sean,feux
      double precision o3un(52),o3r(52),sigea(8,6,4),          ! keep this for array size declaration
     &               sigaa(8,6,4),gaa(8,6,4),sean(8,4),feux(8) ! and potential backward compatibility

      common /height/ zx(n4)
      double precision zx

      common /ozon/ qmo3(n4)
      double precision qmo3

      common /tmp2/ as(mbs),ee(mbir)
      double precision as, ee

!- End of header ---------------------------------------------------------------

      p21(tt)=610.7*dexp(17.15*(tt-273.15)/(tt-38.33))


! albedo of the ground for the six solar wavelength regions of the radiation code
      albedo(:)=0.05
!      albedo(:)=0.8       ! snow
      as = albedo
! emissivity of the ground
      ee(:)=1.0


! read optical parameters
! -----------------------
!     input constants for radiation code
      fname=TRIM(ADJUSTL(inpdir))//'initr_v3.dat'
      open (unit=58, file=fname, status='old')

! input formats to read this file
 6000 format (5e15.7)

      read (58,6000) (((sigea(i,j,k),i=1,8),j=1,6),k=1,4)
      read (58,6000) (((sigaa(i,j,k),i=1,8),j=1,6),k=1,4)
      read (58,6000) (((gaa(i,j,k),i=1,8),j=1,6),k=1,4)
      read (58,6000) ((sean(i,k),i=1,8),k=1,4)
      read (58,6000) (berayl(j),j=1,6)
      read (58,6000) (o3un(i),i=1,52) ! unreduced ozone amount from Craig table
      read (58,6000) (o3r(i),i=1,52)  ! reduced ozone amount from Craig table
      read (58,6000) (feux(i),i=1,8)

      close (58)

! actual meteorological profiles in the model domain z(1)-z(n)
      call load1
      do k=1,n-1
         zx(k)=etw(k)
      enddo
      gamma=0.0065
      do k=n,n+6
         zx(k)=zx(k-1)+1000.
         tx(k)=tx(k-1)-gamma*1000.
         px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
         rf=.3
         xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
         rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
         rnaer(k)=100.
      enddo

      k=n+7
      zx(k)=20000.
      tx(k)=tx(k-1)
      px(k)=px(k-1)*dexp(-g*10000./(r0*tx(k)))
      rf=.1
      xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
      rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
      rnaer(k)=0.

      k=k+1
      zx(k)=30000.
      gamma=-0.001
      tx(k)=tx(k-1)-gamma*10000.
      px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
      rf=.8
      xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
      rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
      rnaer(k)=0.

      k=k+1
      zx(k)=40000.
      gamma=-0.0027
      tx(k)=tx(k-1)-gamma*10000.
      px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
      rf=.05
      xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
      rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
      rnaer(k)=0.

      k=k+1
      zx(k)=50000.
      gamma=-0.0015
      tx(k)=tx(k-1)-gamma*10000.
      px(k)=px(k-1)*(tx(k)/tx(k-1))**(g/(r0*gamma))
      rf=.01
      xm1x(k)=0.62198*rf/(px(k)/p21(tx(k))-0.37802*rf)
      rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
      rnaer(k)=0.

! fictitious level at infinity
      zx(np)=zx(nrfl)+50000.
      tx(np)=tx(nrfl)
      px(np)=0.
      xm1x(np)=0.
      rhox(np)=0.

! Layer thicknesses calculated during initialisation
      do i=1,nrfl
         thk(i)=zx(i+1)-zx(i)
      end do

! aerosol type: 1 rural 2 urban 3 maritme 4 tropospheric
! droplet type: 1-3 cumulus 4 best
      do k=1,nrfl
         ntypa(k)=2
         ntypd(k)=4
      enddo

! interpolate unreduced and reduced ozone amounts from craig
! table, save preliminarily the total amounts in arrays eta
! and xi which will be used for h2o and co2 later.
      do i=1,nrfl
         do j=1,51
            jj=j
            zj=(j-1.)*1000.
            jp=jj+1
            zjp=zj+1000.
            if (j.eq.51) zjp=0.
            if (zx(i).le.zjp) go to 2000
         enddo
         stop 'error in subroutine initr'
 2000    continue
        dd=(zx(i)-zj)/(zjp-zj)
        etax(i)=o3un(jj)+(o3un(jp)-o3un(jj))*dd
      enddo
      etax(np)=0.
! path lengths for every layer
      do i=1,nrfl
         j=np-i
         jp=j+1
         usav(i)=(etax(j)-etax(jp))*0.01
      enddo

! start calculate qmo3
      do i=1,nrfl
         j=np-i
         jp=j+1
         qmo3(i)=usav(j)/(2.3808*(px(j)-px(jp)))
      end do
      qmo3(np)=0.

      do i=n,nrfl
!        ip=i+1         ! jjb potential problem here, ip is defined but not used
         na=ntypa(i)
         if (na.gt.0.and.rnaer(i).gt.0.) then
            rf=xm1x(i)*px(i)/(p21(tx(i))*(.62198+.37802*xm1x(i)))
            xnaer=rnaer(i)*1.d+06
            do nw=1,8
               naf=nw
               if (rf.le.feux(nw)) go to 2030
            enddo
 2030       nafm=naf-1                                ! jjb potential problem here, nafm can be equal to 0 and lead to out of bounds index below
            dd=(rf-feux(nafm))/(feux(naf)-feux(nafm))
            emdd=1.-dd
            xdd=xnaer*dd
            xemdd=xnaer*emdd
            do k=1,mb
               bea(k,i)=xemdd*seanew(nafm,k,na)+xdd*seanew(naf,k,na)
               baa(k,i)=xemdd*saanew(nafm,k,na)+xdd*saanew(naf,k,na)
               ga(k,i)=emdd*ganew(nafm,k,na)+dd*ganew(naf,k,na)
            enddo
         else
            do k=1,mb
               bea(k,i)=0.
               baa(k,i)=0.
               ga(k,i)=0.
            end do
         endif
      enddo

! no clouds
      do i=1,nrfl
         frac(i)=0.
         rho2x(i)=0.
      enddo

      end subroutine initr
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! jjb 10/03/2016: this SR is not called, thus commented for clarity
!     28/07/2016: updated nevertheless, pi taken from module 'constants'
!     05/10/2016: removed unused CB49
!     13/10/2016: there was probably a problem with variable icld, which is not
!                 initialised (only defined to 1 if rho2x(i).ge.1.0e-5)
!
c$$$      subroutine load0
c$$$
c$$$      USE constants, ONLY :
c$$$! Imported Parameters:
c$$$     & pi
c$$$
c$$$      USE global_params, ONLY :
c$$$! Imported Parameters:
c$$$     &     mb,
c$$$     &     mbs,
c$$$     &     n1,
c$$$     &     n2,
c$$$     &     nrfl,
c$$$     &     np,
c$$$     &     n,nf,nka,nkt ! jjb remove unused
c$$$
c$$$      implicit double precision (a-h,o-z)
c$$$! loading actual variables and constant optical parameters
c$$$
c$$$      double precision rnaer(n1) ! jjb removed from cb02, now local variable.
c$$$
c$$$      common /aeag/ seanew(8,mb,4),saanew(8,mb,4),ganew(8,mb,4),ff2(8)
c$$$      common /cb02/ tx(n4),px(n4),rhox(n4),xm1x(n4),rho2x(n1),frac(n1),
c$$$     & ts,ntypa(n1),ntypd(n1)
c$$$      double precision tx,px,rhox,xm1x,rho2x,fracx,ts
c$$$      integer ntypa,ntypd
c$$$
c$$$      common /cb16/ u0,albedo(mbs),thk(n1)
c$$$      double precision u0, albedo, thk
c$$$
c$$$      common /cb18/ alat,declin                ! for the SZA calculation
c$$$      double precision alat,declin
c$$$
c$$$      common /cb19/ berayl(6),bea(mb,n1),baa(mb,n1),ga(mb,n1)
c$$$      double precision berayl, bea, baa, ga
c$$$
c$$$      common /cb40/ time,lday,lst,lmin,it,lcl,lct
c$$$      common /cb41/ detw(n),deta(n),eta(n),etw(n)
c$$$      common /cb44/ r0,r1,g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
c$$$     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw
c$$$      common /cb52/ f(nkt,nka,n),fsum(n),nar(n)
c$$$      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
c$$$      double precision theta, thetl, t, talt, p, rho
c$$$      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
c$$$      common /cb56/ o3un(52),o3r(52),vis(n2),sigea(8,6,4),
c$$$     &              sigaa(8,6,4),gaa(8,6,4),sean(8,4),feux(8)
c$$$      common /neu/ icld
c$$$      integer icld
c$$$
c$$$      common /nox/ tauno2(n1)
c$$$      double precision tauno2
c$$$
c$$$! fsum0 estimated aerosol concentration with radius lower than
c$$$! minimum aerosol radius of current model aerosol distribution
c$$$      data fsum0 /5000./
c$$$      p21(tt)=610.7*dexp(17.15*(tt-273.15)/(tt-38.33))
c$$$      ts=t(1)
c$$$      tx(1)=t(2)
c$$$      px(1)=p(1)
c$$$      xm1x(1)=xm1(2)
c$$$      rhox(1)=px(1)/(r0*tx(1)*(1.+.608*xm1x(1)))
c$$$      rho2x(1)=xm2(2)*p(2)/(r0*t(2))
c$$$      rnaer(1)=fsum(2)+fsum0
c$$$      do 1000 k=2,n-1
c$$$         x0=0.5*detw(k+1)/deta(k)
c$$$         tx(k)=(1.-x0)*t(k+1)+x0*t(k)
c$$$         px(k)=(1.-x0)*p(k+1)+x0*p(k)
c$$$         xm1x(k)=(1.-x0)*xm1(k+1)+x0*xm1(k)
c$$$         rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
c$$$         rho2x(k)=xm2(k+1)*p(k+1)/(r0*t(k+1))
c$$$         rnaer(k)=fsum(k+1)+fsum0
c$$$ 1000    continue
c$$$! calculate u0 from geogr. latitude, declination and hourangle
c$$$! make correction because of spherical surface of the earth
c$$$      zeit=lst*3600.+lmin*60.
c$$$      horang=7.272205e-05*zeit-pi
c$$$! pi/180=1.745329e-02
c$$$      rlat=alat*1.745329e-02
c$$$      rdec=declin*1.745329e-02
c$$$      u00=cos(rdec)*cos(rlat)*cos(horang)+sin(rdec)*sin(rlat)
c$$$      ru0=6371.*u00
c$$$      u0=8./(dsqrt(ru0**2+102000.)-ru0)
c$$$! aerosol type: 1 rural 2 urban 3 maritme 4 tropospheric
c$$$! droplet type; best type, if previously defined as 4
c$$$      do 1010 i=1,nrfl
c$$$         ntypa(i)=2
c$$$         ntypd(i)=4
c$$$         if (ntypd(i).ge.4) goto 2000
c$$$         if (rho2x(i).gt.0.25e-3) ntypd(i)=2
c$$$         if (rho2x(i).gt.0.6e-3) ntypd(i)=3
c$$$ 2000    continue
c$$$         ip=i+1
c$$$         na=ntypa(i)
c$$$         if (na.gt.0.and.rnaer(i).gt.0.) goto 2010
c$$$         do k=1,mb
c$$$            bea(k,i)=0.
c$$$            baa(k,i)=0.
c$$$            ga(k,i)=0.
c$$$         enddo
c$$$         vis(i)=1.e20
c$$$         goto 1010
c$$$ 2010    rf=xm1x(i)*px(i)/(p21(tx(i))*(.62198+.37802*xm1x(i)))
c$$$         xnaer=rnaer(i)*1.d+06
c$$$         do 1030 nw=1,8
c$$$            naf=nw
c$$$ 1030       if (rf.le.feux(nw)) go to 2020
c$$$ 2020       nafm=naf-1
c$$$            dd=(rf-feux(nafm))/(feux(naf)-feux(nafm))
c$$$            emdd=1.-dd
c$$$            xdd=xnaer*dd
c$$$            xemdd=xnaer*emdd
c$$$            do k=1,mb
c$$$               bea(k,i)=xemdd*seanew(nafm,k,na)+xdd*seanew(naf,k,na)
c$$$               baa(k,i)=xemdd*saanew(nafm,k,na)+xdd*saanew(naf,k,na)
c$$$               ga(k,i)=emdd*ganew(nafm,k,na)+dd*ganew(naf,k,na)
c$$$            enddo
c$$$            vis(i)=3.912e-3/(xemdd*sean(nafm,na)+xdd*sean(naf,na)+
c$$$     &           .122666e-4)
c$$$ 1010 continue
c$$$c------------------------------------------------------------------
c$$$! optical depth of qno2 only used in first wavelength bin optical
c$$$c
c$$$      qno=0.0
c$$$      do k=1,nrfl
c$$$!        dzd8=thk(j)*0.125 ! jjb has to be defined AFTER j is defined below
c$$$         j=np-k
c$$$         dzd8=thk(j)*0.125 ! jjb has to be defined AFTER j
c$$$         jp=j+1
c$$$         rp=rho(j)*p(j)
c$$$!         qnu=7.8407708e-12*rp*qno2(j)
c$$$!         qnz=7.8407708e-12*rp*(qno2(jp)+qno2(j))*0.5 ! jjb 21/10/2016 qno2 outdated (PIFM1)
c$$$!         tauno2(k)=dzd8*(qnu+qno+6.*qnz)
c$$$!         qno=qnu
c$$$      enddo
c$$$c
c$$$! fractional cloudiness
c$$$!!      clouds=0.                 ! jjb 12/10/2016 'clouds' removed from /cb20/, was not used
c$$$      do 1050 i=1,nrfl
c$$$         j=nrfl+1-i
c$$$         if(rho2x(i).ge.1.0e-5) then
c$$$            frac(j)=1.0
c$$$            icld=1
c$$$         else
c$$$            frac(j)=0.0
c$$$         endif
c$$$!!         clouds=clouds+frac(j)  ! jjb 12/10/2016 'clouds' removed from /cb20/, was not used
c$$$ 1050 continue
c$$$c------------------------------------------------------------------
c$$$      end subroutine load0
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine load1
!
! Description:
!    loading actual variables and actual optical parameters
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1       10/2016  Removed /nox/ tauno2(n1) which was set to 0    <Josue Bock>
!                      (see SR tau in nard.f for explanations)

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

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     & n,
     & n1,
     & n4,
     & nka,
     & nkt,
     & mb,
     & mbs

      implicit none

! Local scalars:
      integer ia, jt, k, ka, ka0, l1
      integer lu0, luf
      double precision zeit, horang, rlat, rdec, u00, ru0, x0

! Common blocks:
      common /cb02/ tx(n4),px(n4),rhox(n4),xm1x(n4),rho2x(n1),frac(n1),
     & ts,ntypa(n1),ntypd(n1)
      double precision tx,px,rhox,xm1x,rho2x,frac,ts
      integer ntypa,ntypd

      common /cb16/ u0,albedo(mbs),thk(n1)
      double precision u0, albedo, thk

      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb19/ berayl(6),bea(mb,n1),baa(mb,n1),ga(mb,n1)
      double precision berayl, bea, baa, ga

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ r0,r1,g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision r0,r1,g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb49/ qabs(mb,nkt,nka,3),qext(mb,nkt,nka,3),
     &              asym(mb,nkt,nka,3)
      double precision qabs,qext,asym

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      double precision ff,fsum
      integer nar

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      double precision xm1, xm2, feu, dfddt, xm1a, xm2a

      common /neu/ icld ! jjb 13/10/2016: added, see comment at the end of this SR
      integer icld

!- End of header ---------------------------------------------------------------

      ts=t(1)
      tx(1)=t(2)                               ! jjb CHECK this! For me, this should be t(1) (same as pressure below)
      px(1)=p(1)
      xm1x(1)=xm1(2)
      rhox(1)=px(1)/(r0*tx(1)*(1.+.608*xm1x(1)))
      do k=2,n-1
         x0=0.5*detw(k)/deta(k)
         tx(k)=t(k)+(t(k+1)-t(k))*x0
         px(k)=p(k)+(p(k+1)-p(k))*x0
         xm1x(k)=xm1(k)+(xm1(k+1)-xm1(k))*x0
         rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
      enddo

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
            baa(l1,k)=0.
            bea(l1,k)=0.
            ga(l1,k)=0.
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
                  baa(l1,k)=baa(l1,k)+qabs(l1,jt,ia,ka)*x0
                  bea(l1,k)=bea(l1,k)+qext(l1,jt,ia,ka)*x0
                  ga(l1,k)=ga(l1,k)+asym(l1,jt,ia,ka)*x0
     &                 *(qext(l1,jt,ia,ka)-qabs(l1,jt,ia,ka))
               end do
            end do
         end do

         do l1=lu0,luf
            ga(l1,k)=ga(l1,k)/(bea(l1,k)-baa(l1,k)+.1e-15)
         end do

      end do ! k loop over vertical layers


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

      icld = 0

! jjb 13/10/2016 >

      end subroutine load1
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine str (mic)
!
! Description:
!    radiative fluxes and heating rates at each minute

!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
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

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     n1,
     &     n4,
     &     nka,
     &     mbs

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      logical mic

! Local scalars:
      integer k
      double precision zeit, horang, rlat, rdec, u00, ru0, x0

! Common blocks:
      common /cb02/ tx(n4),px(n4),rhox(n4),xm1x(n4),rho2x(n1),frac(n1),
     & ts,ntypa(n1),ntypd(n1)
      double precision tx,px,rhox,xm1x,rho2x,frac,ts
      integer ntypa,ntypd

      common /cb15/ fnseb,flgeg,hr(n1)
      double precision fnseb, flgeg, hr

      common /cb16/ u0,albedo(mbs),thk(n1)
      double precision u0, albedo, thk

      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ r0,r1,g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision r0,r1,g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhow,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      double precision xm1, xm2, feu, dfddt, xm1a, xm2a
!- End of header ---------------------------------------------------------------

      if (mic) then
!         if (lmin/2*2.eq.lmin) call load1
         call load1
         ! difference is nearly unplottable (cloudless case) but saves 20 % CPU time !
      else
! next lines are copied from SR load1
         ts=t(1)
         tx(1)=t(2)                      ! jjb CHECK this! For me, this should be t(1) (same as pressure below)
         px(1)=p(1)
         xm1x(1)=xm1(2)                  ! jjb CHECK this! For me, this should be xm1(1) (same as pressure above)
         rhox(1)=px(1)/(r0*tx(1)*(1.+.608*xm1x(1)))
         do  k=2,n-1
            x0=0.5*detw(k)/deta(k)
            tx(k)=t(k)+(t(k+1)-t(k))*x0
            px(k)=p(k)+(p(k+1)-p(k))*x0
            xm1x(k)=xm1(k)+(xm1(k+1)-xm1(k))*x0
            rhox(k)=px(k)/(r0*tx(k)*(1.+.608*xm1x(k)))
         enddo
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
      endif

      call nstrahl
      sk=fnseb
      sl=flgeg
      dtrad(1)=0.
      do k=2,n
         dtrad(k)=hr(k-1)
      end do

      end subroutine str
