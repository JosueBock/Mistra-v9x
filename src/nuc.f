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
!                  *************************************
!                  * nucleation subroutines for MISTRA *
!                  *************************************
!
! See detailed explanations in the header of SR appnucl
! See MISTRA manual for full references and details.
!
! This file contains the following subroutines and functions:
!     - appnucl2
!     - appnucl
!     - dmean
!     - ternucl
!     - oionucl
!     - J_nuc   (function)
!     - nucout1
!     - nucout2
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      module mod_nuc

! Purpose:
! --------
!     set the number of nucleating vapours, and declare variables used in the nucleation scheme

! History:
! --------
!     30-11-2016   J. Bock   first version of this module

      implicit none
      save

      integer, parameter :: nvap = 1 ! Number of condensible vapors, please adjust!
      integer :: ical(nvap)
      integer :: ivap(nvap)
      integer :: ind_H2SO4, ind_NH3, ind_OIO
      double precision :: m_vap(nvap), concsat(nvap)
      character (len=12) :: nuc_name(nvap)

      end module mod_nuc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine nuc_init (Napari,Lovejoy,iod)

! Purpose:
! --------
!     This SR checks consistency between several options, and gets indexes of
!     relevant species. The user has to declare the names of the species used in
!     the nucleation code.

! History:
! --------
!     30-11-2016   J. Bock   first version of this SR



      USE gas_common, ONLY:
! Imported Parameters:
     &     j1,
     &     j5,
! Imported Array Variables with intent (in):
     &     gas_name, gas_mass,
     &     rad_name, rad_mass

      USE mod_nuc, ONLY:
! Imported Parameters:
     &     nvap,
! Imported Scalar Variables with intent (out):
     &     ind_H2SO4, ind_NH3, ind_OIO,
! Imported Array Variables with intent (out):
     &     ical,
     &     ivap,
     &     m_vap,
     &     concsat,
     &     nuc_name

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      logical, intent(in) :: Napari, Lovejoy, iod

! Local scalars:
      integer :: jvap, jspec

!- End of header ---------------------------------------------------------------

! ==============================================================================
!        USER DEFINED LIST OF GASES HAS TO BE DEFINED BELOW
! ==============================================================================

!   --- Prescribe "initial" concentration(s) (in ppt) in gas_species.csv ---
!   --- also, NUCV (i.e. hypothetical unreactive condensible vapor, gas 75) can be used ---
!              ivap(i): condensible vapors, please adjust kind and number!
!              ical(i): radical(s3) =1, no radical(s1) =0, please adjust!
!   --- saturation vapor density in [molec/cm3] (=0. for non-volatile vapors), please adjust!
!       note: ideally, non-zero concsat should be similar in magnitude to avoid errors!

      jvap = 0

      jvap = jvap+1 ; nuc_name(jvap) = 'OIO'   ; ical(jvap) = 1  ! radical
         concsat(jvap) = 0.
!      jvap = jvap+1 ; nuc_name(jvap) = 'I2O2'  ; ical(jvap) = 0  ! non radical
!         concsat(jvap) = 0.
!      jvap = jvap+1 ; nuc_name(jvap) = 'IO'    ; ical(jvap) = 1  ! radical
!         concsat(jvap) = 0.
!      jvap = jvap+1 ; nuc_name(jvap) = 'HOI'   ; ical(jvap) = 0  ! non radical
!         concsat(jvap) = 0.
!      jvap = jvap+1 ; nuc_name(jvap) = 'H2SO4' ; ical(jvap) = 0  ! non radical
!         concsat(jvap) = 0.

      if (jvap /= nvap) then
         print*,'Error in SR nuc_init:'
         print*,'  the number of condensible vapours (nvap = ',nvap,
     &          ') differs from the number of species names (jvap = ',
     &          jvap,')'
         print*,'  Please adjust.'
         stop 'Stopped by SR nuc_init'
      end if

! ==============================================================================
!   Retrieve the species names from gas (radical or non radical) name lists
! ==============================================================================
      do jvap = 1,nvap
         jspec = 1 ! initialise index
         if (ical(jvap) == 0) then
            do while ( trim(gas_name(jspec)) /= trim(nuc_name(jvap)) )
               jspec = jspec+1
               if (jspec > j1) then
                  print*,'Error in SR nuc_init:'
                  print*,'  ',trim(nuc_name(jvap)),' not found in the '
     &                   //'list of gas names'
                  stop 'Stopped by SR nuc_init'
               end if
            end do
            ivap(jvap) = jspec
            m_vap(jvap) = gas_mass (jspec)

         else if (ical(jvap) == 1) then
            do while ( trim(rad_name(jspec)) /= trim(nuc_name(jvap)) )
               jspec = jspec+1
               if (jspec > j5) then
                  print*,'Error in SR nuc_init:'
                  print*,'  ',trim(nuc_name(jvap)),' not found in the '
     &                   //'list of radicals names'
                  stop 'Stopped by SR nuc_init'
               end if
            end do
            ivap(jvap) = jspec
            m_vap(jvap) = rad_mass (jspec)

         else
            print*,'Error in SR nuc_init:'
            print*,'  ical must be 0 or 1, ical = ',ical(jvap),jvap
            stop 'Stopped by SR nuc_init'
         end if

      end do

! ==============================================================================
!   Retrieve specific indexes depending on the options
! ==============================================================================
      if (Napari) then ! Search for H2SO4 and NH3 indexes
         jspec = 1              ! initialise index
         do while ( trim(gas_name(jspec)) /= 'H2SO4' )
            jspec = jspec+1
            if (jspec > j1) then
               print*,'Error in SR nuc_init:'
               print*,'  H2SO4 not found in the list of gas names'
               stop 'Stopped by SR nuc_init'
            end if
         end do
         ind_H2SO4 = jspec

         jspec = 1              ! initialise index
         do while ( trim(gas_name(jspec)) /= 'NH3' )
            jspec = jspec+1
            if (jspec > j1) then
               print*,'Error in SR nuc_init:'
               print*,'  NH3 not found in the list of gas names'
               stop 'Stopped by SR nuc_init'
            end if
         end do
         ind_NH3 = jspec
      else
         ind_H2SO4 = 0
         ind_NH3   = 0
      end if

      if (Lovejoy) then ! Search for OIO index
         if (.not.iod) then
            print*,'Error: inconsistency in the nucleation module'
            print*,'  Lovejoy scheme require OIO, thus iod must be true'
            stop 'Stopped by SR nuc_init'
         end if

         jspec = 1              ! initialise index
         do while ( trim(rad_name(jspec)) /= 'OIO' )
            jspec = jspec+1
            if (jspec > j5) then
               print*,'Error in SR nuc_init:'
               print*,'  OIO not found in the list of gas names'
               stop 'Stopped by SR nuc_init'
            end if
         end do
         ind_OIO = jspec
      else
         ind_OIO = 0
      end if

!      print*,'ind_H2SO4, ind_NH3, ind_OIO',ind_H2SO4, ind_NH3, ind_OIO
!      write(*,1011)(jvap,ivap(jvap),ical(jvap),m_vap(jvap),jvap=1,nvap)
! 1011 format(3i4,f6.3)


      end subroutine nuc_init
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine appnucl2 (dt,both)
!
! Description:
!     help routine for apparent nucleation if apparent nucleation
!     has to be called twice because of two different real nucleation mechanisms
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header, implicit none         <Josue Bock>
!                    Use module for parameters
!
! 1.0       ?        Original code, Mistra v7.4.0  <Susanne Pechtl (b. Marquart)>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision dt
      logical both

! Local scalars:
      integer k

! Local arrays:
      double precision xn_app1(n), dnucv1(n), grorate1(n), grorate2(n),
     &                 concnuc1(n), concnuc2(n)

! Common blocks:
      common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n),
     &                 bd_mean(n), dnucv(n), grorate(n), concnuc(n)
      double precision xn_app, xn_apacc, xv_apacc, bn_ges,
     &                 bd_mean, dnucv, grorate, concnuc
!- End of header ---------------------------------------------------------------



      call appnucl (dt,.true.,.false.,both)

      do k = 2, n-1
         xn_app1(k)  = xn_app(k)
         dnucv1(k)   = dnucv(k)
         grorate1(k) = grorate(k)
         concnuc1(k) = concnuc(k)
      enddo

      call appnucl (dt,.false.,.true.,both)

      do k = 2, n-1
         xn_app(k) = xn_app(k) + xn_app1(k)
         dnucv(k)  = dnucv(k)  + dnucv1(k)
         grorate2(k) = grorate(k)
         grorate(k)  = (grorate1(k) + grorate2(k)) /2.
         concnuc2(k) = concnuc(k)
         if (concnuc2(k).ge.concnuc1(k)) then
            if ((xn_app(k)-xn_app1(k)).gt.0.01) then
               concnuc(k) = concnuc2(k)*xn_app(k)/(xn_app(k)-xn_app1(k))
            else
               concnuc(k) = concnuc2(k)
            endif
         else
            concnuc(k) = concnuc1(k) * xn_app(k)/xn_app1(k)
         endif
      enddo

      end subroutine appnucl2
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine appnucl (dt,Napari,Lovejoy,both)
!
! Description:
!     Calculation of the "apparent" nucleation rate of non-valotile vapors
!     after Kerminen and Kulmala (2002) @149.
!     With amendments for semivolatile vapors after Kerminen et al. (2004) @429
!
! Method:
!     assumptions of Kerminen and Kulmala (2002):
!     (1) only important sink for nuclei is coagulation to
!     larger pre-existing particles (i.e., no self-coagulation!)
!     (2) nuclei grow by condensation at a constant rate
!     (3) pre-existing particle population remains unchanged during nuclei growth
!     main subsequent restrictions:
!     total nuclei number concentration must remain sufficiently low
!     to prevent effective self-coagulation (< 10^5-10^6 nuclei/cm3)

!     modifications of Kerminen et al. (2004):
!     First, nuclei grow at a constant rate by condensation of non-volatile vapors.
!     After the nucleus has reached a critical size, it grows by non-volatile
!     and semi-volatile vapors (i.e. two different growth rates).
!
!     switches that must be set in str.f:
!       Napari:  optional coupling with ternary H2SO4-H2O-NH3 nucleation
!       Lovejoy: optional coupling with homogeneous OIO nucleation
!       ifeed:   optional feedback of nucleation with microphysics/chemistry
!              ifeed = 0 !no feedback to background particles
!              ifeed = 1 !with feedback to background particles
!              ifeed = 2 !feedback to background particles, but particles do not grow
!                         out of the 1st microphysical bin and chemistry does not "see"
!                         the 1st microphysical bin (i.e. no feedback nucl.->chemistry)
!
!     if chosen (Napari = .true.):
!     "Real" ternary H2SO4-H2O-NH3 nucleation rate and cluster size can be calculated
!     after Napari et al. (2002) @116
!     validity ranges: temp = 240 - 300 K
!                      RH = 0.05 - 0.95
!                      NH3 = 0.1 - 100 ppt
!                      H2SO4 = 1e4 - 1e9 molec./cm3
!                      J_real = 1e-5 - 1e6 /(cm3s)
!     if chosen (Lovejoy = .true.):
!     "Real" homogeneous OIO nucleation rate from data of Lovejoy/Burkholder
!
!     Output: J_real und d_nucini (= Input fuer APPNUCL)
!
!     to adjust in nuc.f: nvap, ivap, ical, m_vap, conc, concsat, (fcs), (xnue),
!                           d_nucini, J_real, uppest nucleation level (k-loop!)
!
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.3      11/2016   bugfix: transfer gas ==> liq was nested in if construct, thus
!                            H2SO4 was not transferred properly
!
! 1.2      07/2016   Header, cleaning                               <Josue Bock>
!                    Use module for parameters
!                    Explicit declaration of Knnuc, Kncrit
!                    Declarations (with comments) and implicit none
!
! 1.1       ?                                                       <Susanne Pechtl (b. Marquart)>
!
! 1.0       ?        Original code, Mistra v7.3.2                   <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE constants, ONLY :
     &     conv1, ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol
     &     pi

      USE gas_common, ONLY:
! Imported Parameters:
     &     j1,
     &     j5,
! Imported Array Variables with intent (inout):
     &     s1,
     &     s3

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nkc,
     &     nka,
     &     nkt,
     &     n

      USE mod_nuc, ONLY:
! Imported Parameters:
     &     nvap,
     &     ical,
     &     ivap,
     &     nuc_name,
     &     m_vap,         ! Molar weight of nucleation vapor (e.g. NUCV=OIO) [kg/mol]
     &     concsat        ! saturation vapor density of condensible vapor [/cm3]
!                           (measure for saturation vapor pressure over particle surface)

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision dt
      logical Napari        ! Napari=true: J_real and d_nucini are calculated from ternary nucleation
      logical Lovejoy       ! Lovejoy=true: J_real and d_nucini are calculated from OIO nucleation
      logical both          ! true if both Napari and Lovejoy are true (needed)


!     --- Local Variables ---
      integer i,ia,icount,iv,j,jt,jts,jtt,k
      integer isum, mult

      double precision alphaa !accomodation coefficient of condensing vapor
      double precision betanuc, betacrit
      double precision ro_nuc !nuclei density [kg/m3]
      double precision temp, press !Temperature [K] and Pressure [Pa]
      double precision rh     ! relative humidity in layer k
      double precision dp(nkt)  !Particle diameter [nm]
      double precision dpmin    !Particle diameter of smallest dry size bin [nm]
      double precision dpmint   !Particle diameter of according total size bin [nm]
      double precision Np(nkt)  !Particle number [/cm3]
      double precision Nges     !total particle number [/cm3] (calculated)
      double precision lambda   !Mean free path [m]
      double precision v_mean(nvap) !Mean molecular speed of vapor [m/s]
      double precision d_nucini !initial particle diameter [nm]; from subroutine ternucl
      double precision J_real   !real nucleation rate [/(cm3*s)]; from subroutine ternucl
      double precision d_mean   !Number mean diameter [nm]; (calculated)
      double precision conc(nvap)    !concentration of condensible vapor [/cm3]
      double precision conc_nuc
      double precision concsemi      !sum of all semi-volatile vapors
      double precision concsatmean   !mean saturation vapor density for semi-volatile vapors
      double precision Scrit !Ratio of supersaturation of semi-volatile vapors
      double precision dcrit !Critical nuclei diameter from which particle grows also
!                             by condensation of semivolatile vapors [nm]
!
      double precision Kn(nkt), beta(nkt), cs,gr,grs,gamma,eta,etages !calculated
      double precision J_app    !apparent nucleation rate [/(cm3*s)] (main output)
      double precision Knnuc, Kncrit ! jjb explicit definition needed for these two variables
!                                          whose names start with a "K"
      double precision factor
      double precision a0, a0mn, b0, b0mn
      double precision deltax, facul, gg, over, rmin, summ 
      double precision rg !equilibrium total aerosol radius
!     --- necessary for MISTRA mass bilances: ---
      double precision rel_vap(nvap) !rel. contribution of vapor to nuclei growth
      double precision m_vapmean     !mean molar weight of vapor
      double precision s1old(j1,n), s3old(j5,n)   ![mol/m3]
      double precision soldsum(n) ![mol/m3]

! Common blocks:
      common /backpart/ partd(nkt,n), partN(nkt,n), partsa(n)
      double precision partd, partN, partsa

      common /blck01/ am3(n),cm3(n)
      double precision am3, cm3

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      double precision ff, fsum
      integer nar

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      double precision xm1, xm2, feu, dfddt, xm1a, xm2a

      common /nucfeed/ ifeed
      integer ifeed          !ifeed = 0/1/2: feedback with background particles no/yes/partly

      common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n),
     &                 bd_mean(n), dnucv(n), grorate(n), concnuc(n)
      double precision xn_app, xn_apacc, xv_apacc, bn_ges,
     &                 bd_mean, dnucv, grorate, concnuc

      common /nucl2/ Jn, dc      ! input from real nucleation
      double precision
     &     Jn,       ! nucleation rate [1/(cm3 s)]
     &     dc        ! diameter of cluster [nm]

      common /nucl3/ Jnio, dcio  ! input from real nucleation
      double precision 
     &     Jnio,                ! nucleation rate [1/(cm3 s)]
     &     dcio                 ! diameter of cluster [nm]

      common /nucsum/ ssum(n)
      double precision ssum

! External function:
      double precision rgl
      external rgl
!- End of header ---------------------------------------------------------------


      open(unit=20,file='nuc.out',status='unknown',form='formatted')

      dpmin = rn(1) * 2000.d0
      ro_nuc = 2000.       !nuclei density [kg/m3] (same as assumed aesosol density rho3)

!     --- alphaa: for values for some species see alpha(k,spec) in subroutine st_coeff_t ---
!         but: in the current version only 1 alphaa can be used for all species together
      alphaa = 1.         !mass absorption coefficient of condensing vapor


      do k=2,n-1  !Loop over vertical levels (does not necessariliy have to extend to n-1!)

        temp = t(k)
        press = p(k)
        rh = feu(k)
        if (rh.ge.1.) rh = 0.999 !for equilibrium growth (see function rgl)

        ! Initialise concentration from non radical (s1) or radical (s3) gas concentration
        do i=1,nvap
           if(ical(i)==0) then
              conc(i)=s1(ivap(i),k) * conv1
           else
              conc(i)=s3(ivap(i),k) * conv1
           end if
        end do

!     --- 1D size distribution of background particles wrt total droplet diameter ---
        write (20,1070)
        partsa(k) = 0.
        do jt = 1,nkt
          dp(jt) = 2000. * rq(jt,1) !first dry aerosol class defines size bins for 1D distribution
          Np(jt) = 0.
          do ia = 1, nka
            if (rn(ia).gt.rw(jt,1)) goto 2001
            do jtt = 1,nkt
              if (rq(jtt,ia).le.rw(jt,1)) then
                if (jt.gt.1) then
                  if (rq(jtt,ia).gt.rw(jt-1,1))
     &               Np(jt) = Np(jt) + ff(jtt,ia,k)

!                else if ((jt.eq.1).and.(rq(jtt+1,ia).gt.rw(jt,1))) then   ! jjb index out of bounds if jtt = nkt
                else if ((jt.eq.1).and.(jtt.lt.nkt).and.                   ! jjb quick fix, CHECK!
     &                    (rq(jtt+1,ia).gt.rw(jt,1))) then
                  Np(jt) = Np(jt) + ff(jtt,ia,k)
                endif
              else
                goto 2002
              endif
            enddo
 2002       continue
          enddo
 2001     continue
          write (20,1080) k, jt, dp(jt), Np(jt)
!         -- background particle surface area --
          partsa(k) = partsa(k) + Np(jt) *dp(jt)**2 *3.1416* 1.d-6 !in um2/cm3
!         -- 1D particle distribution --
          partN(jt,k) = Np(jt)
          partd(jt,k) = dp(jt)
        enddo  !jt
 1070   format ('k  jt  dp(jt)   Np(jt)')
 1080   format (2I4,2F16.5)

!     --- lambda = freep(k) in kpp.f ---
        lambda = 2.28e-5 * temp / press

!     --- Initialize ---
        cs = 0.
        gr = 0.
        grs = 0.

        do i = 1,nvap
!     --- vmean = sqrt((8*R_gas*T/(M*pi)) ---
!     --- sqrt(8*R_gas/pi)=4.60138 ---
          v_mean(i) = sqrt(temp/m_vap(i)) * 4.60138
        enddo

! Note: both nucleation routines are called any time (for output reasons),
!       but nucleation is "physically seen" only if respective switches are .true.
!       I.i., if a switch is .false., nuclei are not seen by the model and
!       concentrations of nucleating species do not change.

!     --- ternary H2SO4-H2O-NH3 nucleation routine: ---
        if ((both).and.(.not.Napari).and.(Lovejoy)) then
!         do not call ternucl
        else
          call ternucl(dt,k,Napari)
        endif

!     --- homogeneous OIO nucleation routine: ---
        if ((both).and.(Napari).and.(.not.Lovejoy)) then
!         do not call oionucl
        else
          call oionucl(dt,k,Lovejoy)
        endif

        if (Napari) then
          d_nucini = dc
          J_real = Jn
        else if (Lovejoy) then
          d_nucini = dcio
          J_real = Jnio
        else !prescribe values for clusters
          d_nucini = 1.0       ![nm]          please adjust!
          J_real = 1000.      ![/(cm3*s)]     please adjust!
        endif

!     --- semi-volatile vapors ---
!     --- concsemi can be used as switch for the existence of semi-volatile vapors ---
        concsemi = 0. !sum of all semi-volatile vapors
        concsatmean = 0. !mean saturation vapor density for semi-volatile vapors
        icount = 0
        do i=1, nvap !mean values for semi-volatile vapors
          if (concsat(i).gt.0.) then
            concsemi = concsemi + conc(i)
            concsatmean = concsatmean + concsat(i)
            icount = icount + 1
          endif
        enddo
        if (icount .ge. 1) concsatmean = concsatmean / real(icount)
!     --- supersauration ratio and critical diameter for semi-volatile vapors after @429 ---
        if (concsatmean.ne.0.)  then
          Scrit = concsemi/concsatmean
          if (Scrit.le.1) then !no growth for subsaturation
            do i=1, nvap
              if (concsat(i).ne.0.) conc(i)=0.
            enddo
              dcrit = dpmin
          else
            dcrit = ( 6.49 - 0.01556*temp + 0.039*log(Scrit) ) /
     $              ( 1. - 0.002*temp + 0.174*log(Scrit) )
            if (dcrit.le.d_nucini) concsemi = 0. !semi-volatiles are defined as non-volatile
          endif
        endif

        call dmean(nkt,dp,Np,d_mean,Nges)

        do j = 1,nkt
!     --- Knudsen number ---
          Kn(j) = 2 * 1.d9 * lambda / dp(j)
!     --- transition correction for condensational mass flux (Fuchs & Sutugin, 1971) ---
          beta(j) = (1.+Kn(j)) /
     $              (1.+0.377*Kn(j)+1.33*Kn(j)*(1.+Kn(j))/alphaa)
!     --- Condensation sink (condensation of vapor on pre-existing particles) ---
          cs = cs + 0.5 * dp(j) * 1.d-7 * beta(j) * Np(j) ![1/cm2]
        enddo

!     --- growth rate of nuclei: change of diameter with time (eq.(20) of @149) ---
!     --- gr: growth rate for all volatile vapors;
!     --- grs: growth rate for all semi-volatile vapors
        Knnuc = 2 * 1.d9 * lambda / d_nucini
        betanuc = (1.+Knnuc) /
     $            (1.+0.377*Knnuc+1.33*Knnuc*(1.+Knnuc)/alphaa)
        if (concsemi.ne.0.) then !semi-volatiles
          Kncrit = 2 * 1.d9 * lambda / dcrit
          betacrit = (1.+Kncrit) /
     $            (1.+0.377*Kncrit+1.33*Kncrit*(1.+Kncrit)/alphaa)
        endif
        do i = 1,nvap
          if ((concsat(i).ne.0.).and.(concsemi.ne.0.)) then !semi-volatiles
            grs = grs + v_mean(i) * m_vap(i) * (conc(i)-concsat(i))
          else !non-volatiles
            gr = gr + v_mean(i) * m_vap(i) * conc(i)
          end if
        enddo
        gr= gr * 7969.45 * lambda * betanuc / d_nucini / ro_nuc   ![nm/h], non-volatiles
        if (concsemi.ne.0.)
     $    grs= grs * 7969.45 * lambda * betacrit / dcrit / ro_nuc   ![nm/h], semi-volatiles

!     ---  relative contribution of vapor i to new aerosol mass ---
!          (calculated according to definition of gr and grs)
!          factor determines mass fraction of semi-volatile vapors
        m_vapmean = 0.
        if (concsemi.ne.0.) then
          factor = grs * (dpmin**3-dcrit**3) /
     $     ( grs * (dpmin**3-dcrit**3) + gr * (dpmin**3-d_nucini**3) )
        else
          factor = 0.
        endif
        do  i = 1,nvap
          if ((concsat(i).ne.0.).and.(concsemi.ne.0.)) then !semi-volatiles
            if (grs .gt. 1.d-2) then
              rel_vap(i) = (v_mean(i) * m_vap(i)* (conc(i)-concsat(i)))
     $                / (grs *ro_nuc *dcrit /betacrit /lambda /7969.45)
     $                * factor
            else
              rel_vap(i) = 1. / REAL(nvap)
            endif
          else
            if (gr .gt. 1.d-2) then
              rel_vap(i) = (v_mean(i) * m_vap(i)* conc(i))
     $               / (gr *ro_nuc *d_nucini /betanuc /lambda /7969.45)
     $               * (1. - factor)
            else
              rel_vap(i) = 1. / REAL(nvap)
            endif
          endif
          m_vapmean = m_vapmean + rel_vap(i) * m_vap(i)
        enddo

!     --- calculate equilibrium with ambient water vapor (RH) ---
!     parameters a0m, b0m of koehler curve of subroutine subkon:
!     sr(ia,jt)=exp(a0m/(r(ia,jt)*t)-b0m(ia)*en(ia)/ew(jt)))
!     a0m see p. 141 pruppacher and klett a0m=2 sigma/(r1*t*rhow*a)
!     152200= 2 sigma*10**6 with sigma is surface tension = 76.1*10**-3
!     see pruppacher and klett p. 104
!     a0mn = a0m=152200./(r1*rhow)
!     b0mn = b0m=fcs*xnue*xmol2/xmol3: fcs(ia) fraction of soluble species
!     xnue number of ions; xmol2 (xmol3) mol masses of water (aerosol)
        a0mn = 152200. / (461.51 * ro_nuc)
        b0mn = 1. * 1. * 0.018/m_vapmean  !fcs=1 and xnue=1 are assumed -> please adjust!
        a0 = a0mn/temp
!       ---  b0=b0m*rho3/rhow; rho3=2000; rhow=1000
        b0 = b0mn*2.
        rmin = dpmin / 2000.
        rg = rgl(rmin,a0,b0,rh) !equilibrium total aerosol radius
        do jt = 1,nkt !find correct droplet size bin jts to smallest dry bin
          if (rw(jt,1) .ge. rg) then
            jts = jt
            goto 2003
          endif
        enddo
 2003   continue
        dpmint = rq(jts,1) * 2000.d0
!     --- Update growth rate (including equilibrium with ambient RH) after Kerminen (pers. comm.)---
        gr = gr * dpmint/dpmin
        grs = grs * dpmint/dpmin

!     --- semi-empirical proportionality factor
!     (describing coagulation with larger pre-existing particles)
        gamma = 2300.d0 *(d_nucini/1.d0)**0.2d0 *(dpmint/3.d0)**0.075d0
     $       * (d_mean/150.d0)**0.048d0 * (ro_nuc/1000.d0)**(-0.33d0)
     $       * (temp/293.d0)**(-0.75d0) ![nm2cm2/h]

!     --- Apparent nucleation rate for above calculated size bin of total diameter ---
        if ((gr .gt. 1.d-2).and.(J_real .gt. 0.01)) then
          eta = gamma * cs / gr   ![nm]
          etages = gamma * cs / (gr + grs)  ![nm]
          if (concsemi.ne.0.) then
            J_app = J_real *
     $        exp (etages/dpmint + (eta-etages)/dcrit - eta/d_nucini) ![molecules/(cm3*s)]
          else
            J_app = J_real * exp (eta/dpmint - eta/d_nucini)
          endif
!         --- diagnose total nuclei concentration < dpmint ---
!             in order to check the validity range (19) in @149
!             conc_nunc should by < 1.e6 nuclei/cm3
!             (exactly valid only if no semi-volatile vapors are present)
          SUMM = 0.
          do isum = 1,20
            facul = 1.
            mult = isum
            do while (mult.GE.1)  !Calculation of k!
               facul = facul * mult
               mult = mult - 1
            enddo !while
            SUMM = SUMM + ( eta**(isum+1) * (d_nucini**(-isum)-
     $             dpmint**(-isum)) ) / isum/facul
          enddo
! avoid floating overflow (exp(700) =(approx) 1.d308 which is highest number on manolito):
! (this is output only)
          over=0.
          if (eta/dpmint.gt.700.) over=1.
          if (eta/d_nucini.gt.700.) over=1.
          if (over.eq.0.)
     $        GG = dpmint*exp(eta/dpmint) - d_nucini*exp(eta/d_nucini)
     $             + eta*log(dpmint/d_nucini) + SUMM
          conc_nuc = J_app/gr * exp (-eta/dpmint) * GG * 3600.  ![nuclei/cm3]
        else
          eta = 0.
          etages = 0.
          J_app = 0.
          conc_nuc = 0.
        endif

!     --- Update of variables ---
        xn_app(k) = 0.
        bn_ges(k) = Nges
        bd_mean(k) = d_mean
        soldsum(k) = 0.
        concnuc(k) = 0.

        do iv = 1,nvap
          if (ical(iv).eq.0) then      ! non radical gas phase species
            s1old(ivap(iv),k) = s1(ivap(iv),k)
            soldsum(k) = soldsum(k) + s1old(ivap(iv),k)
          else if (ical(iv).eq.1) then ! radical gas phase species
            s3old(ivap(iv),k) = s3(ivap(iv),k)
            soldsum(k) = soldsum(k) + s3old(ivap(iv),k)
          endif
        enddo !iv

        ssum(k) = soldsum(k)

        if (J_app .gt. 0.1) then
          xn_app(k) = J_app  !nucleation rate of new "apparent" particles [/(cm3*s)]
          concnuc(k) = conc_nuc ! 1D (real) nuclei concentration for output
          if (ifeed.ne.0) ff(jts,1,k) = ff(jts,1,k) + xn_app(k)*dt
          deltax = xn_app(k)*dt * pi/6. * (dpmin**3-d_nucini**3)
     $           * ro_nuc / m_vapmean * 1.d-21    ![mol/m3(air)] !additional dry material
          ssum(k) = 0.
          do iv = 1,nvap
            if (ical(iv).EQ.0) then
              s1(ivap(iv),k) = s1(ivap(iv),k) - deltax*rel_vap(iv)
              if (s1(ivap(iv),k) .lt. 0.d0) s1(ivap(iv),k) = 0.d0
              ssum(k) = ssum(k) + s1(ivap(iv),k)
            else
              s3(ivap(iv),k) = s3(ivap(iv),k) - deltax*rel_vap(iv)
              if (s3(ivap(iv),k) .lt. 0.d0) s3(ivap(iv),k) = 0.d0
              ssum(k) = ssum(k) + s3(ivap(iv),k) !ssum = sum of nucleating vapors
            end if
! OIO goes into liquid phase (for mass conserving reasons only)
! note: other condensing vapors than OIO and H2SO4 are NOT conserved!
!    OIO goes  to unreactive OIO in the liquid phase
            if (trim(nuc_name(iv)).eq.'OIO')
     &            sl1(j1+12,1,k) = sl1(j1+12,1,k) +
     &                             s3old(ivap(iv),k) - s3(ivap(iv),k)
!    H2SO4 goes to SO4 in the liquid phase (i.e., equilibrium HSO4m and Hp)
            if (trim(nuc_name(iv)).eq.'H2SO4')
     &            sl1(6,1,k) = sl1(6,1,k) +
     &                         s1old(ivap(iv),k) - s1(ivap(iv),k)
          enddo !iv

          xn_apacc(k) = xn_apacc(k) + xn_app(k)*dt !accumulated number [/cm3]
          xv_apacc(k) = xv_apacc(k) +
     $         pi/6. * dpmint**3 *xn_app(k)*dt ! acculumated volume [nm3/cm3]
          write (20,109)
          write (20,110) k,jts,xn_apacc(k),xn_app(k),ssum(k)
        endif
!        dnucv(k) = soldsum(k) - ssum(k)  !change in nucleationg vapors (sum) [mol/m3]
        dnucv(k) = (soldsum(k) - ssum(k)) / am3(k) *1.d12  !change in nucleationg vapors (sum) [ppt]
        grorate(k) = gr   ! 1D growth rate for output (for non-volatiles)

!     --- Output (Test) ---
        write(20,1003)
        do i=1,nvap
          if (j_real.ne.0.) then
             write(20,1004) k,d_nucini,dpmin,dpmint,d_mean,Nges,temp,
     $            press,rh,ro_nuc,m_vap(i),conc(i),cs,gr,eta,gamma,
     $            J_real,J_app,(J_app/J_real)
          else
             write(20,1004) k,d_nucini,dpmin,dpmint,d_mean,Nges,temp,
     $            press,rh,ro_nuc,m_vap(i),conc(i),cs,gr,eta,gamma,
     $            J_real,J_app
          endif
        enddo

      enddo !Loop over vertical levels (k)

 1003 format ('k,d_nucini(nm),dpmin(nm),dpmint(nm),d_mean(nm),Nges, '
     $ // 'temp(K),press(Pa),rh,ro_nuc(kg/m3),m_vap(i)(kg/mol), '
     $ // 'conc(i)(/cm3),cs(/cm2),gr(nm/h),eta(nm),gamma(nm2cm2/h), '
     $ // 'J_real(/cm3s),J_app(/cm3s),(J_app/J_real)')
 109  format ('k,jts,xn_apacc(k),xn_app(k),ssum(k)')
 1004 format (I4,10F14.3,E12.5,7F14.5)
 110  format (2i4,9d16.8)

      close (20)

      end subroutine appnucl
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine dmean (n,dp,Np,d_mean,Nges)
!
! Description:
!     --- This routine calculates the number mean diameter d_mean ---
!         and the total particle number Np_ges of the background particles
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      07/2016   Header, implicit none          <Josue Bock>
!
! 1.1       ?        Original code.                 <Susanne Pechtl (b. Marquart)>
!
! 1.0       ?        Original code, Mistra v7.3.2   <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer n
! Array arguments with intent(in):
      double precision dp(n), Np(n)
! Scalar arguments with intent(out):
      double precision d_mean, Nges
! Local scalars:
      integer i
!- End of header ---------------------------------------------------------------

      d_mean = 0.
      Nges = 0.
      do i=1,n
        d_mean = d_mean + dp(i)*Np(i)
        Nges = Nges + Np(i)
      enddo
      if (Nges.gt.0.) then
         d_mean = d_mean / Nges
      else
         d_mean = 0.
      endif

      end subroutine dmean
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine ternucl (dt,k,Napari)
!
! Description:
!    calculation of nucleation rate and cluster size
!    using the parameterization of Napari et al., JGR, 107, 4381,
!    doi:10.1029/2002JD002132, 2002
!    validity range: h2so4: 1.d4 - 1.d9 molec/cm3
!                    nh3: 0.1 - 100 ppt
!                    J_n: 1.d-5 - 1.d6 /(cm3 s)
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      07/2016   Header, implicit none                     <Josue Bock>
!                    Use modules for constant and parameters
!
! 1.1       ?                                                  <Susanne Pechtl (b. Marquart)>
!
! 1.0       ?        Original code, Mistra v7.1.1              <Roland von Glasow>
!
! 02-Mar-2017  <Josue Bock> bugfix: missing initialisation of nn and nh

! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE constants, ONLY :
     &     conv1 ! multiply by conv1 to convert cm^3(air)/mlc --> m^3(air)/mol

      USE gas_common, ONLY: s1

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      USE mod_nuc, ONLY:
! Imported Parameters:
     &     ind_H2SO4,
     &     ind_NH3

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision dt
      integer k
      logical Napari

! Local scalars:
      double precision
     &     c_h2so4,
     &     c_nh3,
     &     nh,        ! number of H2SO4 molecules in critical cluster
     &     nn,        ! number of NH3 molecules in critical cluster
!     &     n_tot,nt, ! total number of molecules in critical cluster
     &     rc,        ! radius of cluster [nm]
     &     rh,
     &     Temp

! Internal functions:
      double precision 
     &     a,b,       ! arguments
     &     n_h2so4,   ! number of H2SO4 molecules in critical cluster
     &     n_nh3,     ! number of NH3 molecules in critical cluster
     &     r_clust    ! radius of cluster [nm]

! External function:
      double precision, external ::
     &     J_nuc      ! nucleation rate [1/(cm3 s)]

! Common blocks:
      common /blck01/ am3(n),cm3(n)
      double precision am3, cm3

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      double precision xm1, xm2, feu, dfddt, xm1a, xm2a

      common /nucl/ xn_new(n), xn_acc(n), xv_acc(n), dh2so4(n), dnh3(n)
      double precision xn_new, xn_acc, xv_acc, dh2so4, dnh3

      common /nucl2/ Jn, dc      !SM: Jn and dc are used in subroutine appnucl
      double precision
     &     Jn,       ! nucleation rate [1/(cm3 s)]
     &     dc        ! diameter of cluster [nm]
!- End of header ---------------------------------------------------------------

! functions for diagnostics of critical cluster; a = J_nuc, b=Temp
      n_h2so4 (a,b) = 38.1645  + 0.774106*log(a)  + 2.98879d-3*log(a)**2
     &     - 0.357605*b  - 3.66358d-3*log(a)*b + 8.553d-4  *b**2
      n_nh3 (a,b)   = 26.8982  + 0.682905*log(a)  + 3.57521d-3*log(a)**2
     &     - 0.265748*b  - 3.41895d-3*log(a)*b + 6.73454d-4*b**2
!      n_tot (a,b)   = 79.3484  + 1.7384  *log(a)  + 7.11403d-3*log(a)**2
!     &     - 0.744993*b  - 8.20608d-3*log(a)*b + 1.7855d-3 *b**2
      r_clust (a,b) = 0.141027 - 1.22625d-3*log(a)- 7.82211d-6*log(a)**2
     &     - 1.56727d-3*b- 3.076d-5  *log(a)*b + 1.08375d-5*b**2

! convert units
         rh      = feu(k)                  ! relative humidity [0..1]
         c_nh3   = s1(ind_NH3,k) / am3(k) *1.d12 ! NH3 in [pmol/mol]
         c_nh3   = min(100.,c_nh3)         ! NH3 <= 100 ppt
         c_h2so4 = s1(ind_H2SO4,k) * conv1         ! H2SO4 in [molec/cm3]
!        c_h2so4 = c_h2so4 * 1000.   !artificial amplification for test
         write (20,120)
         write (20,110) k,rh,c_nh3,c_h2so4
         Temp    = t(k)                      ! temperature [K]
         Jn = 0.
         nn = 0.
         nh = 0.
         xn_new(k) = 0.
         rc = 1. !initialize
         dc = 2.* rc
         if (c_h2so4.gt.1.d4) then
! calculate the nucleation rate
            Jn=min(1.d6,J_nuc(rh,c_nh3,c_h2so4,Temp))     ! nucleation rate [1/(cm3 s)]
            if (Jn .ge. 0.01) then
               nh=n_h2so4(Jn,Temp) ! number of H2SO4 molecules in critical cluster
               nn=n_nh3(Jn,Temp) ! number of NH3 molecules in critical cluster
!               nt=n_tot(Jn,Temp)   ! total number of molecules in critical cluster
               rc=r_clust(Jn,Temp) ! radius of cluster [nm]
               dc = 2.* rc         ! diameter of cluster [nm]
               xn_new(k) = Jn !nucleation rate of new clusters [/(cm3*s)]
               if (Napari) then
                 s1(ind_NH3,k)   = s1(ind_NH3,k) - xn_new(k)*dt*nn/conv1
                 if (s1(ind_NH3,k) .lt. 0.) s1(ind_NH3,k)=0.

                 s1(ind_H2SO4,k) = s1(ind_H2SO4,k)-xn_new(k)*dt*nh/conv1
                 if (s1(ind_H2SO4,k) .lt. 0.) s1(ind_H2SO4,k)=0.
               endif
               xn_acc(k) = xn_acc(k) + xn_new(k)*dt !accumulated number of new clusters [/cm3]
               xv_acc(k) = xv_acc(k) + 4.18879*rc**3*xn_new(k)*dt ! acculumated volume of new clusters [nm3/cm3]
               write (20,109)
               write (20,110) k,dc,nn,nh,xn_acc(k),xn_new(k)
            endif
         endif
!         dnh3(k) = (c_nh3*am3(k)/1.d12) - s1(ind_NH3,k)  ![mol/m3]
!         dh2so4(k) = (c_h2so4/conv1) - s1(ind_H2SO4,k)   ![mol/m3]
!         dnh3(k) = xn_new(k)*dt*nn/conv1           ![mol/m3]
!         dh2so4(k) = xn_new(k)*dt*nh/conv1         ![mol/m3]
         dnh3(k) = xn_new(k)*dt*nn/conv1 / am3(k) *1.d12 ![pmol/mol]
         dh2so4(k) = xn_new(k)*dt*nh                     ![molec/cm3]

! save xnew to some output array
 109  format ('k,dc,nn,nh,xn_acc(k),xn_new(k)')
 110  format (i4,9d16.8)
 120  format ('k,rh,c_nh3[ppt],c_h2so4[/cm3]')

      end subroutine ternucl

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine oionucl (dt,k,Lovejoy)
!
! Description:
!    calculation of nucleation rate from function
!    inferred from model simulation of homogeneous
!    OIO nucleation in the lab by Lovejoy/Burkholder (pers. comm.)
!
!    function is valid for nuclei with d=2.0nm and density=2g/cm3
!    (i.e. 1 nuclei has 34 OIO molecules)
!    validity range: 0.4 ppt < OIO < 22 ppt, extrapolation seems ok
!
!    it is assumed that only OIO contributes to nucleation
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.1      07/2016   Header                                                  <Josue Bock>
!                    Use module for constant and parameters
!                    Declarations (with comments) and implicit none
!                    Debuging: changed name of internal function J_nuc2
!                             which was conflicting with external fct J_nuc
!
! 1.0       ?        Original code, Mistra v7.4.0                             <Susanne Pechtl (b. Marquart)>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE constants, ONLY :
     &     conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

      USE gas_common, ONLY:
     &     s3

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      USE mod_nuc, ONLY:
     &     ind_OIO

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision dt
      integer k
      logical Lovejoy

! Local parameters:
      double precision n_oio     ! number of OIO molecules in critical cluster
      parameter ( n_oio=34. )
      double precision rc        ! radius of cluster [nm]
      parameter ( rc=1. )

! Local scalar:
      double precision 
     &     c_oio,                ! OIO in [pmol/mol]
     &     Temp                  ! temperature [K]

! Internal function:
      double precision J_nuc2
      double precision a, b      ! internal function arguments

! Common blocks:
      common /blck01/ am3(n),cm3(n)
      double precision am3, cm3

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      common /nuclio/ xn_newio(n), xn_accio(n), xv_accio(n), doio(n)
      double precision xn_newio, xn_accio, xv_accio, doio

      common /nucl3/ Jnio, dcio
      double precision 
     &     Jnio,                ! nucleation rate [1/(cm3 s)]
     &     dcio                 ! diameter of cluster [nm]
!- End of header ---------------------------------------------------------------


! Internal function for diagnostics of nucleation rate; a = oio[ppt], b=temp[K]
      J_nuc2(a,b) = a**(0.030657*b-4.4471) * exp(-0.30947*b+81.097)

! convert units
      c_oio   = s3(ind_OIO,k) / am3(k) *1.d12
      write (20,120)
      write (20,110) k,c_oio

      dcio = 2. * rc
      Temp = t(k)
      Jnio = 0.
      xn_newio(k) = 0.
      if (c_oio .gt. 0.01) then
! calculate the nucleation rate
         Jnio=min(1.d4,J_nuc2(c_oio,Temp))
         if (Jnio .ge. 0.01) then
            xn_newio(k) = Jnio    !nucleation rate of new clusters [/(cm3*s)]
            if (Lovejoy) then
               s3(ind_OIO,k) = s3(ind_OIO,k)-xn_newio(k)*dt*n_oio/conv1
               if (s3(ind_OIO,k) .lt. 0.d0) s3(ind_OIO,k)=0.d0
            endif
            xn_accio(k) = xn_accio(k) + xn_newio(k)*dt !accumulated number of new clusters [/cm3]
            xv_accio(k) = xv_accio(k) + 4.18879*rc**3*xn_newio(k)*dt ! acculumated volume of new clusters [nm3/cm3]
            write (20,109)
            write (20,110) k,dcio,n_oio,xn_accio(k),xn_newio(k)
         endif
      endif

!      doio(k) = (c_oio*am3(k)/1.d12) - s3(ind_OIO,k)   ![mol/m3]
!      doio(k) = xn_newio(k)*dt*n_oio/conv1             ![mol/m3]
      doio(k) = xn_newio(k)*dt*n_oio/conv1 / am3(k) *1.d12 ![pmol/mol]

! save xnew to some output array
 109  format ('k,dcio,n_oio,xn_accio(k),xn_newio(k)')
 110  format (i4,9d16.8)
 120  format ('k,c_oio[ppt]')

      end subroutine oionucl

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      double precision function J_nuc (rH,nh3,h2so4,Temp)
!
! Description:
!    nucleation rate
!    checked and compared to figs in #2158/@116: okay

! the parameterization is valid until for nucleation rates up to 10^6 1/(cm3 s), here
! it's used also for higher rates, so it's clearly an upper limit - magnitude of error
! not checked yet

! ln c  = ln [H2SO4], []=concentration
! ln S  = ln [NH3],   []=mixing ratio
! ln rH = ln rH
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      07/2016   Header, implicit none          <Josue Bock>
!                    Use module for parameters
!
! 1.1       ?                                       <Susanne Pechtl (b. Marquart)>
!
! 1.0       ?        Original code, Mistra v7.3.2   <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:

      implicit none

! Function arguments
! Scalar arguments with intent(in):
      double precision rh,nh3,h2so4,Temp
         
! Local scalars:
      double precision lnc, lnS, lnrH

! Local arrays:
      double precision fpd(4,20)

! Internal function:
      double precision fpol
      double precision T
      integer n

!- End of header ---------------------------------------------------------------

      data fpd  /-0.355297,    -3.38448d+1,   0.34536,     -8.24007d-4,
     &            3.13735,     -0.772861,     5.61204d-3,  -9.74576d-6,
     &            1.90359d+1,  -0.170957,     4.79808d-4,  -4.14699d-7,
     &            1.07605,      1.48932,     -7.96052d-3,   7.61229d-6,
     &            6.0916,      -1.25378,      9.39836d-3,  -1.74927d-5,
     &            0.31176,      1.64009,     -3.43852d-3,  -1.09753d-5,
     &           -2.00735d-2,  -0.752115,     5.25813d-3,  -8.98038d-6,
     &            0.165536,     3.26623,     -4.89703D-2,   1.46967D-4,
     &            6.52645,     -0.258002,     1.43456D-3,  -2.02036D-6,
     &            3.68024,     -0.204098,     1.06259D-3,  -1.26560D-6,
     &           -6.6514D-2,   -7.82382,      1.22938D-2,   6.18554D-5,
     &            0.65874,      0.190542,    -1.65718D-3,   3.41744D-6,
     &            5.99321D-2,   5.96475,     -3.62432D-2,   4.93337D-5,
     &           -0.732731,    -1.84179d-2,   1.47186d-4,  -2.37711d-7,
     &            0.728429,     3.64736,     -2.7422d-2,    4.93478d-5,
     &            4.13016d+1,  -0.35752,      9.04383d-4,  -5.73788d-7,
     &           -0.160336,     8.89881d-3,  -5.39514d-5,   8.39522d-8,
     &            8.57868,     -0.112358,     4.72626d-4,  -6.48365d-7,
     &            5.301767d-2, -1.98815,      1.57827d-2,  -2.93564d-5,
     &           -2.32736,      2.34646d-2,  -7.6519d-5,    8.0459d-8/

      fpol(n,T)=fpd(1,n) + fpd(2,n)*T + fpd(3,n)*T**2 + fpd(4,n)*T**3
      lnc  = log(h2so4)
      lnS  = log(nh3)
      lnrH = log(rH)

      J_nuc = exp(-84.7551 + fpol(1,Temp)/lnc + fpol(2,Temp)*lnc +
     &        fpol(3,Temp)*lnc**2 + fpol(4,Temp)*lnS +
     &        fpol(5,Temp)*lnS**2 + fpol(6,Temp)*rH +
     &        fpol(7,Temp)*lnrH + fpol(8,Temp)*lnS/lnc +
     &        fpol(9,Temp)*lnS*lnc + fpol(10,Temp)*rH*lnc +
     &        fpol(11,Temp)*rH/lnc + fpol(12,Temp)*rH*lnS +
     &        fpol(13,Temp)*lnrH/lnc + fpol(14,Temp)*lnrH*lnS +
     &        fpol(15,Temp)*lnS**2/lnc + fpol(16,Temp)*lnc*lnS**2 +
     &        fpol(17,Temp)*lnc**2*lnS +  fpol(18,Temp)*rH*lnS**2 +
     &        fpol(19,Temp)*rH*lnS/lnc + fpol(20,Temp)*lnc**2*lnS**2)

      end function J_nuc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine nucout1
!
! Description:
!    output nucleation: .asc output for gnu-plotting
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      07/2016   Header, implicit none          <Josue Bock>
!                    Use module for parameters
!
! 1.1       ?                                       <Susanne Pechtl (b. Marquart)>
!
! 1.0       ?        Original code, Mistra v7.3.2   <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE constants, ONLY :
     &     conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

      USE gas_common, ONLY: s1,s3

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      USE mod_nuc, ONLY:
     &     ind_H2SO4,
     &     ind_NH3,
     &     ind_OIO

      implicit none

! Common blocks:
      common /blck01/ am3(n),cm3(n)
      double precision am3,cm3

      common /nucl/ xn_new(n), xn_acc(n), xv_acc(n), dh2so4(n), dnh3(n)
      double precision xn_new, xn_acc, xv_acc, dh2so4, dnh3

      common /nuclio/ xn_newio(n), xn_accio(n), xv_accio(n), doio(n)
      double precision xn_newio, xn_accio, xv_accio, doio

      common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n),
     &                 bd_mean(n), dnucv(n), grorate(n), concnuc(n)
      double precision xn_app, xn_apacc, xv_apacc, bn_ges,
     &                 bd_mean, dnucv, grorate, concnuc

      double precision ssum
      common /nucsum/ ssum(n)
!- End of header ---------------------------------------------------------------

      write (21,106) xn_app(2),xn_apacc(2),bd_mean(2),bn_ges(2),
     &               grorate(2),concnuc(2)
      write (22,104) xn_new(2),xn_acc(2),xn_newio(2),xn_accio(2)
      write (23,106) (s1(ind_NH3,2)/am3(2)*1.d12),
     &               (s1(ind_H2SO4,2)*conv1),
     &               dnh3(2),dh2so4(2),
     &               (s3(ind_OIO,2)/am3(2)*1.d12),doio(2)
      write (24,102) (ssum(2)/am3(2)*1.d12),dnucv(2)

 102  format (2d16.8)
 104  format (4d16.8)
 106  format (6d16.8)

      end subroutine nucout1

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine nucout2
!
! Description:
!    output nucleation
!
! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      07/2016   Header, implicit none          <Josue Bock>
!                    Use module for parameters
!
! 1.1       ?                                       <Susanne Pechtl (b. Marquart)>
!
! 1.0       ?        Original code, Mistra v7.3.2   <Roland von Glasow>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit none

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday,lst,lmin,it,lcl,lct

      common /nucl/ xn_new(n), xn_acc(n), xv_acc(n), dh2so4(n), dnh3(n)
      double precision xn_new, xn_acc, xv_acc, dh2so4, dnh3

      common /nuclapp/ xn_app(n), xn_apacc(n), xv_apacc(n),bn_ges(n),
     &                 bd_mean(n), dnucv(n), grorate(n), concnuc(n)
      double precision xn_app, xn_apacc, xv_apacc, bn_ges,
     &                 bd_mean, dnucv, grorate, concnuc

      common /nuclio/ xn_newio(n), xn_accio(n), xv_accio(n), doio(n)
      double precision xn_newio, xn_accio, xv_accio, doio

! Local scalars:
      integer k
!- End of header ---------------------------------------------------------------


      write (25,104) lday,lst,lmin
      do k=1,n
         write (25,105) k,xn_new(k),xn_acc(k),xv_acc(k)
      enddo
      write (25,104) lday,lst,lmin
      do k=1,n
         write (25,105) k,xn_newio(k),xn_accio(k),xv_accio(k)
      enddo
      write (25,104) lday,lst,lmin
      do k=1,n
         write (25,105) k,xn_app(k),xn_apacc(k),xv_apacc(k)
      enddo

 104  format (3i4)
 105  format (i4,3d16.8)
  
      end subroutine nucout2
