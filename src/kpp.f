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


c additional subroutines for KPP version

      subroutine initc (box,n_bl)
c initialization of chemistry module

!      USE config, ONLY :
!     &     iod

      USE constants, ONLY :
! Imported Parameters:
     &     Avogadro,
     &     m_air

      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5,
! Imported Array Variables with intent(in):
     &     ind_gas,
     &     gas_is_halo,
     &     gas_name,
! Imported Array Variables with intent(inout):
     &     s1,
     &     s1_init_grd,
     &     s1_init_top,
     &     es1,
! Imported Array Variables with intent (out):
     &     s3

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j3,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc,
     &     nlev,
     &     nrxn,
     &     nmax_chem_aer

      USE kpp_aer_Parameters, ONLY :
     &     nspec_a=>NSPEC
      USE kpp_tot_Parameters, ONLY :
     &     nspec_t=>NSPEC

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      logical box
      integer n_bl

! Local scalars:
      double precision x0
      integer ia,j,k,kc
      logical tr_ue
      ! sea salt initialisation
      double precision xso42m,xhco3m,xno3m,xbrm,xclm,xim,xio3m,xiod

! Local arrays:
      double precision
     &     freep(nf),
     &     is4(n,3),
     &     xm(n),
     &     x4(n),
     &     x2(j1)

! Common blocks:
      common /blck01/ am3(n),cm3(n)
      double precision am3, cm3

      common /blck11/ rc(nkc,n)
      double precision rc

      common /blck12/ cw(nkc,n),cm(nkc,n)
      double precision cw, cm

      common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
      double precision conv2

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision :: sl1, sion1

      common /blck78/ sa1(nka,j2),sac1(nka,j2)
      double precision sa1,sac1

      common /budg/ bg(2,nrxn,nlev),il(nlev)
      double precision bg
      integer il

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      common /cb63/ fcs(nka),xmol3(nka)
      double precision fcs, xmol3

      common /kinv_i/ kinv
      integer kinv

      common /kpp_l1/ cloudt(nkc,n)
      logical cloudt

      common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss
      double precision xcryssulf,xcrysss,xdelisulf,xdeliss

      common /kpp_laer/ henry_la(NSPEC_a,nf),xkmt_la(nf,nkc,NSPEC_a),
     &     xkef_la(nf,nkc,NSPEC_a),xkeb_la(nf,nkc,NSPEC_a)
      double precision henry_la, xkmt_la, xkef_la, xkeb_la

      common /kpp_ltot/ henry_lt(NSPEC_t,nf),xkmt_lt(nf,nkc,NSPEC_t),
     &     xkef_lt(nf,nkc,NSPEC_t),xkeb_lt(nf,nkc,NSPEC_t)
      double precision henry_lt, xkmt_lt, xkef_lt, xkeb_lt
!- End of header ---------------------------------------------------------------


! print input concentrations and emission rates (user values):
! ------------------------------------------------------------
! initial mixing ratio of gas phase species in nmol mol-1 (=ppb)
      ! mixing ratio at ground
      write (60,6010) 
 6010 format (6x,'initial gas concentration at the surface [ppb]')
      write (60,6020) (s1_init_grd(j),j=1,j1)
      ! mixing ratio at top
      write (60,6012) 
 6012 format (6x,'initial gas concentration at the top [ppb]')
      write (60,6020) (s1_init_top(j),j=1,j1)
! emission rates of gas phase species in molecules/cm**2/s
      write (60,6025)
 6025 format (6x,'emission rates [molecules/cm**2/s]')
      write (60,6020) (es1(j),j=1,j1)
 6020 format (1x,10e12.5)


! jjb
! Initialisation to 0 of some arrays, which are not necessarily updated for layers nf+1 to n
         cw(:,:) = 0.d0
         rc(:,:) = 0.d0
         cm(:,:) = 0.d0
         conv2(:,:) = 0.d0

! jjb: partly from Peter Brauer version
! initiate all ions with marginal concentration to avoid computational problems
      sion1(:,:,:) = 0.d0
      sl1(:,:,:) = 0.d0

! jjb 14/02/2017
!     initialise the reaction rate arrays
      henry_la(:,:) = 0.d0
      henry_lt(:,:) = 0.d0
      xkmt_la(:,:,:) = 0.d0
      xkmt_lt(:,:,:) = 0.d0
      xkef_la(:,:,:) = 0.d0
      xkef_lt(:,:,:) = 0.d0
      xkeb_la(:,:,:) = 0.d0
      xkeb_lt(:,:,:) = 0.d0

c conversion of gaseous species and air density
c air density: [rho]=kg/m^3
      cm3(1)=rho(1)*Avogadro/m_air
      am3(1)=rho(1)/m_air

! Initialize arrays(k)
      do k=2,n
         cm3(k)=rho(k)*Avogadro/m_air        ! [air] in mlc/cm^3
         am3(k)=rho(k)/m_air                 ! [air] in mol/m^3
c conversion of gaseous species in ppb to mol/m**3(air)
         xm(k)=am3(k)*1.d-9                     ! ppb --> mol/m^3

c exp. decrease of concentrations from surface to free troposphere
         x4(k)=eta(k)/1900.
         x4(k)=min(1.d0,x4(k))
      end do

! Initialize arrays(j)
      do j=1,j1
         if(s1_init_grd(j).ne.0.d0) then
            ! avoid log(0) by adding a small value to s1_init_top
            x2(j) = -log(s1_init_grd(j))+log(s1_init_top(j)+1.d-10)
         else
            x2(j) = 0.d0
! The interpolation method does not allow to calculate a gradient
         ! (see below: s1(:,:) = s1_init_grd(:)*... )
         ! with a top concentration > 0 and a ground concentration = 0
         ! Warn the user if this case arise
            if(s1_init_top(j).ne.0.d0) then
               print*,"Warning with gas species nb. ",ind_gas(j)
               print*,"  Its top concentration is > 0 while its ground"
               print*,"  concentration is = 0"
               print*,"  The interpolation method does not allow this"
               print*,"  See SR initc"
            end if
         end if
      end do

! ......................................................................
! Initialise gas concentration in the whole column
      s1(:,1) = 0.d0
      do k=2,n
         do j=1,j1

            ! halogen only in BL, but there no gradient
            if(gas_is_halo(j) .and. k<kinv .and. k>2
     &         .and.trim(gas_name(j))/='HCl') then
               s1(j,k) = s1(j,k-1)
            else if(gas_is_halo(j) .and. k>=kinv
     &              .and. trim(gas_name(j))/='HCl' ) then
               s1(j,k) = 0.d0

            ! general case
            else
               s1(j,k)=s1_init_grd(j)*exp(x4(k)*x2(j))*xm(k)
            end if
         end do
      end do
! ......................................................................

c initial radical concentrations in mol/m**3(air)
      do k=1,n
         do j=1,j5
            s3(j,k)=0.
         end do
      end do
      xiod=0.
      !if (iod) xiod=1.


c define crystallization and deliquescene rel humidities (Seinfeld and Pandis, 
c Fig 9.4, p. 519)
c usually several cloud cycles have been made therefore sulfate aerosol is 
c usually shrinking; for sea salt crystallization point is used as well, because 
c sea salt particle were produced as droplets and shrank, but did not get "dry"
c if rel hum is below crys rH in FT, and it's getting more humid then the deli rH
c has to be taken for reactivation of aerosol chemistry
      xcryssulf= 0.4  ! crystallization humities 
      xcrysss  = 0.42  
      xdelisulf= 0.7  ! assumed based on mixing between different salts - should be
                      ! calculated explicitly if this starts to be critical (ie for
                      ! non-marine cases)
      xdeliss  = 0.75
      print *,'deliquescence rH  ',xdelisulf,xdeliss
      print *,'crystallization rH',xcryssulf,xcrysss
c array cloudt: was there a cloud in this layer in previous timestep?
c initialize with "true" to assure that aerosol chemistry is on in FT for
c layers in which rH > xcrystallization
      do k=1,n
         do kc=1,nkc
            cloudt(kc,k)=.true.
         enddo
      enddo
c initial loading of aerosols with nh3,so4,fe(3),mn(2) (x0=mole/particle) 
c watch out: sa1 is defined as sa1(..,j2,..) but addressed in j6 (=ion, sion1) terms
c            except for DOM which is in sl1 (therefore it is in j2)!!
      do ia=1,nka
cc ocean aerosol: particles with rn(ia)<.5 mum: 32% (NH4)2SO4, 64% NH4HSO4, 4% NH4NO3
c ocean aerosol: particles with rn(ia)<.5 mum: 34% (NH4)2SO4, 65.6% NH4HSO4, 0.4% NH4NO3
         if (rn(ia).lt.0.5) then  
            x0=en(ia)*1.d-03*fcs(ia)/xmol3(ia) 
c            sa1(ia,13)=x0*0.04 !NO3-
            sa1(ia,13)=x0*0.004 !NO3-
            sa1(ia,2)=x0*1.34   !NH4+
            sa1(ia,8)=x0*0.34   !SO4=    
            sa1(ia,19)=x0*0.656 !HSO4-    
c larger particles: pure nacl
         else
            x0=en(ia)*1.d-03*fcs(ia)/xmol3(ia) 
c sea salt particle
c x0 = mol / particle
c all the xiii are scaled to the sum of all negative ions in seawater, 
c Na+ is the sum of all positive ions; to get the correct molar ratios of
c Cl- or Br- to Na+, the lumped Na+ has to be multiplied by 0.806
            xso42m=0.0485
            xhco3m=4.2d-3
            xno3m =1.0d-7
            xbrm  =1.45d-3
            xim   =7.4d-8/.545*xiod
            xio3m =2.64d-7/.545*xiod
            xclm  =1-(xso42m+xhco3m+xno3m+xbrm+xim+xio3m)
            sa1(ia, 8)=xso42m*x0    ! SO4=    
            sa1(ia, 9)=xhco3m*x0    ! HCO3-
            sa1(ia,13)=xno3m*x0     ! NO3-
            sa1(ia,14)=xclm*x0      ! Cl-
            sa1(ia,20)=x0           ! "Na+" -  eletronegativity
            sa1(ia,24)=xbrm*x0      ! Br-
            sa1(ia,34)=xim*x0       ! I-
            sa1(ia,36)=xio3m*x0     ! IO3-
            sa1(ia,j2-j3+4)=0.27*xbrm*x0 !unspecified DOM
                                    !according to #2210: 0.27*[Br-]; enriched compared to ocean water ratio
         endif

cc urban aerosol:
cc 2/3*xm = mole mass nh4no3; 1/3*xm = mole mass (nh4)2so4;
cc xm=en(ia)*fcs(ia)*1d-3: total soluble mole mass
ccc --> 4/3*xm mole nh3, 2/3*xm mole no3, 1/3*xm mole so4
c            if (xmol3(ia).lt.130.) then
c               x0=en(ia)*1.d-03*fcs(ia)/(3.*xmol3(ia))
c               sa1(ia,3)=x0*2.
c               sa1(ia,4)=x0*4.
c               sa1(ia,6)=x0
c            else
          end do

c print initial concentrations (continued)
      write (60,6030)
 6030 format (6x,'sa1(nka,4)')
      write (60,6020) (sa1(ia,4),ia=1,nka)
      write (60,6040)
 6040 format (6x,'sa1(nka,6)')
      write (60,6020) (sa1(ia,6),ia=1,nka)
      write (60,6050)
c 6050 format (6x,'sa1(nka,j2-j3+4)')
c      write (60,6020) (sa1(ia,j2-j3+4),ia=1,nka)
c      write (60,6060)
c 6060 format (6x,'sa1(nka,j2-j3+5)')
c      write (60,6020) (sa1(ia,j2-j3+5),ia=1,nka)
 6050 format (6x,'sa1(nka,14)')
      write (60,6020) (sa1(ia,14),ia=1,nka)
      write (60,6060)
 6060 format (6x,'sa1(nka,24)')
      write (60,6020) (sa1(ia,24),ia=1,nka)

c levels for  rate output 
      il(1) =  5
      if (box) il(1)=n_bl
      il(2) = 15
      il(3) = 25
      il(4) = 35
      il(5) = 45
      il(6) = 55
      il(7) = 65
      il(8) = 75
      il(9) = 85
      il(10)= 95
      il(11)=105
      il(12)=115
      il(13)=125
      il(14)=135
      il(15)=145
c initial output for plotting
      do k=1,n
            is4(k,1)=am3(k) ! stay consistent with plot routine array size!
            is4(k,2)=0.
            is4(k,3)=0.
c            write (542,12) k,am3(k,1),am3(k,2)
c            write (543,12) k,cm3(k,1),cm3(k,2)
      enddo
! 12   format (i3,2d16.8)
      write (61) is4
      close (61)
      write (64) is4
      close (64)

! Initialise vmean (constant factor calculation)
      call v_mean_init
! ... and make the first calculation
      call v_mean (t(:nmax_chem_aer))

c Init calc of v_mean and henry
c initialize all variables also for box run
!      call v_mean_a  (t,nf)
!     call henry_a (t,p,nf) ! jjb second argument (p) not used
      call henry_a (t,nf)   ! jjb removed
      call st_coeff_a
!      call v_mean_t  (t,nf)
!     call henry_t (t,p,nf) ! jjb second argument (p) not used
      call henry_t (t,nf)   ! jjb removed
      call st_coeff_t

c initialize all levels also in box run
c free path length (lambda=freep):
      do k=1,nf
         freep(k)=2.28d-5 * t(k) / p(k)
      enddo
      tr_ue = .false. !.true.
      call init_konc
      call fast_k_mt_a(freep,tr_ue,nf)
      call activ_init

      end subroutine initc

c
c-----------------------------------------------------
c

      subroutine liq_parm (xra,box,n_bl)
c parameter for liquid phase chemistry
c LWC, henry, k_mt needed for calculation of gas <--> liquid transfer
c (see cb kpp_l)
c nkc=4; 1: sulfate aerosol, 2: seasalt aerosol, 3: sulfate droplets
c        4: seasalt droplets
c xra: aerodynamic resistence, needed for calculation of dry deposition velocities

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nkc,
     &     nmax_chem_aer

      implicit double precision (a-h,o-z)

      common /blck01/ am3(n),cm3(n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kpp_l1/ cloudt(nkc,n)
      logical cloudt

      logical box
      dimension freep(nf)
      logical update_now

      itime = int(time)

      nmin  = 1
      nmin2 = 2
      nmaxf = nf
      nmax  = n
      if (box) then
         nmin  = n_bl
         nmin2 = n_bl
         nmaxf = n_bl
         nmax  = n_bl
      endif

c free path length (lambda=freep):
      do k=nmin,nmaxf
         freep(k)=2.28d-5 * t(k) / p(k)
      enddo

c conversion of gaseous species and air density
c air density: [rho]=kg/m^3
      cm3(1)=rho(1)*6.022d20/29.
      am3(1)=rho(1)/29.d-3
      cm3(nmin2:nmax)=rho(nmin2:nmax)*6.022d20/29.        ! [air] in mlc/cm^3
      am3(nmin2:nmax)=rho(nmin2:nmax)/29.d-3              ! [air] in mol/m^3


c dry deposition velocities for gas phase
      call gasdrydep (xra,t,rho,freep)

      call cw_rc (nmaxf)

! Check if the tot mechanism will be called or not (xph3/4 = 1.) and if a new bin is activated somewhere (update_now)
      xph3=0.
      xph4=0.
      update_now = .false.

      do k=nmin,nmaxf
         do kc=1,nkc
            if (cm(kc,k).gt.0.) then
               if (.not.cloudt(kc,k)) update_now = .true.
               if (kc.eq.3) xph3=1.
               if (kc.eq.4) xph4=1.
            endif
         enddo
      enddo


      call v_mean (t(:nmax_chem_aer))

! Call all subroutines needed for aerosols (bins 1 & 2) (aer mechanism)
      call henry_a (t,nmaxf)
      call st_coeff_a
      call equil_co_a (t,nmaxf)
      if (itime/120*120.eq.itime .or. update_now )
     &   call fast_k_mt_a(freep,box,n_bl)

! If necessary, call all subroutines needed for droplets (bins 3 & 4) (tot mechanism)
      if (xph3.eq.1..or.xph4.eq.1.) then
         call henry_t (t,nmaxf)
         call st_coeff_t
         call equil_co_t (t,nmaxf)
         if(itime/120*120.eq.itime .or. update_now)
     &       call fast_k_mt_t(freep,box,n_bl)
      endif

c calculate rate for surface reaction OH + Cl-
!      call gamma_surf (box,n_bl) ! jjb not used
      
c calculate rates for DRY heterogeneous reactions
      call dry_cw_rc (nmax)
      call dry_rates_g (t,p,nmax)
      call dry_rates_a (freep,nmaxf)
      call dry_rates_t (freep,nmaxf)

      call activ (box,n_bl)

      end subroutine liq_parm


c
c--------------------------------------------------------------------------------
c

      subroutine st_coeff_t
c sticking coefficients, needed for calculation of k_mt

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n

      implicit none

      include 'tot_Parameters.h' !additional common blocks and other definitions

      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      double precision alpha, vmean

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      double precision, parameter :: RGAS= 8.314  ! gas constant [J/(mol*K)]
      double precision, parameter :: CAL = 4.1868 ! conversion calorie->Joule [J/cal]

      double precision, parameter :: C_o_R = CAL/RGAS

      double precision :: RT               ! RGAS * T(k)
      double precision :: C_o_RT           ! CAL / (RGAS*T(k))
      double precision :: tcorr

      integer :: j, k

! Initialisation (default value = 0.1)
      alpha(:,:) = 0.1d0

      do k=2,nf

         tcorr=1./t(k)-1./298.15
         RT = RGAS*t(k)
         C_o_RT = CAL/RT

c for some species a temperature dependence is calculated/estimated
c using:
c alpha = 1./(1.+1./(1./(1./alpha(T0)-1.)*EXP((-DeltaH/RGAS)*TCORR)))
c standard value
      alpha(ind_H2SO4,k) =0.65
      alpha(ind_CH4,k) =0.1
      alpha(ind_HNO4,k) =0.1
c experimentally determined values (from MOCCA)
      alpha(ind_O3P,k) =1.0d-6  ! JPL #15
      alpha(ind_O1D,k) =1.0d-6  ! JPL #15
      alpha(ind_O3,k) =  2.0D-03
      alpha(ind_O2,k) =
     &     1./(1.+1./(1./(1./1.0d-2 -1.)*exp(2000.d0*TCORR))) ! 06.04.00
c      alpha(ind_OH,k) =  4.0D-03 
      alpha(ind_OH,k) =  1.0D-2
      alpha(ind_HO2,k) =  2.0D-01 
c      alpha(ind_H2O2,k) =  9.1D-02 
      alpha(ind_H2O2,k) =
     &     1./(exp(-26.d3/RT+107.8456/RGAS) + 1.) 
      alpha(ind_NO,k) =  5.0D-05 
      alpha(ind_NO2,k) =  1.5D-03 
c      alpha(ind_NO3,k) =  2.5D-03 
      alpha(ind_NO3,k) =  4.0D-02 
c      alpha(ind_N2O5,k) =  2.7D-02 
      alpha(ind_N2O5,k) =  1.0D-01
      alpha(ind_HONO,k) =  4.0D-02 
c      alpha(ind_HNO3,k) =  8.6D-02 
      alpha(ind_HNO3,k) =  5.0D-01
      alpha(ind_NH3,k) =  6.0D-02 
      alpha(ind_MO2,k) =
     &     1./(1.+1./(1./(1./1.0d-2-1.)*exp(2000.d0*TCORR)))   ! 06.04.00
c+      alpha(ind_ROOH,k) =  5.5D-03 
c      alpha(ind_ROOH,k) =  5.5D-03 
c      alpha(ind_ROOH,k) =  0.01
      alpha(ind_ROOH,k)=
     &     1./(exp(-6.5D3*C_o_RT+32.5*C_o_R)+1.)
      alpha(ind_HCHO,k) =  4.0D-02 
c      alpha(ind_ACO2,k) =  1.8D-02  !HCOOH
      alpha(ind_ACO2,k)=
     &     1./(exp(-7.9E3*C_o_RT+34.9*C_o_R)+1.)   !HCOOH
      alpha(ind_ACTA,k) =  6.7D-02  ! at 273 K, Jayne et al., 1991
      alpha(ind_CH3OH,k) =  5.6D-02  ! at 273 K, Jayne et al., 1991
      alpha(ind_C2H5OH,k) =  4.8D-02  ! at 273 K, Jayne et al., 1991
      alpha(ind_CO2,k) =
     &     1./(1.+1./(1./(1./1.0d-2-1.)*exp(2000.d0*TCORR))) !06.04.00
c      alpha(ind_HCl,k) =  7.2D-02 
c      alpha(ind_HCl,k) = 0.1
      alpha(ind_HCl,k) =  1./(exp(-3.072d3/t(k) + 1.283d1)+1.) !T=290: 0.096, T=270: 0.190
c      alpha(ind_HOCl,k) =  7.2D-02 
c      alpha(ind_HOCl,k) =  5.0D-01 see below
      alpha(ind_ClNO3,k) =  1.0D-01 
c      alpha(ind_Cl2,k) =  5.5D-02 
      alpha(ind_Cl2,k) =
     &     1./(exp(-1.3D4*C_o_RT+50.*C_o_R)+1.)
c      alpha(ind_HBr,k) =  7.2D-02 
c      alpha(ind_HBr,k) = 0.05
      alpha(ind_HBr,k) = 1./(exp(-3.94d3/t(k) + 1.664d1) + 1.) !T=290K: 0.017, T=270K: 0.130
      alpha(ind_HOBr,k) =  6.0D-01  ! #1077
      alpha(ind_HOCl,k) = alpha(ind_HOBr,k)
      alpha(ind_BrNO3,k) =  8.0D-01 
c      alpha(ind_Br2,k) =  5.5D-02 
      alpha(ind_Br2,k) =
     &     1./(exp(-1.3D4*C_o_RT+50.*C_o_R)+1.)
c      alpha(ind_BrCl,k) =  5.5D-02 
      alpha(ind_BrCl,k) = 0.33 !#840 alpha(ind_Cl2,k)
      alpha(ind_SO2,k) =  1.1D-01 
c      alpha(ind_CH3SO3H,k) =  8.4D-02
      alpha(ind_CH3SO3H,k)=
     &     1./(exp(-3.50D3*C_o_RT+16.7*C_o_R)
     &     +1.) ! MSA, #955
      alpha(ind_DMS,k) = 1.0D-2 ! assumed
c      alpha(ind_DMSO,k) =  5.6D-02 
      alpha(ind_DMSO,k) =1./(exp(-5.12D3*C_o_RT+23.1*C_o_R)
     &     +1.)  ! #955
      alpha(ind_DMSO2,k) =
     &     1./(exp(-10.7D3*C_o_RT+43.0*C_o_R)
     &     +1.) ! #955
      alpha(ind_CH3SO2H,k) = 2.0D-4 ! assumed #2123, MSIA
c no uptake of other DMS products like CH3SCH2OO, CH3S, CH3SO, CH3SO2, CH3SO3
      alpha(ind_INO3,k) = 1./(1.+1./(1./(1./1.0d-1 -1.)*
     &                   exp(2000.d0*TCORR)))!06.04.00
cc      alpha(ind_HOI,k) =  7.2D-02 
cc      alpha(ind_HOI,k) =  5.0D-01 
      alpha(ind_HOI,k) = alpha(ind_HOBr,k) 
cc      alpha(ind_HI,k) =  7.2D-02 
      alpha(ind_HI,k) =  1./(exp(-4.13d3/t(k) + 1.715d1)+1.)
      alpha(ind_I2,k) =  1./(1.+1./(1./(1./ 1.0d-2 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_IO,k) =  1./(1.+1./(1./(1./ 5.0d-1 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_I2O2,k) = 1./(1.+1./(1./(1./ 1.0d-1 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_ICl,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_IBr,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_INO2,k) = 1./(1.+1./(1./(1./ 1.0d-1 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
c      alpha(ind_ClCHOk,k) =  1./(1.+1./(1./(1./ 1.0d-1 -1.)*
c     &                 exp(2000.d0*TCORR)))!06.04.00
c      alpha(ind_BrCHOk,k) =  1./(1.+1./(1./(1./ 1.0d-1 -1.)*
c     &                  exp(2000.d0*TCORR)))!06.04.00
c      alpha(ind_OIO,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*exp(2000.d0*
c     &     TCORR)))             ! assumed, #980
      alpha(ind_OIO,k) = 1.
      alpha(ind_HIO3,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*exp(2000.d0*
     &     TCORR)))             ! assumed, #980
      alpha(ind_XOR,k) = 7.0d-2  ! same as bromoethanol, Jayne et al., 1991
c     alpha(ind_I2O,k) = ?
c     alpha(ind_I2O3,k) = ?

c      alpha(ind_Hg,k)    =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgO,k)   =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgCl,k)  =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgCl2,k) =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgBr,k)  =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgBr2,k) =  1.d-1       ! caution - wild guess, no information found!!

      enddo

c check that alpha <= 1 to avoid that the T-dependencies mess up the
c  numbers
      do k=2,nf
         do j=1,NSPEC
            alpha(j,k)=min(1.d0,alpha(j,k))
         end do
      end do

      end subroutine st_coeff_t

c
c--------------------------------------------------------------------------------
c

      subroutine st_coeff_a
c sticking coefficients, needed for calculation of k_mt

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n

      implicit none

      include 'aer_Parameters.h' !additional common blocks and other definitions

      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      double precision alpha, vmean

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

      double precision, parameter :: RGAS= 8.314  ! gas constant [J/(mol*K)]
      double precision, parameter :: CAL = 4.1868 ! conversion calorie->Joule [J/cal]

      double precision, parameter :: C_o_R = CAL/RGAS

      double precision :: RT               ! RGAS * T(k)
      double precision :: C_o_RT           ! CAL / (RGAS*T(k))
      double precision :: tcorr

      integer :: j, k

! Initialisation (default value = 0.1)
      alpha(:,:) = 0.1d0

      do k=2,nf

         tcorr=1./t(k)-1./298.15
         RT = RGAS*t(k)
         C_o_RT = CAL/RT

c for some species a temperature dependence is calculated/estimated
c using:
c alpha = 1./(1.+1./(1./(1./alpha(T0)-1.)*EXP((-DeltaH/RGAS)*TCORR)))
c standard value
      alpha(ind_H2SO4,k) =0.65
      alpha(ind_CH4,k) =0.1
      alpha(ind_HNO4,k) =0.1
c experimentally determined values (from MOCCA)
      alpha(ind_O3P,k) =1.0d-6  ! JPL #15
      alpha(ind_O1D,k) =1.0d-6  ! JPL #15
      alpha(ind_O3,k) =  2.0D-03
      alpha(ind_O2,k) =
     &     1./(1.+1./(1./(1./1.0d-2 -1.)*exp(2000.d0*TCORR))) ! 06.04.00
c      alpha(ind_OH,k) =  4.0D-03 
      alpha(ind_OH,k) =  1.0D-2
      alpha(ind_HO2,k) =  2.0D-01 
c      alpha(ind_H2O2,k) =  9.1D-02 
      alpha(ind_H2O2,k) =
     &     1./(exp(-26.d3/RT+107.8456/RGAS) + 1.) 
      alpha(ind_NO,k) =  5.0D-05 
      alpha(ind_NO2,k) =  1.5D-03 
c      alpha(ind_NO3,k) =  2.5D-03 
      alpha(ind_NO3,k) =  4.0D-02 
c      alpha(ind_N2O5,k) =  2.7D-02 
      alpha(ind_N2O5,k) =  1.0D-01
      alpha(ind_HONO,k) =  4.0D-02 
c      alpha(ind_HNO3,k) =  8.6D-02 
      alpha(ind_HNO3,k) =  5.0D-01
      alpha(ind_NH3,k) =  6.0D-02 
      alpha(ind_MO2,k) =
     &     1./(1.+1./(1./(1./1.0d-2-1.)*exp(2000.d0*TCORR)))   ! 06.04.00
c+      alpha(ind_ROOH,k) =  5.5D-03 
c      alpha(ind_ROOH,k) =  5.5D-03 
c      alpha(ind_ROOH,k) =  0.01
      alpha(ind_ROOH,k)=
     &     1./(exp(-6.5D3*C_o_RT+32.5*C_o_R)+1.)
      alpha(ind_HCHO,k) =  4.0D-02 
c      alpha(ind_ACO2,k) =  1.8D-02  !HCOOH
      alpha(ind_ACO2,k)=
     &     1./(exp(-7.9E3*C_o_RT+34.9*C_o_R)+1.)   !HCOOH
      alpha(ind_ACTA,k) =  6.7D-02  ! at 273 K, Jayne et al., 1991
      alpha(ind_CH3OH,k) =  5.6D-02  ! at 273 K, Jayne et al., 1991
      alpha(ind_C2H5OH,k) =  4.8D-02  ! at 273 K, Jayne et al., 1991
      alpha(ind_CO2,k) =
     &     1./(1.+1./(1./(1./1.0d-2-1.)*exp(2000.d0*TCORR))) !06.04.00
c      alpha(ind_HCl,k) =  7.2D-02 
c      alpha(ind_HCl,k) = 0.1
      alpha(ind_HCl,k) =  1./(exp(-3.072d3/t(k) + 1.283d1)+1.) !T=290: 0.096, T=270: 0.190
c      alpha(ind_HOCl,k) =  7.2D-02 
c      alpha(ind_HOCl,k) =  5.0D-01 see below
      alpha(ind_ClNO3,k) =  1.0D-01 
c      alpha(ind_Cl2,k) =  5.5D-02 
      alpha(ind_Cl2,k) =
     &     1./(exp(-1.3D4*C_o_RT+50.*C_o_R)+1.)
c      alpha(ind_HBr,k) =  7.2D-02 
c      alpha(ind_HBr,k) = 0.05
      alpha(ind_HBr,k) = 1./(exp(-3.94d3/t(k) + 1.664d1) + 1.) !T=290K: 0.017, T=270K: 0.130
      alpha(ind_HOBr,k) =  6.0D-01  ! #1077
      alpha(ind_HOCl,k) = alpha(ind_HOBr,k)
      alpha(ind_BrNO3,k) =  8.0D-01 
c      alpha(ind_Br2,k) =  5.5D-02 
      alpha(ind_Br2,k) =
     &     1./(exp(-1.3D4*C_o_RT+50.*C_o_R)+1.)
c      alpha(ind_BrCl,k) =  5.5D-02 
      alpha(ind_BrCl,k) = 0.33 !#840 alpha(ind_Cl2,k)
      alpha(ind_SO2,k) =  1.1D-01 
c      alpha(ind_CH3SO3H,k) =  8.4D-02
      alpha(ind_CH3SO3H,k)=
     &     1./(exp(-3.50D3*C_o_RT+16.7*C_o_R)
     &     +1.) ! MSA, #955
      alpha(ind_DMS,k) = 1.0D-2 ! assumed
c      alpha(ind_DMSO,k) =  5.6D-02 
      alpha(ind_DMSO,k) =1./(exp(-5.12D3*C_o_RT+23.1*C_o_R)
     &     +1.)  ! #955
      alpha(ind_DMSO2,k) =
     &     1./(exp(-10.7D3*C_o_RT+43.0*C_o_R)
     &     +1.) ! #955
      alpha(ind_CH3SO2H,k) = 2.0D-4 ! assumed #2123, MSIA
c no uptake of other DMS products like CH3SCH2OO, CH3S, CH3SO, CH3SO2, CH3SO3
      alpha(ind_INO3,k) = 1./(1.+1./(1./(1./1.0d-1 -1.)*
     &                   exp(2000.d0*TCORR)))!06.04.00
cc      alpha(ind_HOI,k) =  7.2D-02 
cc      alpha(ind_HOI,k) =  5.0D-01 
      alpha(ind_HOI,k) = alpha(ind_HOBr,k) 
cc      alpha(ind_HI,k) =  7.2D-02 
      alpha(ind_HI,k) =  1./(exp(-4.13d3/t(k) + 1.715d1)+1.)
      alpha(ind_I2,k) =  1./(1.+1./(1./(1./ 1.0d-2 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_IO,k) =  1./(1.+1./(1./(1./ 5.0d-1 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_I2O2,k) = 1./(1.+1./(1./(1./ 1.0d-1 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_ICl,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_IBr,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
      alpha(ind_INO2,k) = 1./(1.+1./(1./(1./ 1.0d-1 -1.)*
     &                 exp(2000.d0*TCORR)))!06.04.00
c      alpha(ind_ClCHOk,k) =  1./(1.+1./(1./(1./ 1.0d-1 -1.)*
c     &                 exp(2000.d0*TCORR)))!06.04.00
c      alpha(ind_BrCHOk,k) =  1./(1.+1./(1./(1./ 1.0d-1 -1.)*
c     &                  exp(2000.d0*TCORR)))!06.04.00
c      alpha(ind_OIO,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*exp(2000.d0*
c     &     TCORR)))             ! assumed, #980
      alpha(ind_OIO,k) = 1.
      alpha(ind_HIO3,k) = 1./(1.+1./(1./(1./ 1.0d-2 -1.)*exp(2000.d0*
     &     TCORR)))             ! assumed, #980
      alpha(ind_XOR,k) = 7.0d-2  ! same as bromoethanol, Jayne et al., 1991
c     alpha(ind_I2O,k) = ?
c     alpha(ind_I2O3,k) = ?

c      alpha(ind_Hg,k)    =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgO,k)   =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgCl,k)  =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgCl2,k) =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgBr,k)  =  1.d-1       ! caution - wild guess, no information found!!
c      alpha(ind_HgBr2,k) =  1.d-1       ! caution - wild guess, no information found!!

      enddo

c check that alpha <= 1 to avoid that the T-dependencies mess up the
c  numbers
      do k=2,nf
         do j=1,NSPEC
            alpha(j,k)=min(1.d0,alpha(j,k))
         end do
      end do

      end subroutine st_coeff_a


c
c-----------------------------------------------------
c

      subroutine v_mean_init

      USE constants, ONLY :
! Imported Parameters:
     &     gas_const,
     &     pi

      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5,
     &     j4,
! Imported Array Variables with intent (in):
     &     gas_mass,
     &     rad_mass,
     &     fix_mass,
! Imported Array Variables with intent (out):
     &     vmean_init,
     &     vmean

      USE global_params, ONLY :
! Imported Parameters:
     &     nmax_chem_aer

      implicit none

      integer :: jtot
      integer :: jspec
      double precision :: const_fact

      double precision :: sqrt_mass (j1 + j5 + j4)

      jtot = j1 + j5 + j4

      allocate ( vmean_init(jtot) )
      allocate ( vmean(jtot,nmax_chem_aer) )

      const_fact = sqrt(8.d0*gas_const/pi)

      do jspec = 1,j1
         sqrt_mass(jspec) = sqrt(gas_mass(jspec))
      end do
      do jspec = 1,j5
         sqrt_mass(j1+jspec) = sqrt(rad_mass(jspec))
      end do
      do jspec = 1,j4
         sqrt_mass(j1+j5+jspec) = sqrt(fix_mass(jspec))
      end do

      do jspec = 1,jtot
         vmean_init(jspec) = const_fact / sqrt_mass(jspec)
      end do

      end subroutine v_mean_init

c
c-----------------------------------------------------
c

      subroutine v_mean (temperature)

      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5,
     &     j4,
! Imported Array Variables with intent (in):
     &     vmean_init,
! Imported Array Variables with intent (out):
     &     vmean

      USE global_params, ONLY :
! Imported Parameters:
     &     nmax_chem_aer

      implicit none

      double precision, intent(in) :: temperature (nmax_chem_aer)
      integer :: jtot
      integer :: j,k
      double precision :: sqrtt(nmax_chem_aer)

      jtot = j1 + j5 + j4

      sqrtt = sqrt(temperature)

      do k=1,nmax_chem_aer
         do j=1,jtot
            vmean(j,k) = vmean_init(j) * sqrtt(k)
         end do
      end do

      end subroutine v_mean
c
c-----------------------------------------------------
c

      subroutine henry_t (tt,nmaxf)

c Temp dependent Henry constants
c inverse dimensionless Henry constant: k_(H,inv)^cc:=1/(k_H^cp*RT)
c in equilibrium: XXXaq = k_H^cp * LWC * XXXg

! jjb work done:
!     - removed pp (pressure) from the argument list: unused
!     - removed hard coded parameters, use modules instead
!     - improved computation efficiency: Tfact calculated only once per layer
!     - missing declarations and implicit none
!     - final conversion (inverse if henry /= 0.) optimised
!     - reindexed henry(NSPEC,nf) instead of (nf,NSPEC) for computing efficiency
!     - added initialisation of henry!
!     - added test: nmaxf must be <= nf

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nkc

      implicit none

      include 'tot_Parameters.h'    !additional common blocks and other definitions   

      double precision, intent(in) :: tt(n) ! temperature array
      integer         , intent(in) :: nmaxf ! max layer index where henry has to be computed

      double precision :: func3, a0, b0 ! temperature dependency function, and its arguments
      double precision :: FCT           ! conversion factor, see below
      double precision :: Tfact         ! Tfact, local variable [K**-1]

      integer :: j,k                    ! loop indexes

      common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      double precision henry, xkmt, xkef, xkeb

c 0.082=8.3145*10^3/101325=R/p_0*10^3 : R includes conversion from M/atm --> mol/(m^3*Pa)
c so k_H^cp is taken as M/atm (the most common literature unit): SEE END OF SR

c      func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))*0.082*tt(k0)
!      func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))

! jjb 11/02/2017 slight improvement in calculation efficiency: compute Tfact only once per layer
!     Tfact = 1/T - 1/Tref, see below
      func3(a0,b0)=a0*exp(b0*Tfact)

      if (nmaxf > nf) stop 'Error in henry_t: nmaxf must be <= nf'

      henry(:,:) = 0.d0

      do k=1,nmaxf
         
c         henry(ind_H2SO4,k)=4.1d-12   !???
         henry(ind_H2SO4,k)=1.d+16    !NIST --> RS_Gmitro_Vermeulen
         henry(ind_CH4,k)=1.3d-3    !RS_Mackay
         henry(ind_C2H6,k)=2.0d-3   !RS_Mackay
c         henry(ind_C3H8,k)=1.4d-3   !RS_Mackay
c         henry(ind_ALKA,k)=9.5d-4   !ok
         henry(ind_ETHE,k)=4.9d-3   !ok
c         henry(ind_ALKE,k)=4.9d-3   !ok
c         henry(ind_AROM,k)=1.5d-1   !RS_Mackay: methylbenzene
c         henry(ind_KET,k)=10.d0     !RS_higher ketone
c         henry(ind_CRES,k)=8.2d2    !RS_Hine
c         henry(ind_DIAL,k)=10.0d0   !RS_Snider:propenal #?
c         henry(ind_GLYX,k)=3.6d5    !RS_Zhou
c        henry(ind_NH4NO3,k)=4.1d-12  !# was soll dieses Salz in Gasphase??
c         henry(ind_RAN2,k)=1.2d0    !RS_Kames: pentyl-nitrate
c         henry(ind_RAN1,k)=1.d0     !RS_hine
c         henry(ind_N2O5,k)=0.0      !RS_Sander  infinity!
c         henry(ind_ClNO2,k)=0.0     !  10^-2
c         henry(ind_ClNO3,k)=0.0     !  infinity
c         henry(ind_BrNO2,k)=0.0     !  10^-1
c         henry(ind_BrNO3,k)=0.0     !  infinity
         henry(ind_HI,k)=0.0d0        !   (dissociation)
         henry(ind_I2O2,k)=0.0d0      !
         henry(ind_INO2,k)=0.0d0      !
         henry(ind_INO3,k)=0.0d0      !
!         henry(ind_I2O,k)=0.0d0      !
!         henry(ind_I2O3,k)=0.0d0      !
c         henry(ind_OIO,k)  =   unknown but not needed as only surface reaction and no reversible uptake
c         henry(ind_HIO3,k) =                      -"-
         henry(ind_C3H7I,k)=1.1d-1  !RS_Hine
c        henry(ind_CH3SO2,k)= !RS_MOCCA
c        henry(ind_CH3SO3,k)= !RS_MOCCA
c         henry(ind_Hg,k)    = 1.3d-1                   ! #233 in #3127
c         henry(ind_HgO,k)   = 2.69d12                  ! #233 in #3127
c         henry(ind_HgCl,k)  = 2.75d6                   ! assumed HgCl2 (poss. higher??)
c         henry(ind_HgCl2,k) = 2.75d6                   ! #233 in #3127
c         henry(ind_HgBr,k)  = 2.75d6                   ! assumed HgBr2 (poss. higher??)
c         henry(ind_HgBr2,k) = 2.75d6                   ! #3127

c explicitly Temp dependent
         ! Tfact = 1/T - 1/Tref  with Tref = 298.15
         ! 1/298.15 = 3.3540d-3
         Tfact = 1.d0/tt(k) - 3.3540d-3
         
         henry(ind_NO,k)=func3(1.9d-03,1480.d0) !RS_Lide
c         henry(ind_NO2,k)=func3(1.d-02,2500.d0) !RS_Chameides
         henry(ind_NO2,k)=func3(6.4d-03,2500.d0) !RS_MOCCA
c         henry(ind_HNO3,k)=func3(1.66d5,8694.d0) !RS_MOCCA
         henry(ind_HNO3,k)=func3(2.5d6/1.5d1,8694.d0) !RS_MOCCA_exakt
c         henry(ind_HNO4,k)=1.4d4     !Goetz, 1996, cited in Warneck, 1999, #695
         henry(ind_HNO4,k)=func3(1.2d4,6900.d0) !06.04.00
         henry(ind_NH3,k)=func3(58.d0,4085.d0) !RS_MOCCA
         henry(ind_SO2,k)=func3(1.2d0,3120.d0) !RS_MOCCA
         henry(ind_O3,k)=func3(1.2d-02,2560.d0) !RS_MOCCA
         henry(ind_ACO2,k)=func3(3.7d+03,5700.d0) !RS_MOCCA
         henry(ind_ACTA,k)=func3(4.1d+03,6300.d0) !RS_Johnson
         henry(ind_HCHO,k)=func3(7.0d+03,6425.d0) !RS_MOCCA
         henry(ind_ALD2,k)=func3(1.3d+01,5700.d0) !RS_Benkelberg: acetaldehyde
         henry(ind_H2O2,k)=func3(1.d+05,6338.d0) !RS_MOCCA
c         henry(ind_ROOH,k)=func3(7.45d+04,6620.d0) ! # aktualisieren
         henry(ind_ROOH,k)=func3(3.0d+02,5322.d0) !RS_MOCCA
c         henry(ind_HONO,k)=func3(5.0d+01,4900.d0) !RS_Becker
         henry(ind_HONO,k)=func3(4.9d+01,4780.d0) !RS_MOCCA
         henry(ind_PAN,k)=func3(2.8d0,6500.d0) !RS_Kames
c         henry(ind_TPAN,k)=henry(22,k)
c         henry(ind_MGLY,k)=func3(3.7d+03,7553.d0) !RS_Betterton
c         henry(ind_HCl,k)=func3(1.17d0,9001.d0) !RS_MOCCA
         henry(ind_HCl,k)=func3(2.d0/1.7d0,9001.d0) !RS_MOCCA_exakt
c         henry(ind_R3N2,k)=func3(1.d0,5450.d0) !RS_kames: propyl nitrate
         henry(ind_NO3,k)=func3(2.d0,2000.d0) !RS_MOCCA
         henry(ind_DMS,k)=func3(4.8d-1,3100.d0) !RS_deBruyn 
c         henry(ind_DMSO,k)=5.d4     !RS_MOCCA
         henry(ind_DMSO,k)=func3(5.d4,6425.d0)     !RS_MOCCA 06.04.00
         henry(ind_DMSO2,k)= 1.d+16    !DMSO2=H2SO4, assumed
         henry(ind_CH3SO2H,k)=1.d+16     !MSIA=H2SO4, assumed
         henry(ind_CH3SO3H,k)=1.d+16    !MSA=H2SO4, assumed
         henry(ind_HOCl,k)=func3(6.7d2,5862.d0) !RS_MOCCA
c         henry(ind_Cl2,k)=9.2d-2    !RS_MOCCA
         henry(ind_Cl2,k)=func3(9.1d-2,2500.d0)    !RS_MOCCA 06.04.00
         henry(ind_HBr,k)=func3(1.3d0,10239.d0) !RS_MOCCA
         henry(ind_Br2,k)=func3(7.6d-1,4094.d0) !RS_MOCCA
         henry(ind_BrCl,k)=func3(9.4d-1,5600.d0)   !RS_MOCCA
c         henry(ind_HOBr,k)=9.3d1    !RS_MOCCA
         henry(ind_HOBr,k)=func3(9.3d1,5862.d0)    !RS_MOCCA 06.04.00
cc         henry(ind_HI,k)=func3(2.5d9/K_a,9800.d0) !RS_Brimblecombe #K_a
         henry(ind_I2,k)=func3(3.d0,4431.d0) !RS_MOCCA
cc         henry(ind_HOI,k)=4.5d2     !RS_MOCCA
         henry(ind_HOI,k)=func3(4.5d2,5862.d0)     !RS_MOCCA 06.04.00
cc         henry(ind_ICl,k)=1.1d2     !RS_MOCCA
         henry(ind_ICl,k)=func3(1.1d2,5600.d0)     !RS_MOCCA 06.04.00
cc         henry(ind_IBr,k)=2.4d1     !RS_MOCCA
         henry(ind_IBr,k)=func3(2.4d1,5600.d0)     !RS_MOCCA 06.04.00
         henry(ind_CH3I,k)=func3(1.4d-1,4300.d0) !RS_Moore
         henry(ind_CH2I2,k)=func3(2.3d0,5000.d0) !RS_Moore
         henry(ind_CH2ClI,k)=func3(8.9d-1,4300.d0) !RS_Moore

c for radicals only OH, HO2, MO2 are defined ! #
c         henry(ind_OH,k)=25.d0      !RS_MOCCA
         henry(ind_OH,k)=func3(3.0d1,4300.d0)      !RS_MOCCA 06.04.00
c         henry(ind_HO2,k)=9.d3      !RS_MOCCA
         henry(ind_HO2,k)=func3(3.9d3,5900.d0)      !RS_MOCCA 06.04.00
         henry(ind_MO2,k)=func3(6.d0,5600.d0) !RS_Jacob
c         henry(ind_MO2,k)=6.d0      !RS_MOCCA
cc         henry(ind_IO,k)=4.5d2      !RS_MOCCA
         henry(ind_IO,k)=func3(4.5d2,5862.d0)      !RS_MOCCA 06.04.00
         henry(ind_CO2,k)=func3(3.1d-02,2423.d0) !RS_MOCCA
         henry(ind_CO,k)=func3(9.9d-04,1300.d0) !RS_MOCCA
         henry(ind_O2,k)=func3(1.3d-3,1500.d0) !RS_Lide95 06.04.00
c         henry(ind_O2,k)=1.7d-3     !RS_MOCCA
c         henry(ind_I2O4,k)=func3()
c         henry(ind_I2O5,k)=func3()
c         henry(ind_INO,k)=func3()
c         henry(ind_Br2O,k)=func3()
         henry(ind_ClONO,k)=4.6d-2  !RS_MOCCA, same as ClNO2
c         henry(ind_ClO3,k)=func3()
c         henry(ind_Cl2O3,k)=func3()
         henry(ind_CH3OH,k)=func3(1.6d2,5600.d0)  !RS_MOCCA
         henry(ind_C2H5OH,k)=func3(1.5d2,6400.d0)  !RS_MOCCA
         henry(ind_H2,k)=func3(7.8d-4,500.d0)  !RS_MOCCA
c         henry(ind_NHS,k)=func3()
c         henry(ind_RCl,k)=func3()
c         henry(ind_RBr,k)=func3()
         henry(ind_XOR,k)=func3(1.5d2,6400.d0) ! same as ethanol
         henry(ind_SOR,k)=func3(1.5d2,6400.d0) !  same as ethanol
c         henry(ind_SPAN,k)=func3()


c include HERE call to get effective henry coeffs (if needed), 
c before inverse k_H are calculated
c         xCO2=func3(4.3d-7,-919.d0)
c         x3CO2=xCO2*func3(4.7d-11,-1787.d0)
c         xNH3=func3(1.71d-5,-4325.d0)/func3(1.d-14,-6716.d0)
c         if (sion1(1,1,k).gt.0.) then         !anhaengig von tropfenklasse! (--> pH)
c            xhp=(dmax(sion1(1,1,k),1.d-6))*1.d-3                   ! in mol/l !
c            henry(ind_CO2,k)=henry(ind_CO2,k)*(1+xCO2/xhp+x3CO2/(xhp**2))
c            henry(ind_NH3,k)=henry(ind_NH3,k)*(1+xNH3*xhp)
c         endif

      enddo

c unit: mol/(l*atm) --> mol(aq)/m3(aq) / mol(g)/m3(g) (i.e. dimensionless)
c FCT=1.d3*8.3145*T/p_0=0.082*T
c PLUS conversion to inverse Henry constant:  h_(H,inv)^cc = 1/k_H^cc
c i.e. mol(aq)/m3(aq) / mol(g)/m3(air) -->  mol(g)/m3(air) / mol(aq)/m3(aq)

      do k=1,nmaxf
         FCT=0.0820577d0*tt(k)
         do j=1,NSPEC
            if (henry(j,k).ne.0.d0) then
               henry(j,k)=1.d0/(henry(j,k)*FCT)
! "else": henry=0 <=> k_H^cc=infinity
            end if
         end do
      end do

      end subroutine henry_t


c
c------------------------------------------------------
c

      subroutine henry_a (tt,nmaxf)

c Temp dependent Henry constants
c inverse dimensionless Henry constant: k_(H,inv)^cc:=1/(k_H^cp*RT)
c in equilibrium: XXXaq = k_H^cp * LWC * XXXg

! jjb work done:
!     - removed pp (pressure) from the argument list: unused
!     - removed hard coded parameters, use modules instead
!     - improved computation efficiency: Tfact calculated only once per layer
!     - missing declarations and implicit none
!     - final conversion (inverse if henry /= 0.) optimised
!     - reindexed henry(NSPEC,nf) instead of (nf,NSPEC) for computing efficiency
!     - added initialisation of henry!
!     - added test: nmaxf must be <= nf

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nkc

      implicit none

      include 'aer_Parameters.h'    !additional common blocks and other definitions   

      double precision, intent(in) :: tt(n) ! temperature array
      integer         , intent(in) :: nmaxf ! max layer index where henry has to be computed

      double precision :: func3, a0, b0 ! temperature dependency function, and its arguments
      double precision :: FCT           ! conversion factor, see below
      double precision :: Tfact         ! Tfact, local variable [K**-1]

      integer :: j,k                    ! loop indexes

      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      double precision henry, xkmt, xkef, xkeb

c 0.082=8.3145*10^3/101325=R/p_0*10^3 : R includes conversion from M/atm --> mol/(m^3*Pa)
c so k_H^cp is taken as M/atm (the most common literature unit): SEE END OF SR

c      func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))*0.082*tt(k0)
!      func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))

! jjb 11/02/2017 slight improvement in calculation efficiency: compute Tfact only once per layer
!     Tfact = 1/T - 1/Tref, see below
      func3(a0,b0)=a0*exp(b0*Tfact)

      if (nmaxf > nf) stop 'Error in henry_t: nmaxf must be <= nf'

      henry(:,:) = 0.d0

      do k=1,nmaxf

c         henry(ind_H2SO4,k)=4.1d-12   !???
         henry(ind_H2SO4,k)=1.d+16    !NIST --> RS_Gmitro_Vermeulen
         henry(ind_CH4,k)=1.3d-3    !RS_Mackay
         henry(ind_C2H6,k)=2.0d-3   !RS_Mackay
c         henry(ind_C3H8,k)=1.4d-3   !RS_Mackay
c         henry(ind_ALKA,k)=9.5d-4   !ok
         henry(ind_ETHE,k)=4.9d-3   !ok
c         henry(ind_ALKE,k)=4.9d-3   !ok
c         henry(ind_AROM,k)=1.5d-1   !RS_Mackay: methylbenzene
c         henry(ind_KET,k)=10.d0     !RS_higher ketone
c         henry(ind_CRES,k)=8.2d2    !RS_Hine
c         henry(ind_DIAL,k)=10.0d0   !RS_Snider:propenal #?
c         henry(ind_GLYX,k)=3.6d5    !RS_Zhou
c        henry(ind_NH4NO3,k)=4.1d-12  !# was soll dieses Salz in Gasphase??
c         henry(ind_RAN2,k)=1.2d0    !RS_Kames: pentyl-nitrate
c         henry(ind_RAN1,k)=1.d0     !RS_hine
c         henry(ind_N2O5,k)=0.0      !RS_Sander  infinity!
c         henry(ind_ClNO2,k)=0.0     !  10^-2
c         henry(ind_ClNO3,k)=0.0     !  infinity
c         henry(ind_BrNO2,k)=0.0     !  10^-1
c         henry(ind_BrNO3,k)=0.0     !  infinity
         henry(ind_HI,k)=0.0d0        !   (dissociation)
         henry(ind_I2O2,k)=0.0d0      !
         henry(ind_INO2,k)=0.0d0      !
         henry(ind_INO3,k)=0.0d0      !
!         henry(ind_I2O,k)=0.0d0      !
!         henry(ind_I2O3,k)=0.0d0      !
c         henry(ind_OIO,k)  =   unknown but not needed as only surface reaction and no reversible uptake
c         henry(ind_HIO3,k) =                      -"-
         henry(ind_C3H7I,k)=1.1d-1  !RS_Hine
c        henry(ind_CH3SO2,k)= !RS_MOCCA
c        henry(ind_CH3SO3,k)= !RS_MOCCA
c         henry(ind_Hg,k)    = 1.3d-1                   ! #233 in #3127
c         henry(ind_HgO,k)   = 2.69d12                  ! #233 in #3127
c         henry(ind_HgCl,k)  = 2.75d6                   ! assumed HgCl2 (poss. higher??)
c         henry(ind_HgCl2,k) = 2.75d6                   ! #233 in #3127
c         henry(ind_HgBr,k)  = 2.75d6                   ! assumed HgBr2 (poss. higher??)
c         henry(ind_HgBr2,k) = 2.75d6                   ! #3127

c explicitly Temp dependent
         ! Tfact = 1/T - 1/Tref  with Tref = 298.15 K
         ! 1/298.15 = 3.3540d-3
         Tfact = 1.d0/tt(k) - 3.3540d-3
         
         henry(ind_NO,k)=func3(1.9d-03,1480.d0) !RS_Lide
c         henry(ind_NO2,k)=func3(1.d-02,2500.d0) !RS_Chameides
         henry(ind_NO2,k)=func3(6.4d-03,2500.d0) !RS_MOCCA
c         henry(ind_HNO3,k)=func3(1.66d5,8694.d0) !RS_MOCCA
         henry(ind_HNO3,k)=func3(2.5d6/1.5d1,8694.d0) !RS_MOCCA_exakt
c         henry(ind_HNO4,k)=1.4d4     !Goetz, 1996, cited in Warneck, 1999, #695
         henry(ind_HNO4,k)=func3(1.2d4,6900.d0) !06.04.00
         henry(ind_NH3,k)=func3(58.d0,4085.d0) !RS_MOCCA
         henry(ind_SO2,k)=func3(1.2d0,3120.d0) !RS_MOCCA
         henry(ind_O3,k)=func3(1.2d-02,2560.d0) !RS_MOCCA
         henry(ind_ACO2,k)=func3(3.7d+03,5700.d0) !RS_MOCCA
         henry(ind_ACTA,k)=func3(4.1d+03,6300.d0) !RS_Johnson
         henry(ind_HCHO,k)=func3(7.0d+03,6425.d0) !RS_MOCCA
         henry(ind_ALD2,k)=func3(1.3d+01,5700.d0) !RS_Benkelberg: acetaldehyde
         henry(ind_H2O2,k)=func3(1.d+05,6338.d0) !RS_MOCCA
c         henry(ind_ROOH,k)=func3(7.45d+04,6620.d0) ! # aktualisieren
         henry(ind_ROOH,k)=func3(3.0d+02,5322.d0) !RS_MOCCA
c         henry(ind_HONO,k)=func3(5.0d+01,4900.d0) !RS_Becker
         henry(ind_HONO,k)=func3(4.9d+01,4780.d0) !RS_MOCCA
         henry(ind_PAN,k)=func3(2.8d0,6500.d0) !RS_Kames
c         henry(ind_TPAN,k)=henry(22,k)
c         henry(ind_MGLY,k)=func3(3.7d+03,7553.d0) !RS_Betterton
c         henry(ind_HCl,k)=func3(1.17d0,9001.d0) !RS_MOCCA
         henry(ind_HCl,k)=func3(2.d0/1.7d0,9001.d0) !RS_MOCCA_exakt
c         henry(ind_R3N2,k)=func3(1.d0,5450.d0) !RS_kames: propyl nitrate
         henry(ind_NO3,k)=func3(2.d0,2000.d0) !RS_MOCCA
         henry(ind_DMS,k)=func3(4.8d-1,3100.d0) !RS_deBruyn 
c         henry(ind_DMSO,k)=5.d4     !RS_MOCCA
         henry(ind_DMSO,k)=func3(5.d4,6425.d0)     !RS_MOCCA 06.04.00
         henry(ind_DMSO2,k)= 1.d+16    !DMSO2=H2SO4, assumed
         henry(ind_CH3SO2H,k)=1.d+16     !MSIA=H2SO4, assumed
         henry(ind_CH3SO3H,k)=1.d+16    !MSA=H2SO4, assumed
         henry(ind_HOCl,k)=func3(6.7d2,5862.d0) !RS_MOCCA
c         henry(ind_Cl2,k)=9.2d-2    !RS_MOCCA
         henry(ind_Cl2,k)=func3(9.1d-2,2500.d0)    !RS_MOCCA 06.04.00
         henry(ind_HBr,k)=func3(1.3d0,10239.d0) !RS_MOCCA
         henry(ind_Br2,k)=func3(7.6d-1,4094.d0) !RS_MOCCA
         henry(ind_BrCl,k)=func3(9.4d-1,5600.d0)   !RS_MOCCA
c         henry(ind_HOBr,k)=9.3d1    !RS_MOCCA
         henry(ind_HOBr,k)=func3(9.3d1,5862.d0)    !RS_MOCCA 06.04.00
cc         henry(ind_HI,k)=func3(2.5d9/K_a,9800.d0) !RS_Brimblecombe #K_a
         henry(ind_I2,k)=func3(3.d0,4431.d0) !RS_MOCCA
cc         henry(ind_HOI,k)=4.5d2     !RS_MOCCA
         henry(ind_HOI,k)=func3(4.5d2,5862.d0)     !RS_MOCCA 06.04.00
cc         henry(ind_ICl,k)=1.1d2     !RS_MOCCA
         henry(ind_ICl,k)=func3(1.1d2,5600.d0)     !RS_MOCCA 06.04.00
cc         henry(ind_IBr,k)=2.4d1     !RS_MOCCA
         henry(ind_IBr,k)=func3(2.4d1,5600.d0)     !RS_MOCCA 06.04.00
         henry(ind_CH3I,k)=func3(1.4d-1,4300.d0) !RS_Moore
         henry(ind_CH2I2,k)=func3(2.3d0,5000.d0) !RS_Moore
         henry(ind_CH2ClI,k)=func3(8.9d-1,4300.d0) !RS_Moore

c for radicals only OH, HO2, MO2 are defined ! #
c         henry(ind_OH,k)=25.d0      !RS_MOCCA
         henry(ind_OH,k)=func3(3.0d1,4300.d0)      !RS_MOCCA 06.04.00
c         henry(ind_HO2,k)=9.d3      !RS_MOCCA
         henry(ind_HO2,k)=func3(3.9d3,5900.d0)      !RS_MOCCA 06.04.00
         henry(ind_MO2,k)=func3(6.d0,5600.d0) !RS_Jacob
c         henry(ind_MO2,k)=6.d0      !RS_MOCCA
cc         henry(ind_IO,k)=4.5d2      !RS_MOCCA
         henry(ind_IO,k)=func3(4.5d2,5862.d0)      !RS_MOCCA 06.04.00
         henry(ind_CO2,k)=func3(3.1d-02,2423.d0) !RS_MOCCA
         henry(ind_CO,k)=func3(9.9d-04,1300.d0) !RS_MOCCA
         henry(ind_O2,k)=func3(1.3d-3,1500.d0) !RS_Lide95 06.04.00
c         henry(ind_O2,k)=1.7d-3     !RS_MOCCA
c         henry(ind_I2O4,k)=func3()
c         henry(ind_I2O5,k)=func3()
c         henry(ind_INO,k)=func3()
c         henry(ind_Br2O,k)=func3()
         henry(ind_ClONO,k)=4.6d-2  !RS_MOCCA, same as ClNO2
c         henry(ind_ClO3,k)=func3()
c         henry(ind_Cl2O3,k)=func3()
         henry(ind_CH3OH,k)=func3(1.6d2,5600.d0)  !RS_MOCCA
         henry(ind_C2H5OH,k)=func3(1.5d2,6400.d0)  !RS_MOCCA
         henry(ind_H2,k)=func3(7.8d-4,500.d0)  !RS_MOCCA
c         henry(ind_NHS,k)=func3()
c         henry(ind_RCl,k)=func3()
c         henry(ind_RBr,k)=func3()
         henry(ind_XOR,k)=func3(1.5d2,6400.d0) ! same as ethanol
         henry(ind_SOR,k)=func3(1.5d2,6400.d0) !  same as ethanol
c         henry(ind_SPAN,k)=func3()

         
c include HERE call to get effective henry coeffs (if needed), 
c before inverse k_H are calculated
c         xCO2=func3(4.3d-7,-919.d0)
c         x3CO2=xCO2*func3(4.7d-11,-1787.d0)
c         xNH3=func3(1.71d-5,-4325.d0)/func3(1.d-14,-6716.d0)
c         if (sion1(1,1,k).gt.0.) then         !anhaengig von tropfenklasse! (--> pH)
c            xhp=(dmax(sion1(1,1,k),1.d-6))*1.d-3                   ! in mol/l !
c            henry(ind_CO2,k)=henry(ind_CO2,k)*(1+xCO2/xhp+x3CO2/(xhp**2))
c            henry(ind_NH3,k)=henry(ind_NH3,k)*(1+xNH3*xhp)
c         endif

      enddo

c unit: mol/(l*atm) --> mol(aq)/m3(aq) / mol(g)/m3(g) (i.e. dimensionless)
c FCT=1.d3*8.3145*T/p_0=0.082*T
c PLUS conversion to inverse Henry constant:  h_(H,inv)^cc = 1/k_H^cc
c i.e. mol(aq)/m3(aq) / mol(g)/m3(air) -->  mol(g)/m3(air) / mol(aq)/m3(aq)

      do k=1,nmaxf
         FCT=0.0820577d0*tt(k)
         do j=1,NSPEC
            if (henry(j,k).ne.0.d0) then
               henry(j,k)=1.d0/(henry(j,k)*FCT)
! "else": henry=0 <=> k_H^cc=infinity
            end if
         end do
      end do

!      write (1717,*)henry(:,10)

      end subroutine henry_a



c
c------------------------------------------------------
c

      subroutine cw_rc (nmaxf)
c mean radius and LWC for "chemical" particles size bins

! jjb work done
!     cleaning
!     modules, instead of hard coded parameters
!     rc has to be always defined
!     useless test removed to set cm (case rH>rH_deliq : whatever cloud value is ok to define cm)
!     inconsistency between .lt. and .gt. : case "==" was thus missing. corrected.
!     implicit none, and all missing declarations
!     ial definition only once, out of the main do loop

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka,
     &     nkt,
     &     nkc

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: nmaxf ! max index for calculation  ! jjb might add a test on it, error if > n, warning if > nf

! Local parameters:
      double precision, parameter :: xpi = 4./3.*pi

      double precision, parameter :: cwm = 1.d-1 ! Threshold for switching chemistry in bins 1 & 2 ("aerosols")
      double precision, parameter :: cwmd = 1.d2 ! Threshold for switching chemistry in bins 3 & 4 ("droplets") 
                                                 !   here "d" stands for droplets, not dry !

! Local scalars:
      double precision :: cm1, cm2, cm3, cm4
      double precision :: cw1, cw2, cw3, cw4
      double precision :: rc1, rc2, rc3, rc4
      double precision :: x0, x1
      integer :: ia, ial, jt ! loop indexes for 2D particle grid
      integer :: k  ! index for grid number in do loops

! Common blocks:
      common /blck06/ kw(nka),ka
      integer kw, ka

      common /blck11/ rc(nkc,n)
      double precision rc

      common /blck12/ cw(nkc,n),cm(nkc,n)
      double precision cw, cm

      common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
      double precision conv2

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), ! only rq and e are used
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw, ew, rn, rw, en, e, dew, rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n) ! only ff is used
      double precision ff, fsum
      integer nar

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n) ! only feu is used
      double precision xm1, xm2, feu, dfddt, xm1a, xm2a

      common /kpp_l1/ cloud(nkc,n)
      logical cloud

      common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss
      double precision xcryssulf,xcrysss,xdelisulf,xdeliss

      common /kinv_i/ kinv
      integer kinv

      common /nucfeed/ ifeed
      integer ifeed

c rc(nkc,n): mean radius of droplets in m
c cw(nkc,n): LWC in given radius range in m^3(aq)g/m^3(air) 
c      cwm=1.d-05

! Define lower bound depending on nucleation settings
      if (ifeed.eq.2) then
         ial = 2
      else
         ial = 1
      endif

! Main loop
      do k=2,nmaxf

c         if (k.lt.lcl.or.k.gt.lct) go to 1000
c         if (feu(k).lt.xcryssulf.and.feu(k).lt.xcrysss) go to 1000
c see SR kon: if rH>0.7 microphysics is calculated, this
c is also start for liquid aerosol chemistry

! Initialisation
         rc1=0. ; rc2=0. ; rc3=0. ; rc4=0.
         cw1=0. ; cw2=0. ; cw3=0. ; cw4=0.
         cm1=0. ; cm2=0. ; cm3=0. ; cm4=0.

c here TOTAL particle volume is used for calculating LWC
c this is correct only for completely soluble aerosol
c
c the cw variables and therefore also cvv?/conv2 use
c the volume of solution (cvv? converts to mol/l_solution)
c for the calculation of the activity coefficients the
c molality is needed (mol/kg_solvent) ==> cm is only water
c volume (m^3 water/m^3 air)


c small aerosol
         do ia=ial,ka
            do jt=1,kw(ia)
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cw1=cw1+x0
               rc1=rc1+x0*rq(jt,ia)
               x1=ff(jt,ia,k)*e(jt)
               cm1=cm1+x1
            enddo
c small droplets
            do jt=kw(ia)+1,nkt
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cw3=cw3+x0
               rc3=rc3+x0*rq(jt,ia)
               x1=ff(jt,ia,k)*e(jt)
               cm3=cm3+x1
            enddo
         enddo
c large aerosol
         do ia=ka+1,nka
            do jt=1,kw(ia)
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cw2=cw2+x0
               rc2=rc2+x0*rq(jt,ia)
               x1=ff(jt,ia,k)*e(jt)
               cm2=cm2+x1
            enddo
c large droplets
            do jt=kw(ia)+1,nkt
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cw4=cw4+x0
               rc4=rc4+x0*rq(jt,ia)
               x1=ff(jt,ia,k)*e(jt)
               cm4=cm4+x1
            enddo
         enddo

c conversion: um^3/cm^3 --> m^3(aq)/m^3(air):10^-12
c           : um        --> m               :10^-6
c aerosol: cw must be greater than 1.d-13  => cwm
c droplet: cw must be greater than 1.d-10  => cwmd
c conversion for cm: mg(water)/cm^3(air)
c                       --> m^3(wat)/m^3(air):10^-3

! define rc for all rH because rc is used in calculation of particle
! chemistry sedimentation
         if (cw1.gt.0.) then
            rc(1,k)=rc1/cw1*1.d-6
         else
            rc(1,k) = 0.d0
         end if
         if (cw2.gt.0.) then
            rc(2,k)=rc2/cw2*1.d-6
         else
            rc(2,k) = 0.d0
         end if
         if (cw3.gt.0.) then
            rc(3,k)=rc3/cw3*1.d-6
         else
            rc(3,k) = 0.d0
         end if
         if (cw4.gt.0.) then
            rc(4,k)=rc4/cw4*1.d-6
         else
            rc(4,k) = 0.d0
         end if

         cw(1,k)=cw1*1.d-12
         cw(2,k)=cw2*1.d-12
         cw(3,k)=cw3*1.d-12
         cw(4,k)=cw4*1.d-12


! cm (old:cw) is used as switch for aerosol chemistry therefore:
!     off: below crys point
!     on : only if above threshold value and crys rH (if cloud was true) or
!          deli rH (if cloud was false)

         if (feu(k).lt.min(xcryssulf,xcrysss)) then
            if (k.le.kinv) print*,k,feu(k),' below both crystal. points'
            cm(:,k)    = 0.d0
            conv2(:,k) = 0.d0

         else

            if ( (cw1.ge.cwm) .and.                       ! cw1 has to be > cwm, and...
     &          ((cloud(1,k).and.feu(k).ge.xcryssulf).or. !  ... if cloud is true, rH >= rH_crys is enough
     &           (feu(k).ge.xdelisulf)) ) then            !  ... else (if cloud is false), rH >= rH_deliq is necessary
                                                          !         (cloud can be true as well in the later case)
               cm(1,k)=cm1*1.d-3
               conv2(1,k) = 1.d9/cw1 ! 1.d-3/cw(1,k)
            else
               cm(1,k)=0.d0
               conv2(1,k) = 0.d0
            endif

            if ( (cw2.ge.cwm) .and.
     &          ((cloud(2,k).and.feu(k).ge.xcrysss).or.    ! same comments
     &           (feu(k).ge.xdeliss)) ) then 
               cm(2,k)=cm2*1.d-3
               conv2(2,k) = 1.d9/cw2
            else
               cm(2,k)=0.d0
               conv2(2,k) = 0.d0
            endif

            if (cw3.ge.cwmd) then
               cm(3,k)=cm3*1.d-3
               conv2(3,k) = 1.d9/cw3
            else
               cm(3,k)=0.d0
               conv2(3,k) = 0.d0
            endif

            if (cw4.ge.cwmd) then
               cm(4,k)=cm4*1.d-3
               conv2(4,k) = 1.d9/cw4
            else
               cm(4,k)=0.d0
               conv2(4,k) = 0.d0
            endif

         end if

      end do

      end subroutine cw_rc


c
c--------------------------------------------------------
c

      subroutine fast_k_mt_t (freep,box,n_bl)  !_1D
c transfer coefficient after Schwarz, 1986 (see Sander & Crutzen '96, JGR, 9127)
c but no mean values used (like in SR k_mt_a/t) but integrated values

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc


      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5,
! Imported Array Variables with intent (in):
     &     gas_k2m_t,
     &     rad_k2m_t,
     &     vm=>vmean


      USE kpp_tot_Global, ONLY :
     &     SPC_NAMES


      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h' !additional common blocks and other definitions
      parameter (nx=50)
      logical box
      common /blck06/ kw(nka),ka
      common /blck12/ cw(nkc,n),cm(nkc,n)

      common /cb40/ time,lday,lst,lmin,it,lcl,lct ! jjb only time is used, for potential error message
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)

      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
      common /nucfeed/ ifeed
!     dimension cw(nf,nkc),freep(nf),dndlogr(nkt),rqm(nkt,nka),lex(nx) ! jjb dndlogr not used
!     dimension cw(nf,nkc),freep(nf),rqm(nkt,nka),lex(nx)              ! jjb thus removed
      dimension freep(nf),rqm(nkt,nka),lex(nx)              ! jjb cw removed from argument list

      data lex /ind_NO2,ind_HNO3,ind_NH3,ind_SO2,ind_H2SO4,ind_O3,
     &    ind_ACO2,ind_HCHO,ind_H2O2,ind_HONO,ind_HCl,ind_N2O5,ind_HNO4,
     &    ind_NO3,ind_OH,ind_HO2,ind_MO2,ind_CO2,ind_O2,ind_ROOH,
     &    ind_HOCl,ind_Cl2,ind_HBr,ind_HOBr,ind_Br2,ind_BrCl,ind_DMSO,
     &    ind_ClNO3,ind_BrNO3,ind_CH3SO3H,ind_DMS,ind_CH3SO2H,ind_DMSO2,
     &    ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr,ind_OIO,ind_INO2,
     &    ind_INO3,ind_HI,ind_I2O2,ind_HIO3,ind_NO,ind_ACTA,ind_CH3OH,
     &    ind_C2H5OH,ind_XOR,ind_SOR/
c    &    1,2,3,4,5/           !ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr/
c ind_Hg,ind_HgO,ind_HgCl,ind_HgCl2,ind_HgBr,ind_HgBr2; also change/check nx

c change rq in um to rqm in m
!     do ia=1,nka
!        do jt=1,nkt
!           rqm(jt,ia)=rq(jt,ia)*1.d-6  ! jjb matrix below
!        enddo
!     enddo
      rqm(:,:) = rq(:,:) * 1.d-6        ! jjb matrix

      nmin=2
      nmax=nf
      if (box) then
         nmin=n_bl
         nmax=n_bl
      endif
c loop over vertical grid
      do 1000 k=nmin,nmax
c loop over the nkc different chemical bins
         do kc=1,nkc
            if (cm(kc,k).eq.0.) goto 1001  ! switch changed from cw
c loop over the species to be exchanged between gas and aqueous phase---
            do l=1,nx

! jjb search for desired species index, from lex table, in gas_m2k_t table
               jspec = 1
               ! search in non radical gas list
               do while(gas_k2m_t(jspec) /= lex(l) .and. jspec+1 <= j1)
                  jspec = jspec+1
               end do

               ! search in radical gas list
               if(gas_k2m_t(jspec) /= lex(l)) then
               jspec = 1
               do while(rad_k2m_t(jspec) /= lex(l) .and. jspec+1 <= j5)
                  jspec = jspec+1
               end do

               if(rad_k2m_t(jspec) /= lex(l)) then
                  if(trim(SPC_NAMES(lex(l))) == 'O2') then
                     jspec = j1+j5+1
                  else
                     if(time<121.d0 .and. k==nmin.and.kc==1) then
                        print*,"in fast_k_mt_t, error"
                        print*,spc_names(lex(l))," kpp index ",lex(l)
                     end if
                     cycle
                     !stop "stopped by SR fast_k_mt_t"
                  end if
               else
                  jspec = jspec + j1 ! add offset in vmean array (radical case)
               end if

               end if
! end jjb

c define summation limits (1) ---
               if (kc.eq.1.or.kc.eq.3) then
                  if (ifeed.eq.2) then
                    iia_0=2
                  else
                    iia_0=1
                  endif
                  iia_e=ka
               endif
               if (kc.eq.2.or.kc.eq.4) then
                  iia_0=ka+1
                  iia_e=nka
               endif

c fast version without logarithmic integration
! Initialisation of local variables
               x1=0.
               xk1=0.
               if (l.eq.1) xx1=0.
               if (alpha(lex(l),k).gt.0.) x1=4./(3.*alpha(lex(l),k))

               do ia=iia_0,iia_e
c define summation limits (2)
                  if (kc.eq.1.or.kc.eq.2) then
                     jjt_0=1
                     jjt_e=kw(ia)
                  endif
                  if (kc.eq.3.or.kc.eq.4) then
                     jjt_0=kw(ia)+1
                     jjt_e=nkt
                  endif

                  do jt=jjt_0,jjt_e
c conversion: um      --> m               :10^-6
                     rqq=rqm(jt,ia)
c kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
c     with D_g=lambda*vmean/3.
c here a volume weighted value is calculated, therefore weighting with r^3:
c kmt=4/3*pi/L*sum(a)*sum(r){r^3*N*kt}

! < jjb 21-12-2016
!                     x2=vmean(lex(l),k)/(rqq/freep(k)+x1) ![1/s]
                     x2=vm(jspec,k)/(rqq/freep(k)+x1) ![1/s]
! jjb 21-12-2016 >

c conversion: 1/cm^3 --> 1/m^3(air):10^6
                     xk1=xk1+x2*rqq*rqq*ff(jt,ia,k)*1.d6
c LWC weighted sedimentation velocity
                     if (l.eq.1) then
                        xvs=vterm(rqq,t(k),p(k))
                        xx1=xx1+rqq*rqq*rqq*xvs*ff(jt,ia,k)*1.d6
                     endif
                  enddo !jt
               enddo !ia
c k_mt=4*pi/(3*LWC)*sum
               if (cw(kc,k).gt.0.d0) then
!                  xkmt(k,kc,lex(l))=4.*3.1415927/(3.*cw(kc,k))*xk1 ![1/s]
                  xkmt(k,kc,lex(l))=4.*pi/(3.*cw(kc,k))*xk1 ![1/s]
c sedimentation velocity:
!                  if (l.eq.1) vt(kc,k)=4.*3.1415927/(3.*cw(kc,k))*xx1
                  if (l.eq.1) vt(kc,k)=4.*pi/(3.*cw(kc,k))*xx1
               end if
            enddo !l 

            goto 1002

c this part only for calculation of sedimentation velocity in "dry" layers
 1001       continue
            xx1=0.
c define summation limits (1) ---
            if (kc.eq.1.or.kc.eq.3) then
               if (ifeed.eq.2) then
                 iia_0=2
               else    
                 iia_0=1 
               endif
               iia_e=ka
            endif
            if (kc.eq.2.or.kc.eq.4) then
               iia_0=ka+1
               iia_e=nka
            endif
            do ia=iia_0,iia_e
c define summation limits (2)
               if (kc.eq.1.or.kc.eq.2) then
                  jjt_0=1
                  jjt_e=kw(ia)
               endif
               if (kc.eq.3.or.kc.eq.4) then
                  jjt_0=kw(ia)+1
                  jjt_e=nkt
               endif
               do jt=jjt_0,jjt_e
c conversion: um      --> m               :10^-6
                  rqq=rqm(jt,ia)
                  xvs=vterm(rqq,t(k),p(k))
                  xx1=xx1+rqq*rqq*rqq*xvs*ff(jt,ia,k)*1.d6
               enddo            !jt
            enddo               !ia
c sedimentation velocity:
!            if (cw(kc,k).gt.0.) vt(kc,k)=4.*3.1415927/(3.*cw(kc,k))*xx1
            if (cw(kc,k).gt.0.) vt(kc,k)=4.*pi/(3.*cw(kc,k))*xx1

 1002       continue
         enddo                  ! kc
         
 1000 continue  !k

      end subroutine fast_k_mt_t


c
c--------------------------------------------------------
c

!     subroutine fast_k_mt_a (freep,cw,box,n_bl)  !_1D ! jjb cw now in /blck12/
      subroutine fast_k_mt_a (freep,box,n_bl)  !_1D
c transfer coefficient after Schwarz, 1986 (see Sander & Crutzen '96, JGR, 9127)
c but no mean values used (like in SR k_mt_a/t) but integrated values

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc


      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
     &     j5,
! Imported Array Variables with intent (in):
     &     gas_k2m_a,
     &     rad_k2m_a,
     &     vm=>vmean

      USE kpp_aer_Global, ONLY :
     &     SPC_NAMES


 
      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions
!      parameter (nf=100,n=nf+50,nka=70,nkt=70,nkc=4,nx=50)
      parameter (nx=50)
      logical box
      common /blck06/ kw(nka),ka
      common /blck12/ cw(nkc,n),cm(nkc,n)

      common /cb40/ time,lday,lst,lmin,it,lcl,lct ! jjb only time is used, for potential error message
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated
      common /nucfeed/ ifeed
!     dimension cw(nf,nkc),freep(nf),dndlogr(nkt),rqm(nkt,nka),lex(nx) ! jjb dndlogr not used
!     dimension cw(nf,nkc),freep(nf),rqm(nkt,nka),lex(nx)              ! jjb thus removed
      dimension freep(nf),rqm(nkt,nka),lex(nx)              ! jjb cw removed from argument list

      data lex /ind_NO2,ind_HNO3,ind_NH3,ind_SO2,ind_H2SO4,ind_O3,
     &    ind_ACO2,ind_HCHO,ind_H2O2,ind_HONO,ind_HCl,ind_N2O5,ind_HNO4,
     &    ind_NO3,ind_OH,ind_HO2,ind_MO2,ind_CO2,ind_O2,ind_ROOH,
     &    ind_HOCl,ind_Cl2,ind_HBr,ind_HOBr,ind_Br2,ind_BrCl,ind_DMSO,
     &    ind_ClNO3,ind_BrNO3,ind_CH3SO3H,ind_DMS,ind_CH3SO2H,ind_DMSO2,
     &    ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr,ind_OIO,ind_INO2,
     &    ind_INO3,ind_HI,ind_I2O2,ind_HIO3,ind_NO,ind_ACTA,ind_CH3OH,
     &    ind_C2H5OH,ind_XOR,ind_SOR/
c    &    1,2,3,4,5/           !ind_HOI,ind_IO,ind_I2,ind_ICl,ind_IBr/
c ind_Hg,ind_HgO,ind_HgCl,ind_HgCl2,ind_HgBr,ind_HgBr2; also change/check nx

c change rq in um to rqm in m
!     do ia=1,nka
!        do jt=1,nkt
!           rqm(jt,ia)=rq(jt,ia)*1.d-6  ! jjb matrix below
!        enddo
!     enddo
      rqm(:,:) = rq(:,:) * 1.d-6        ! jjb matrix

      nmin=2
      nmax=nf
      if (box) then
         nmin=n_bl
         nmax=n_bl
      endif

c loop over vertical grid
      do 1000 k=nmin,nmax
c loop over the nkc different chemical bins
         do kc=1,2!nkc
            if (cm(kc,k).eq.0.) goto 1001 ! switch changed from cw
c loop over the species to be exchanged between gas and aqueous phase---
            do l=1,nx

! jjb search for desired species index, from lex table, in gas_m2k_a table
               !print*,spc_names(lex(l)),j1,j5
               jspec = 1
               ! search in non radical gas list
               do while(gas_k2m_a(jspec) /= lex(l) .and. jspec+1 <= j1)
                  jspec = jspec+1
               end do

               ! search in radical gas list
               if(gas_k2m_a(jspec) /= lex(l)) then
               jspec = 1
               do while(rad_k2m_a(jspec) /= lex(l) .and. jspec+1 <= j5)
                  jspec = jspec+1
               end do

               if(rad_k2m_a(jspec) /= lex(l)) then
                  
                  if(trim(SPC_NAMES(lex(l))) == 'O2') then
                     jspec = j1+j5+1
                  else
                     if(time<121.d0 .and. k==nmin.and.kc==1) then
                        print*,"in fast_k_mt_t, error"
                        print*,spc_names(lex(l))," kpp index ",lex(l)
                     end if
                     cycle
                     !stop "stopped by SR fast_k_mt_a"
                  end if
               else
                  jspec = jspec + j1 ! add offset in vmean array (radical case)
               end if

               end if
               !print*,"final jspec =",jspec
! end jjb

c define summation limits (1) ---
               if (kc.eq.1) then
                  if (ifeed.eq.2) then
                    iia_0=2
                  else    
                    iia_0=1 
                  endif
                  iia_e=ka
               endif
               if (kc.eq.2) then
                  iia_0=ka+1
                  iia_e=nka
               endif
c fast version without logarithmic integration
               x1=0.
               xk1=0.
               if (l.eq.1) xx1=0.
               if (alpha(lex(l),k).gt.0.) x1=4./(3.*alpha(lex(l),k))
               do ia=iia_0,iia_e
c define summation limits (2)
                  if (kc.eq.1) then
                     jjt_0=1
                     jjt_e=kw(ia)
                  endif
                  if (kc.eq.2) then
                     jjt_0=1
                     jjt_e=kw(ia)
                  endif
                  do jt=jjt_0,jjt_e
c conversion: um      --> m               :10^-6
                     rqq=rqm(jt,ia)
c kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
c     with D_g=lambda*vmean/3.
c here a volume weighted value is calculated, therefore weighting with r^3:
c kmt=4/3*pi/L*sum(a)*sum(r){r^3*N*kt}

! < jjb 21-12-2016
!                     x2=vmean(lex(l),k)/(rqq/freep(k)+x1) ![1/s]
                     x2=vm(jspec,k)/(rqq/freep(k)+x1) ![1/s]
! jjb 21-12-2016 >

c conversion: 1/cm^3 --> 1/m^3(air):10^6
                     xk1=xk1+x2*rqq*rqq*ff(jt,ia,k)*1.d6
c LWC weighted sedimentation velocity
                     if (l.eq.1) then
                        xvs=vterm(rqq,t(k),p(k))
                        xx1=xx1+rqq*rqq*rqq*xvs*ff(jt,ia,k)*1.d6
                     endif
                  enddo !jt
               enddo !ia
c k_mt=4*pi/(3*LWC)*sum
               if (cw(kc,k).gt.0.d0) then
!                  xkmt(k,kc,lex(l))=4.*3.1415927/(3.*cw(kc,k))*xk1 ![1/s]
                  xkmt(k,kc,lex(l))=4.*pi/(3.*cw(kc,k))*xk1 ![1/s]
c sedimentation velocity:
!                  if (l.eq.1) vt(kc,k)=4.*3.1415927/(3.*cw(kc,k))*xx1
                  if (l.eq.1) vt(kc,k)=4.*pi/(3.*cw(kc,k))*xx1
               end if
            enddo !l 

            goto 1002

c this part only for calculation of sedimentation velocity in "dry" layers
 1001       continue
            xx1=0.
c define summation limits (1) ---
            if (kc.eq.1) then
               if (ifeed.eq.2) then
                 iia_0=2
               else    
                 iia_0=1 
               endif
               iia_e=ka
            endif
            if (kc.eq.2) then
               iia_0=ka+1
               iia_e=nka
            endif
            do ia=iia_0,iia_e
c define summation limits (2)
               if (kc.eq.1) then
                  jjt_0=1
                  jjt_e=kw(ia)
               endif
               if (kc.eq.2) then
                  jjt_0=1
                  jjt_e=kw(ia)
               endif

               do jt=jjt_0,jjt_e
c conversion: um      --> m               :10^-6
                  rqq=rqm(jt,ia)
                  xvs=vterm(rqq,t(k),p(k))
                  xx1=xx1+rqq*rqq*rqq*xvs*ff(jt,ia,k)*1.d6
               enddo            !jt
            enddo               !ia
c sedimentation velocity:
!            if (cw(kc,k).gt.0.) vt(kc,k)=4.*3.1415927/(3.*cw(kc,k))*xx1
            if (cw(kc,k).gt.0.) vt(kc,k)=4.*pi/(3.*cw(kc,k))*xx1

 1002       continue
         enddo                  ! kc
         
 1000 continue                  !k


c      do k=2,nf
c            write (534, 1101) k,xkmt(k,1,ind_HNO3),1./xkmt(k,1,ind_HNO3),
c     &        xkmt(k,2,ind_HNO3),1./xkmt(k,2,ind_HNO3),
c     &        xkmt(k,1,ind_HNO3)*cw(1,k),xkmt(k,2,ind_HNO3)*cw(2,k)
c      enddo
! 1101 format(i4, 6d16.8)

      end subroutine fast_k_mt_a


c
c------------------------------------------------------
c

!     subroutine equil_co_t (cw,tt,nmaxf) ! jjb cw now in /blck12/ ! conv2 used directly
      subroutine equil_co_t (tt,nmaxf)
c equilibrium constant (see MOCCA)
c xkef: forward
c xkeb: backward
c activity coefficients included via xgamma 
c acidity constants Ka = XXXaf/XXXab [mol/l];
c and other equilibrium constants (f=forward, b=backward reaction);
c converted to [mol/m3(air)];
c absolute values are chosen arbitrarily to ensure fast equilibration;

      USE global_params, ONLY :
! Imported Parameters:
     &     j6,
     &     nf,
     &     n,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h' !additional common blocks and other definitions          
!      parameter (j6=55,nf=100,n=nf+50,nkc=4)

      common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
      common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated
      common /kpp_mol/ xgamma(nf,j6,nkc) 
!     dimension cw(nf,nkc),tt(n) ! jjb cw removed from argument list
      dimension tt(n)

      funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))
      do k=2,nmaxf
         do kc=1,nkc
c            if (cw(k,kc).gt.1.d-9) then
c            if (cw(k,kc).gt.0.) then  !cm is switch now
!            if (cm(kc,k).gt.0.) then
!               cv2=1.d-3/cw(kc,k)
!            else
!               cv2=0.
!            endif
            cv2 = conv2(kc,k)
            xsw=1.
!            if (cv2.eq.0.) xsw=0. ! jjb no need to compute all the xkef & xkeb which will be 0. !
            if (cv2.gt.0.) then ! but still better than initialising all NSPEC list !

c absolute values of back- and forward reactions are choosen arbitrarily to 
c ensure a quick (but numerically stable) equilibrium
c if illogical error messages arise, try fiddling with the absolute values
c of equilibrium reactions
            xkef(k,kc,ind_H2O) = xsw*funa(1.0d-5,-6716.d0,k) !schneller wird schlechter!
            xkeb(k,kc,ind_H2O) =
     &           1.0D9* cv2*xgamma(k,1,kc)*xgamma(k,3,kc)

c            xkef(k,kc,ind_HO2) = xsw*2.d5  !Chameides '84
            xkef(k,kc,ind_HO2) = xsw*1.6d5  !Weinstein-loyd and Schwartz, '91
            xkeb(k,kc,ind_HO2) =
     &           1.d10*cv2*xgamma(k,1,kc)*xgamma(k,11,kc)

            xkef(k,kc,ind_ACO2) = xsw*1.8D0
            xkeb(k,kc,ind_ACO2) =
     &           1.0D4*cv2*xgamma(k,1,kc)*xgamma(k,16,kc)

c            xkef(k,kc,ind_CO2) = xsw*funa(4.3d3,-913.d0,k)
c            xkeb(k,kc,ind_CO2) = 1.0D10* cv2   
            xkef(k,kc,ind_CO2) = xsw*funa(4.3d-2,-913.d0,k)
            xkeb(k,kc,ind_CO2) = 1.0D5*cv2*xgamma(k,1,kc)*xgamma(k,9,kc)   

            xkef(k,kc,ind_HONO) = xsw*funa(5.1d+3,-1260.d0,k)
            xkeb(k,kc,ind_HONO) =
     &           1.0D7*cv2*xgamma(k,1,kc)*xgamma(k,12,kc)

            xkef(k,kc,ind_HNO3) = xsw*funa(1.54d+10,8700.d0,k)
            xkeb(k,kc,ind_HNO3) =
     &           1.0D9*cv2*xgamma(k,1,kc)*xgamma(k,13,kc)

            xkef(k,kc,ind_HNO4) = xsw*2.0D3
            xkeb(k,kc,ind_HNO4) = 2.0D8 * cv2

            xkef(k,kc,ind_NH3) = xsw*funa(1.7d5,-4325.d0,k) 
            xkeb(k,kc,ind_NH3) =
     &           1.0D10* cv2*xgamma(k,3,kc)*xgamma(k,2,kc)

            xkef(k,kc,ind_HSO3ml1) = xsw*funa(6.0d2,1120.d0,k) 
     &           *xgamma(k,5,kc)
            xkeb(k,kc,ind_HSO3ml1) = 1.0D10* cv2*xgamma(k,1,kc)
     &           *xgamma(k,6,kc)

            xkef(k,kc,ind_H2SO4) = xsw*1.0d12 !Seinfeld, Pandis (1998), p.391
            xkeb(k,kc,ind_H2SO4)=
     &           1.0d9* cv2*xgamma(k,1,kc)*xgamma(k,19,kc)

            xkef(k,kc,ind_HSO4ml1) = xsw*funa(1.02d+6,2720.d0,k)
     &           *xgamma(k,19,kc)
            xkeb(k,kc,ind_HSO4ml1)=
     &           1.0D8*cv2*xgamma(k,1,kc)*xgamma(k,8,kc)

            xkef(k,kc,ind_SO2) = xsw*funa(1.7d8,2090.d0,k)
            xkeb(k,kc,ind_SO2) =
     &           1.0D10*cv2*xgamma(k,1,kc)*xgamma(k,5,kc)

!            xkef(k,kc,ind_HCHO) = 1.d10*cv2  !Chameides '84   !mechanism changed
!            xkeb(k,kc,ind_HCHO) = xsw*1.d5                    !mechanism changed

            xkef(k,kc,ind_HCl) = xsw*funa(1.7d10,6896.d0,k)
            xkeb(k,kc,ind_HCl) =
     &           1.0D4*cv2*xgamma(k,1,kc)*xgamma(k,14,kc)

            xkef(k,kc,ind_Cl2ml1) = xsw*5.2d4*xgamma(k,15,kc)    !Chameides '84
            xkeb(k,kc,ind_Cl2ml1) =  1.d10*cv2*xgamma(k,14,kc)

c            xkef(k,kc,ind_Cl2) = 1.1D-3   
c            xkeb(k,kc,ind_Cl2) = 2.1D2* cv2 

            xkef(k,kc,ind_HOCl) = xsw*3.2D2   
            xkeb(k,kc,ind_HOCl) =
     &           1.0D10*cv2*xgamma(k,1,kc)*xgamma(k,22,kc)

            xkef(k,kc,ind_HBr) = xsw*1.0D13   
            xkeb(k,kc,ind_HBr) =
     &           1.0D4* cv2*xgamma(k,1,kc)*xgamma(k,24,kc)

            xkef(k,kc,ind_Br2) = xsw*funa(2.95d4,-4068.d0,k)* !Liu et al, 2002, #2109
     &           xgamma(k,25,kc)
            xkeb(k,kc,ind_Br2) = funa(1.17d10,-1812.d0,k)* cv2 *
     &           xgamma(k,24,kc)

            xkef(k,kc,ind_HOBr) = xsw*funa(2.3d1,-3091.d0,k)
            xkeb(k,kc,ind_HOBr) =
     &           1.0D10*cv2*xgamma(k,1,kc)*xgamma(k,26,kc)

            xkef(k,kc,ind_BrCl2ml1) = funa(5.d9,1143.d0,k) * cv2  
     &           *xgamma(k,14,kc) !#894
            xkeb(k,kc,ind_BrCl2ml1) = xsw*1.3D9   *xgamma(k,28,kc)
c no activities used, they "cancel out"
            xkef(k,kc,ind_Br2Clml1) = 5.d9 * cv2!*xgamma(k,24,kc) !if original equilibrium 
            xkeb(k,kc,ind_Br2Clml1) = xsw*2.8D5!*xgamma(k,29,kc)  !       -"-  
            xkef(k,kc,ind_Br2l1) = 5.d9 * cv2!*xgamma(k,14,kc)    !if original equilibrium
            xkeb(k,kc,ind_Br2l1) = xsw*3.85D9!*xgamma(k,29,kc)    !       -"- 
            xkef(k,kc,ind_ICl) = 1.0D11  * cv2 *xgamma(k,14,kc)
            xkeb(k,kc,ind_ICl) = xsw*1.3D9   *xgamma(k,37,kc)
            xkef(k,kc,ind_IBr) = 1.0D11  * cv2  *xgamma(k,24,kc)
            xkeb(k,kc,ind_IBr) = xsw*3.5D8  *xgamma(k,38,kc)
c new (speculative) ICl <--> IBr equilibria, assumed to yield the same
c ICl/IBr ratio as in the BrCl <--> Br2 eqilibria (BrCl/Br2) 
c no activities used, they are not known!! 11.04.01 
           xkef(k,kc,ind_IClBrml1) = 5.d9 * cv2!*xgamma(k,24,kc)
           xkeb(k,kc,ind_IClBrml1) = xsw*2.8D5!*xgamma(k,29,kc)
           xkef(k,kc,ind_I2) = 5.d9 * cv2!*xgamma(k,14,kc)    
           xkeb(k,kc,ind_I2) = xsw*3.85D9!*xgamma(k,29,kc)
!            xkef(k,kc,ind_HIO2) = xsw*2.0D3
!            xkeb(k,kc,ind_HIO2) = 2.0D9 * cv2
            xkef(k,kc,ind_HIO3) = xsw*1.57D4
            xkeb(k,kc,ind_HIO3) = 1.0D5 * cv2

c mercury; no activity coefficients; forward not faster than 1.d14
c
c            xkef(k,kc,ind_HgOHpl1)    = 4.27d14 * cv2           ! #493,  K_eq = 4.27d10
c            xkeb(k,kc,ind_HgOHpl1)    = xsw * 1.d4
c            xkef(k,kc,ind_HgOH2l1)    = 2.6d14 * cv2            ! #4174, K_eq = 2.6d11
c            xkeb(k,kc,ind_HgOH2l1)    = xsw * 1.d3
c            xkef(k,kc,ind_HgSO3l1)    = 5.01d14 * cv2           ! #493,  K_eq = 5.01d12
c            xkeb(k,kc,ind_HgSO3l1)    = xsw * 1.d2
c            xkef(k,kc,ind_HgSO322ml1) = 2.5d14 * cv2            ! #4158, K_eq = 2.5d11
c            xkeb(k,kc,ind_HgSO322ml1) = xsw * 1.d3
c            xkef(k,kc,ind_HgOHCll1)   = 2.69d14 * cv2            ! #4158, K_eq = 2.69d7
c            xkeb(k,kc,ind_HgOHCll1)   = xsw * 1.d7
c            xkef(k,kc,ind_HgClpl1)    = 2.0d14 * cv2             ! #4174, K_eq = 2.0d7
c            xkeb(k,kc,ind_HgClpl1)    = xsw * 1.d7
c            xkef(k,kc,ind_HgCl2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d6
c            xkeb(k,kc,ind_HgCl2l1)    = xsw * 1.d8
c            xkef(k,kc,ind_HgCl3ml1)   = 6.7d8 * cv2             ! #4174, K_eq = 6.7d0
c            xkeb(k,kc,ind_HgCl3ml1)   = xsw * 1.d8
c            xkef(k,kc,ind_HgCl42ml1)  = 1.3d9 * cv2             ! #4174, K_eq = 1.3d1
c            xkeb(k,kc,ind_HgCl42ml1)  = xsw * 1.d8 
c            xkef(k,kc,ind_HgOHBrl1)   = 2.69d14 * cv2            ! #4158, assumed
c            xkeb(k,kc,ind_HgOHBrl1)   = xsw * 1.d7
c            xkef(k,kc,ind_HgBrpl1)    = 1.1d14 * cv2             ! #4174, K_eq = 1.1d9
c            xkeb(k,kc,ind_HgBrpl1)    = xsw * 1.d5
c            xkef(k,kc,ind_HgBr2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d8
c            xkeb(k,kc,ind_HgBr2l1)    = xsw * 1.d6
c            xkef(k,kc,ind_HgBr3ml1)   = 1.5d10 * cv2             ! #4174, K_eq = 1.5d2
c            xkeb(k,kc,ind_HgBr3ml1)   = xsw * 1.d8
c            xkef(k,kc,ind_HgBr42ml1)  = 2.3d10 * cv2             ! #4174, K_eq = 2.3d1
c            xkeb(k,kc,ind_HgBr42ml1)  = xsw * 1.d9

            else ! conv2(kc,k) = 0
               xkef(k,kc,:) = 0. ! jjb xsw was set to 0. in this case, leading to xkef = 0.
               xkeb(k,kc,:) = 0. ! jjb cv2 was tested equal to 0 in this case, leading to xkeb = 0.
            end if

         enddo
      enddo

      end subroutine equil_co_t


c
c------------------------------------------------------
c

!     subroutine equil_co_a (cw,tt,nmaxf) ! jjb cw now in /blck12/ ! conv2 used directly
      subroutine equil_co_a (tt,nmaxf)
c equilibrium constant (see MOCCA)
c xkef: forward
c xkeb: backward
c activity coefficients included via xgamma
c acidity constants Ka = XXXaf/XXXab [mol/l];
c and other equilibrium constants (f=forward, b=backward reaction);
c converted to [mol/m3(air)];
c absolute values are chosen arbitrarily to ensure fast equilibration;

      USE global_params, ONLY :
! Imported Parameters:
     &     j6,
     &     nf,
     &     n,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions          
!      parameter (j6=55,nf=100,n=nf+50,nkc=4)

      common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated
      common /kpp_mol/ xgamma(nf,j6,nkc) 
!     dimension cw(nf,nkc),tt(n) ! jjb cw removed from argument list
      dimension tt(n)

      funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))
      do k=2,nmaxf
         do kc=1,2 !nkc
c            if (cw(k,kc).gt.1.d-9) then
c            if (cw(k,kc).gt.0.) then !cm is switch now
!            if (cm(kc,k).gt.0.) then
!               cv2=1.d-3/cw(kc,k)
!            else
!               cv2=0.
!            endif
            cv2 = conv2(kc,k)
            xsw=1.
!           if (cv2.eq.0.) xsw=0. ! jjb no need to compute all the xkef & xkeb which will be 0. !
            if (cv2.gt.0.) then ! but still better than initialising all NSPEC list !

c absolute values of back- and forward reactions are choosen arbitrarily to 
c ensure a quick (but numerically stable) equilibrium
c if illogical error messages arise, try fiddling with the absolute values
c of equilibrium reactions

            xkef(k,kc,ind_H2O) = xsw*funa(1.0d-5,-6716.d0,k) !schneller wird schlechter!
            xkeb(k,kc,ind_H2O) =
     &           1.0D9* cv2*xgamma(k,1,kc)*xgamma(k,3,kc)

c            xkef(k,kc,ind_HO2) = xsw*2.d5  !Chameides '84
            xkef(k,kc,ind_HO2) = xsw*1.6d5  !Weinstein-loyd and Schwartz, '91
            xkeb(k,kc,ind_HO2) =
     &           1.d10*cv2*xgamma(k,1,kc)*xgamma(k,11,kc)

            xkef(k,kc,ind_ACO2) = xsw*1.8D0
            xkeb(k,kc,ind_ACO2) =
     &           1.0D4*cv2*xgamma(k,1,kc)*xgamma(k,16,kc)   

c            xkef(k,kc,ind_CO2) = xsw*funa(4.3d3,-913.d0,k)
c            xkeb(k,kc,ind_CO2) = 1.0D10* cv2   
            xkef(k,kc,ind_CO2) = xsw*funa(4.3d-2,-913.d0,k)
            xkeb(k,kc,ind_CO2) = 1.0D5*cv2*xgamma(k,1,kc)*xgamma(k,9,kc)  

            xkef(k,kc,ind_HONO) = xsw*funa(5.1d+3,-1260.d0,k)
            xkeb(k,kc,ind_HONO) =
     &           1.0D7*cv2*xgamma(k,1,kc)*xgamma(k,12,kc)

            xkef(k,kc,ind_HNO3) = xsw*funa(1.54d+10,8700.d0,k)
            xkeb(k,kc,ind_HNO3) =
     &           1.0D9*cv2*xgamma(k,1,kc)*xgamma(k,13,kc)  

            xkef(k,kc,ind_HNO4) = xsw*2.0D3
            xkeb(k,kc,ind_HNO4) = 2.0D8 * cv2

            xkef(k,kc,ind_NH3) = xsw*funa(1.7d5,-4325.d0,k) 
            xkeb(k,kc,ind_NH3) =
     &           1.0D10* cv2*xgamma(k,3,kc)*xgamma(k,2,kc)

            xkef(k,kc,ind_HSO3ml1) = xsw*funa(6.0d2,1120.d0,k)
     &           *xgamma(k,5,kc)
            xkeb(k,kc,ind_HSO3ml1) = 1.0D10* cv2 *xgamma(k,1,kc)
     &           *xgamma(k,6,kc)

            xkef(k,kc,ind_H2SO4) = xsw*1.0d12 !Seinfeld, Pandis (1998), p.391
            xkeb(k,kc,ind_H2SO4)=
     &           1.0d9* cv2*xgamma(k,1,kc)*xgamma(k,19,kc)

            xkef(k,kc,ind_HSO4ml1) = xsw*funa(1.02d+6,2720.d0,k) 
     &           *xgamma(k,19,kc)
            xkeb(k,kc,ind_HSO4ml1)=
     &           1.0D8*cv2*xgamma(k,1,kc)*xgamma(k,8,kc)

            xkef(k,kc,ind_SO2) = xsw*funa(1.7d8,2090.d0,k)
            xkeb(k,kc,ind_SO2) =
     &           1.0D10*cv2*xgamma(k,1,kc)*xgamma(k,5,kc)

!           xkef(k,kc,ind_HCHO) = 1.d10*cv2  !Chameides '84   !mechanism changed
!           xkeb(k,kc,ind_HCHO) = xsw*1.d5                    !mechanism changed

            xkef(k,kc,ind_HCl) = xsw*funa(1.7d10,6896.d0,k)
            xkeb(k,kc,ind_HCl) =
     &           1.0D4* cv2*xgamma(k,1,kc)*xgamma(k,14,kc)

            xkef(k,kc,ind_Cl2ml1) = xsw*5.2d4*xgamma(k,15,kc)    !Chameides '84
            xkeb(k,kc,ind_Cl2ml1) =  1.d10*cv2*xgamma(k,14,kc)

c            xkef(k,kc,ind_Cl2) = 1.1D-3   
c            xkeb(k,kc,ind_Cl2) = 2.1D2* cv2 

            xkef(k,kc,ind_HOCl) = xsw*3.2D2   
            xkeb(k,kc,ind_HOCl) =
     &           1.0D10*cv2*xgamma(k,1,kc)*xgamma(k,22,kc)

            xkef(k,kc,ind_HBr) = xsw*1.0D13   
            xkeb(k,kc,ind_HBr) =
     &           1.0D4* cv2*xgamma(k,1,kc)*xgamma(k,24,kc)

            xkef(k,kc,ind_Br2) = xsw*funa(2.95d4,-4068.d0,k)* !Liu et al, 2002, #2109
     &           xgamma(k,25,kc)
            xkeb(k,kc,ind_Br2) = funa(1.17d10,-1812.d0,k)* cv2 *
     &           xgamma(k,24,kc)

            xkef(k,kc,ind_HOBr) = xsw*funa(2.3d1,-3091.d0,k)
            xkeb(k,kc,ind_HOBr) =
     &           1.0D10*cv2*xgamma(k,1,kc)*xgamma(k,26,kc)

            xkef(k,kc,ind_BrCl2ml1) = funa(5.d9,1143.d0,k) * cv2 
     &           *xgamma(k,14,kc)
            xkeb(k,kc,ind_BrCl2ml1) = xsw*1.3D9    *xgamma(k,28,kc)
c no activities used, they "cancel out"
            xkef(k,kc,ind_Br2Clml1) = 5.d9 * cv2!*xgamma(k,24,kc) !if original equilibrium 
            xkeb(k,kc,ind_Br2Clml1) = xsw*2.8D5!*xgamma(k,29,kc)  !       -"-  
            xkef(k,kc,ind_Br2l1) = 5.d9 * cv2 !*xgamma(k,14,kc)   !if original equilibrium
            xkeb(k,kc,ind_Br2l1) = xsw*3.85D9!*xgamma(k,29,kc)    !       -"-
            xkef(k,kc,ind_ICl) = 1.0D11  * cv2 *xgamma(k,14,kc)
            xkeb(k,kc,ind_ICl) = xsw*1.3D9   *xgamma(k,37,kc)
            xkef(k,kc,ind_IBr) = 1.0D11  * cv2  *xgamma(k,24,kc)
            xkeb(k,kc,ind_IBr) = xsw*3.5D8  *xgamma(k,38,kc)
c new (speculative) ICl <--> IBr equilibria, assumed to yield the same
c ICl/IBr ratio as in the BrCl <--> Br2 eqilibria (BrCl/Br2) 
c no activities used, they are not known!! 11.04.01 
           xkef(k,kc,ind_IClBrml1) = 5.d9 * cv2!*xgamma(k,24,kc)
           xkeb(k,kc,ind_IClBrml1) = xsw*2.8D5!*xgamma(k,29,kc)
           xkef(k,kc,ind_I2) = 5.d9 * cv2!*xgamma(k,14,kc)    
           xkeb(k,kc,ind_I2) = xsw*3.85D9!*xgamma(k,29,kc)
!            xkef(k,kc,ind_HIO2) = xsw*2.0D3
!            xkeb(k,kc,ind_HIO2) = 2.0D9 * cv2
            xkef(k,kc,ind_HIO3) = xsw*1.57D4
            xkeb(k,kc,ind_HIO3) = 1.0D5 * cv2

c mercury; no activity coefficients; forward not faster than 1.d14
c
c            xkef(k,kc,ind_HgOHpl1)    = 4.27d14 * cv2           ! #493,  K_eq = 4.27d10
c            xkeb(k,kc,ind_HgOHpl1)    = xsw * 1.d4
c            xkef(k,kc,ind_HgOH2l1)    = 2.6d14 * cv2            ! #4174, K_eq = 2.6d11
c            xkeb(k,kc,ind_HgOH2l1)    = xsw * 1.d3
c            xkef(k,kc,ind_HgSO3l1)    = 5.01d14 * cv2           ! #493,  K_eq = 5.01d12
c            xkeb(k,kc,ind_HgSO3l1)    = xsw * 1.d2
c            xkef(k,kc,ind_HgSO322ml1) = 2.5d14 * cv2            ! #4158, K_eq = 2.5d11
c            xkeb(k,kc,ind_HgSO322ml1) = xsw * 1.d3
c            xkef(k,kc,ind_HgOHCll1)   = 2.69d14 * cv2            ! #4158, K_eq = 2.69d7
c            xkeb(k,kc,ind_HgOHCll1)   = xsw * 1.d7
c            xkef(k,kc,ind_HgClpl1)    = 2.0d14 * cv2             ! #4174, K_eq = 2.0d7
c            xkeb(k,kc,ind_HgClpl1)    = xsw * 1.d7
c            xkef(k,kc,ind_HgCl2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d6
c            xkeb(k,kc,ind_HgCl2l1)    = xsw * 1.d8
c            xkef(k,kc,ind_HgCl3ml1)   = 6.7d8 * cv2             ! #4174, K_eq = 6.7d0
c            xkeb(k,kc,ind_HgCl3ml1)   = xsw * 1.d8
c            xkef(k,kc,ind_HgCl42ml1)  = 1.3d9 * cv2             ! #4174, K_eq = 1.3d1
c            xkeb(k,kc,ind_HgCl42ml1)  = xsw * 1.d8 
c            xkef(k,kc,ind_HgOHBrl1)   = 2.69d14 * cv2            ! #4158, assumed
c            xkeb(k,kc,ind_HgOHBrl1)   = xsw * 1.d7
c            xkef(k,kc,ind_HgBrpl1)    = 1.1d14 * cv2             ! #4174, K_eq = 1.1d9
c            xkeb(k,kc,ind_HgBrpl1)    = xsw * 1.d5
c            xkef(k,kc,ind_HgBr2l1)    = 2.5d14 * cv2             ! #4174, K_eq = 2.5d8
c            xkeb(k,kc,ind_HgBr2l1)    = xsw * 1.d6
c            xkef(k,kc,ind_HgBr3ml1)   = 1.5d10 * cv2             ! #4174, K_eq = 1.5d2
c            xkeb(k,kc,ind_HgBr3ml1)   = xsw * 1.d8
c            xkef(k,kc,ind_HgBr42ml1)  = 2.3d10 * cv2             ! #4174, K_eq = 2.3d1
c            xkeb(k,kc,ind_HgBr42ml1)  = xsw * 1.d9

            else ! conv2(kc,k) = 0
               xkef(k,kc,:) = 0. ! jjb xsw was set to 0. in this case, leading to xkef = 0.
               xkeb(k,kc,:) = 0. ! jjb cv2 was tested equal to 0 in this case, leading to xkeb = 0.
            end if

         enddo
      enddo

      end subroutine equil_co_a



c
c------------------------------------------------------------------------
c


!     subroutine konc (ij) ! jjb dummy argument not used
      subroutine konc      ! jjb removed
c new concentrations of chemical species in liquid phases
c due to changes in particle size distribution
c
c all changes are calculated via changes of the volume of
c the nkc size bins

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkc

      implicit double precision (a-h,o-z)

      common /blck06/ kw(nka),ka
      common /blck07/ part_o_a(nka,n),part_o_d(nka,n),
     &                part_n_a(nka,n),part_n_d(nka,n),pntot(nkc,n)
      common /blck08/ vol1_a(nka,n),vol1_d(nka,n),vol2(nkc,n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)

!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka

c vol2 and vol1 old liquid volume in class (1-nkc) and row/class (1-nkc,1-nka)
c (um^3/cm^3), used in SR konc to shift moles from aerosol to drop
c and vice versa. 
c part_o and part_n old and new part. conc. in row/class (1-nkc,1-nka) (cm^-3)

      do k=2,nf

c small (dry) particles
         do ia=1,ka
            dp_1=part_o_a(ia,k)-part_n_a(ia,k)
            dp_3=part_o_d(ia,k)-part_n_d(ia,k)
            if (dabs(dp_1+dp_3).gt.1.d-10) print 
     &           *,'Warning SR konc dp_1 > dp_3',k,ia,dp_1,dp_3
c search bin that loses moles:
            ii=1
c            if (dp_1.ge.1.d-10) ii=1
            if (dp_3.ge.1.d-10) ii=3
c no change if number of growing particles is too small:
            xs=1.d0
            if (dabs(dp_1).lt.1.d-10) xs=0.d0 
c mol/m^3(air) -1-> mol/l(aq) -2-> mol/m^3(air,row) -3-> mol/part(row)
c -4-> mol/m^3(air,change,row)
c mol/m^3(air) -1-> mol/l(aq):
c     mol/m^3(air)/(vol2*1.d-12) = mol/m^3(aq)
c     mol/m^3(aq)*1.d-3 --> mol/l(aq)
c mol/l(aq) -2-> mol/m^3(air,row):
c     mol/l(aq) *vol1_@(ia,k)*1.d-9 = mol/m^3(air,row)
c mol/m^3(air,row) -3-> mol/part(row):
c     mol/m^3(air,row)*1.d-6 cm^3(air,row)/part = mol/part(row) 
c mol/part(row) -4-> mol/m^3(air,change,row):
c     mol/part(row)*1.d6*part/cm^3(air,row,change) = mol/m^3(air,row,change)
            if (ii.eq.1) then
               jj=3
               if (vol2(ii,k).gt.0..and.part_o_a(ia,k).gt.0.) then
                  delta=vol1_a(ia,k)/vol2(ii,k)*
     &              dp_1/part_o_a(ia,k)*xs
               else
                  delta=0.d0
               end if
            else ! ii == 3
               jj=1
               if (vol2(ii,k).gt.0..and.part_o_d(ia,k).gt.0.) then
                  delta=vol1_d(ia,k)/vol2(ii,k)*
     &              dp_3/part_o_d(ia,k)*xs
               else
                  delta=0.d0
               end if
            endif
            if (delta.lt.0.) print *,k,ia,delta,'Warning SR konc s <'
            if (delta.gt.1.) print *,k,ia,delta,'Warning SR konc s >'
            if (delta.eq.0.) goto 1111
            do l=1,j2
               del=sl1(l,ii,k)*delta
               sl1(l,ii,k)=dmax1(0.d0,sl1(l,ii,k)-del)
               sl1(l,jj,k)=dmax1(0.d0,sl1(l,jj,k)+del)
            enddo
            do l=1,j6
               del=sion1(l,ii,k)*delta
               sion1(l,ii,k)=dmax1(0.d0,sion1(l,ii,k)-del)
               sion1(l,jj,k)=dmax1(0.d0,sion1(l,jj,k)+del)
            enddo
 1111       continue
         enddo
c large (dry) particles
         do ia=ka+1,nka
            dp_2=part_o_a(ia,k)-part_n_a(ia,k)
            dp_4=part_o_d(ia,k)-part_n_d(ia,k)
            if (dabs(dp_2+dp_4).gt.1.d-10) print 
     &           *,'Warning SR konc dp_2 > dp_4',k,ia,dp_2,dp_4
            ii=2
c            if (dp_2.ge.1.d-10) ii=2
            if (dp_4.ge.1.d-10) ii=4
            xs=1.d0
            if (dabs(dp_2).lt.1.d-10) xs=0.d0

            if (ii.eq.2) then
               jj=4
               if (vol2(ii,k).gt.0..and.part_o_a(ia,k).gt.0.) then
                  delta=vol1_a(ia,k)/vol2(ii,k)*
     &              dp_2/part_o_a(ia,k)*xs
               else
                  delta=0.d0
               end if
            else ! ii == 4
               jj=2
               if (vol2(ii,k).gt.0..and.part_o_d(ia,k).gt.0.) then
                  delta=vol1_d(ia,k)/vol2(ii,k)*
     &              dp_4/part_o_d(ia,k)*xs
               else
                  delta=0.d0
               end if
            endif

            if (delta.lt.0.) print *,k,ia,delta,'SR konc l <'
            if (delta.gt.1.) print *,k,ia,delta,'SR konc l >'
            if (delta.eq.0.) goto 1112

            do l=1,j2
               del=sl1(l,ii,k)*delta
               sl1(l,ii,k)=dmax1(0.d0,sl1(l,ii,k)-del)
               sl1(l,jj,k)=dmax1(0.d0,sl1(l,jj,k)+del)
            enddo
            do l=1,j6
               del=sion1(l,ii,k)*delta
               sion1(l,ii,k)=dmax1(0.d0,sion1(l,ii,k)-del)
               sion1(l,jj,k)=dmax1(0.d0,sion1(l,jj,k)+del)
            enddo
 1112       continue
         enddo ! kc
         

         
         do l=1,j2
            if (pntot(3,k).lt.1.d-7) then
               sl1(l,1,k)=sl1(l,1,k)+dmax1(0.d0,sl1(l,3,k))
               sl1(l,3,k)=0.
            endif
            if (pntot(4,k).lt.1.d-7) then
               sl1(l,2,k)=sl1(l,2,k)+dmax1(0.d0,sl1(l,4,k))
               sl1(l,4,k)=0.
            endif
         enddo

         do l=1,j6
            if (pntot(3,k).lt.1.d-7) then
               sion1(l,1,k)=sion1(l,1,k)+dmax1(0.d0,sion1(l,3,k))
               sion1(l,3,k)=0.
            endif
            if (pntot(4,k).lt.1.d-7) then
               sion1(l,2,k)=sion1(l,2,k)+dmax1(0.d0,sion1(l,4,k))
               sion1(l,4,k)=0.
            endif
         enddo

      enddo ! k

      end subroutine konc

c
c-----------------------------------------------------------------------------------
c

      subroutine init_konc 
c Init loading of aerosols with ions.
c If rH>0.7 aerosols get activated and liquid chemistry is started.
c The ions are transported vertically by turbulence and sedimentation
c independent of activation.
c      double precision ap1o,ap2o,ap2n,apo,apn

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j3,
     &     j6,
     &     n,
     &     nka,
     &     nkt,
     &     nkc

      implicit double precision (a-h,o-z)

      common /blck06/ kw(nka),ka
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /blck78/ sa1(nka,j2),sac1(nka,j2)
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /nucfeed/ ifeed

      dimension ap(n,nka)

c sa1 in mole/particle; sl1, sion1 in mole m**-3
c ap total number of aerosols in cm**-3 
c of sulfate and sea salt aerosol regime

      do k=2,n-1

         if (ifeed.eq.2) then
           ial = 2
         else
           ial = 1
         endif
         do ia=ial,nka
            ap(k,ia)=0.
c            do jt=1,kg
            do jt=1,nkt
               ap(k,ia)=ap(k,ia)+ff(jt,ia,k)
            enddo
            if (ap(k,ia).eq.0.) print *,ia,k,'init(ap)=0'
         enddo
c no3-, nh4+ and so4= due to nucleation scavenging and evaporation
         do ia=1,ka
            sion1( 2,1,k)=sion1( 2,1,k)+ap(k,ia)*sa1(ia, 2)*1.d6
            sion1( 8,1,k)=sion1( 8,1,k)+ap(k,ia)*sa1(ia, 8)*1.d6
            sion1(13,1,k)=sion1(13,1,k)+ap(k,ia)*sa1(ia,13)*1.d6
            sion1(17,1,k)=sion1(19,1,k)                         !mixing tracer
            sion1(19,1,k)=sion1(19,1,k)+ap(k,ia)*sa1(ia,19)*1.d6
c NO3-, NH4-, SO4=, HSO4- : j6 used in sion1, sa1 
c Br-, HCO3-, I-, IO3-, Cl-: j6 used in sion1, sa1
         enddo
         if (sion1(8,1,k).eq.0.) print *,k,'init(so4=)=0'

         do ia=ka+1,nka
            sion1( 8,2,k)=sion1( 8,2,k)+ap(k,ia)*sa1(ia, 8)*1.d6 !SO4=
            sion1( 9,2,k)=sion1( 9,2,k)+ap(k,ia)*sa1(ia, 9)*1.d6 !HCO3-
            sion1(13,2,k)=sion1(13,2,k)+ap(k,ia)*sa1(ia,13)*1.d6 !NO3-
            sion1(14,2,k)=sion1(14,2,k)+ap(k,ia)*sa1(ia,14)*1.d6 !Cl-
            sion1(17,2,k)=sion1(14,2,k)                          !mixing tracer
            sion1(20,2,k)=sion1(20,2,k)+ap(k,ia)*sa1(ia,20)*1.d6 !Na+; inert
            sion1(24,2,k)=sion1(24,2,k)+ap(k,ia)*sa1(ia,24)*1.d6 !Br-
            sion1(34,2,k)=sion1(34,2,k)+ap(k,ia)*sa1(ia,34)*1.d6 !I-
            sion1(36,2,k)=sion1(36,2,k)+ap(k,ia)*sa1(ia,36)*1.d6 !IO3-
            sl1(j2-j3+4,2,k)=sl1(j2-j3+4,2,k)+ap(k,ia)*
     &           sa1(ia,j2-j3+4)*1.d6 !DOM
         enddo
         if (sion1(14,2,k).eq.0.) print *,k,'init(cl-)=0'
      enddo

      end subroutine init_konc


c
c-----------------------------------------------------------------------------------
c

      subroutine aer_source (box,dd,z_mbl,n_bl)
c calculation of sea salt aerosol source

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j3,
     &     j6,
     &     n,
     &     nka,
     &     nkt,
     &     nkc

      implicit double precision (a-h,o-z)
      logical box

      logical mona,smith
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /blck06/ kw(nka),ka
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
     &       /blck78/ sa1(nka,j2),sac1(nka,j2)
      common /sss/ brsss,clsss,xnasss
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
!     dimension   rnw(nka) ! jjb not used

c choose which version to pick
      mona=.true.      ! Monahan et al., 1986
      smith=.false.    ! Smith et al, 1993

      if (box) then
         k_in = n_bl
         d_z  = z_mbl
      else
         k_in = 2
         d_z  = deta(2)
      endif
c this is also defined for box:
      u10  = 0.5*(dsqrt(u(3)**2+v(3)**2)+dsqrt(u(2)**2+v(2)**2))

      if (smith) then
         a1=10**(0.0676*u10+2.43)
         a2=10**(0.959*u10**0.5-1.476)
         f1=3.1
         f2=3.3
         r01=2.1
         r02=9.2
      endif

      do ia=ka+1,nka
c compare SR vgleich in str.f
c aerosols in equilibrium with ambient rH:
         feu(k_in)=dmin1(feu(k_in),0.99999d0)
c equilibrium radius
         a0=a0m/t(k_in)
         b0=b0m(ia)*2.
c b0=b0m*rho3/rhow; rho3=2000; rhow=1000
         rg=rgl(rn(ia),a0,b0,feu(k_in))
         eg=4.d-09*pi/3.*(rg**3-rn(ia)**3)
         do jt=1,nkt
c find jt that corresponds to particle size at ambient rH
c            if (rq(jt,ia).ge.rg) then
            if (eg.le.ew(jt)) then
c              rr=rq(jt,ia)
c source functions are for rH=80%
               rr=rgl(rn(ia),a0,b0,0.8d0)
c find jt-index that corresponds to rr for this dry aerosol radius
               jt_low = -1
               do jtt=1,nkt
                  if (rq(jtt,ia).le.rr) jt_low = jtt
               enddo
c               print *,ia,jt_low,rr,rq(jt_low,ia)
              if (mona) then
c aerosol source after Monahan et al. 86 cited in Gong et al, 97, JGR, 102, p3805-3818
c  "log" is log10
c               bb=(0.380-log(rr))/0.65
               bb=(0.380-log10(rr))/0.65
c               print *,'bb',bb
               df1=1.373*u10**3.41*rr**(-3.)*(1. + 0.057 *rr**1.05) 
     &              * 10**(1.19*exp(-bb**2))
!              if (rr.gt.10..and.rr.lt.75.) df22=8.6d-6*exp(2.08*u10)*
!    &              rr**(-2)
!              if (rr.gt.75..and.rr.lt.100.) df23=4.83d-2*exp(2.08*u10)*
!    &              rr**(-4)
!              if (rr.gt.100.) df24=8.6d6*exp(2.08*u10)*rr**(-8)
               df=df1!+df22+df23+df24
              endif
              if (smith) then
c aerosol source after Smith et al., 93, QJRMS,119,809-824
               df1=a1*dexp(-f1*(dlog(rr/r01)**2))
               df2=a2*dexp(-f2*(dlog(rr/r02)**2))
               df=df1+df2
              endif
c               print *,df,rr,rq(jt,ia)-rq(jt,ia-1)
c df in m^-2 um^-1 s^-1, convert to: cm^-3 s^-1
c m^-2 --> m^-3: 1/dz, m^-3 --> cm^-3: 10-6, integrate over r
!              if (jt.eq.1) then
              if (jt_low.eq.1) then
c                  df=df*(rq(jt+1,ia)-rq(jt,ia))/d_z*1.d-6
                 df=df*(rq(jt_low+1,ia)-rq(jt_low,ia))/d_z*1.d-6
              else
c                  df=df*(rw(jt,ia)-rw(jt-1,ia))/d_z*1.d-6
                 df=df*(rw(jt_low,ia)-rw(jt_low-1,ia))/d_z*1.d-6
              endif
c               print *,ia,df,df*86400
c new distribution
c change number conc and aerosol conc (Br-, I-, IO3-, Cl-, Na+, DOM)
              ff(jt,ia,k_in)=ff(jt,ia,k_in)+df*dd
              sion1( 8,2,k_in)=sion1( 8,2,k_in)+df*dd*sa1(ia, 8)*1.d6 !SO4=
              sion1( 9,2,k_in)=sion1( 9,2,k_in)+df*dd*sa1(ia, 9)*1.d6 !HCO3-
              sion1(13,2,k_in)=sion1(13,2,k_in)+df*dd*sa1(ia,13)*1.d6 !NO3-
              sion1(14,2,k_in)=sion1(14,2,k_in)+df*dd*sa1(ia,14)*1.d6 !Cl-
              sion1(20,2,k_in)=sion1(20,2,k_in)+df*dd*sa1(ia,20)*1.d6 !Na+, ion balance
              sion1(24,2,k_in)=sion1(24,2,k_in)+df*dd*sa1(ia,24)*1.d6 !Br-
              sion1(34,2,k_in)=sion1(34,2,k_in)+df*dd*sa1(ia,34)*1.d6 !I-
              sion1(36,2,k_in)=sion1(36,2,k_in)+df*dd*sa1(ia,36)*1.d6 !IO3-
              sl1(j2-j3+4,2,k_in)=sl1(j2-j3+4,2,k_in)+df*dd*
     &             sa1(ia,j2-j3+4)*1.d6 !DOM
              brsss=brsss+df*dd*sa1(ia,24)*1.d6
              clsss=clsss+df*dd*sa1(ia,14)*1.d6
              xnasss=xnasss+df*dd*sa1(ia,20)*1.d6
c               print *,jt,f(2,ia,jt),df
              goto 1000
            endif
         enddo
 1000    continue

      enddo

      end subroutine aer_source


c
c-----------------------------------------------------------------------
c

      subroutine plo_ppH
c print NH3,NH4+,HCl,Cl-,HNO3,NO3-,SO2,H2SO4,HSO4-,SO4=
c to calculate the potential pH
c use for parameterisation development for global models

      USE gas_common, ONLY :
! Imported Array Variables with intent (in):
     &     s1

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nkc


      implicit none

! Local scalars:
      character (len=8),parameter :: fname = 'ppHa.out'
      integer k

! Local arrays:
      double precision i0
      dimension i0(n,25)

! Common blocks:
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct
!- End of header ---------------------------------------------------------------

      do k=2,nf
            i0(k,1)=s1(4,k)
            i0(k,2)=s1(30,k)
            i0(k,3)=s1(3,k)
            i0(k,4)=s1(5,k)
            i0(k,5)=s1(6,k)

            i0(k,6)=sl1(4,1,k)
            i0(k,7)=sl1(30,1,k)
            i0(k,8)=sl1(3,1,k)
            i0(k,9)=sl1(5,1,k)
            i0(k,10)=sl1(6,1,k)
            i0(k,11)=sl1(4,2,k)
            i0(k,12)=sl1(30,2,k)
            i0(k,13)=sl1(3,2,k)
            i0(k,14)=sl1(5,2,k)
            i0(k,15)=sl1(6,2,k)

            i0(k,16)=sion1(2,1,k)
            i0(k,17)=sion1(14,1,k)
            i0(k,18)=sion1(13,1,k)
            i0(k,19)=sion1(19,1,k)
            i0(k,20)=sion1(8,1,k)
            i0(k,21)=sion1(2,2,k)
            i0(k,22)=sion1(14,2,k)
            i0(k,23)=sion1(13,2,k)
            i0(k,24)=sion1(19,2,k)
            i0(k,25)=sion1(8,2,k)
      enddo

 3000 continue
      open (73, file=fname,status='unknown',form='unformatted',
     & position='append',err=3000)
      write (73) lday,lst,lmin,i0
      close (73)

      end subroutine plo_ppH


c
c-----------------------------------------------------------
c

      subroutine kpp_driver (box,dd_ch,n_bl)
c interface between MISTRA and the KPP gas phase chemistry

      USE config, ONLY :
     &     halo,
     &     iod,
     &     neula

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE gas_common, ONLY :
! Imported Array Variables with intent (inout):
     &     s1,
     &     s3

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nkc,
     &     nphrxn

      implicit double precision (a-h,o-z) 
!     logical halo,iod,box,new_a,new_c,cloud,short,ros3 ! jjb ros3 removed, thus new_a, new_c & short as well
      logical box,cloud

      common /blck01/ am3(n),cm3(n)

      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck13/ conv2(nkc,n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /band_rat/ photol_j(nphrxn,n)
      common /kinv_i/ kinv
      common /cb_1/ air_cc,te,h2oppm,pk
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_l1/ cloud(nkc,n)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated
!     common /ph_r/ ph_rat(47) ! jjb ! jjb test include in _a_g_t common blocks
      dimension ph_rat(nphrxn)
      common /kpp_eul/ xadv(10),nspec(10)

      airmolec=6.022d+20/18.0
c calculate u0 for prelim photolysis
      rlat=alat*1.745329d-02
      rdec=declin*1.745329d-02
      zeit=lst*3600.+dfloat(lmin-1)*60.
c      zeit=zeit+dtg  !quite exact, but photolysis rates are calculated in 
c greater intervals and variable dtg is known only in the old chemical module
!      horang=7.272205d-05*zeit-3.1415927
      horang=7.272205d-05*zeit-pi
      u00=dcos(rdec)*dcos(rlat)*dcos(horang)+dsin(rdec)*dsin(rlat)
      ru0=6371.*u00
      u0=8./(dsqrt(ru0**2+102000.)-ru0)

      call mass_ch
c loop over all vertical layers
c      n_max=n
      n_max=n-1
      n_min=2
      if (box) then
         n_max=n_bl
         n_min=n_bl
      endif

c eliminate negative values
      where (s1 < 0.d0) s1=0.d0 ! jjb new using where construct for array s1
      where (s3 < 0.d0) s3=0.d0 ! jjb new using where construct for array s3
      !here (sl1  < 0.d0)   sl1=0.d0
      !here (sion1 < 0.d0) sion1=0.d0 ! at the moment, done in each mechanism

      do k=n_min,n_max !c aer#
!         if (k.eq.2) short=.false. !c aer# ! jjb ros2 only
c cloud#      do k=n_max,2,-1 !c cloud#
c cloud#         if (k.eq.n_max) short=.false. !c cloud#

c define temp, H2O, air, ..
c air, h2o in mlc/cm^3 are used for 3rd order reactions that are
c formulated as quasi-2nd order
         te=t(k)
         air_cc=cm3(k)                      ![air] in mlc/cm^3
         air=am3(k)                         ![air] in mol/m^3
!        co=am3(k,2)                        ![co]  in mol/m^3 ! jjb removed, unused
         h2o=xm1(k)*rho(k)/1.8d-2           ![h2o] in mol/m^3
         h2o_cc=xm1(k)*airmolec*rho(k)      ![h2o] in mlc/cm^3
         h2oppm=h2o_cc*1.d6/air_cc        ![h2o] in ppm=umol/mol
         pk=p(k)
         dt_ch=dd_ch
         tkpp=time
!         if (k.lt.nf) then ! conv2 is dimension n now, and initialised to 0
            cvv1=conv2(1,k)
            cvv2=conv2(2,k)
            cvv3=conv2(3,k)
            cvv4=conv2(4,k)
!         else ! jjb check that there is no pb with index nf
!            cvv1=0.
!            cvv2=0.
!            cvv3=0.
!            cvv4=0.
!         endif

c photolysis rates
         if (u0.ge.3.48d-2) then
            do i=1,nphrxn
               ph_rat(i)=( photol_j(i,k)+photol_j(i,k-1) ) / 2 ! jjb photol_j defined for levels, not for layers
            enddo
c end of prelim j rates 
         else
            do i=1,nphrxn
               ph_rat(i)=0.
            enddo
         end if

c set halogen rates to zero (if wanted)
         xhal=1.
         xiod=1.
         if (.not.halo) then
            xhal=0.
            xiod=0.
         endif
         if (.not.iod) xiod=0.
c set liquid rates to zero
         xliq1=1.
         xliq2=1.
         xliq3=1.
         xliq4=1.
c cm is switch now (not cw)
c avoid index out of bounds in cm:
         if (k.ge.nf) then
            xliq1=0.
            xliq2=0.
            xliq3=0.
            xliq4=0.
         else
            if (cm(1,k).eq.0.) xliq1=0.
            if (cm(2,k).eq.0.) xliq2=0.
            if (cm(3,k).eq.0.) xliq3=0.
            if (cm(4,k).eq.0.) xliq4=0.
         endif

         if (xliq1.eq.1..and..not.cloud(1,k)) then
            cloud(1,k)=.true.
            print *,'new aerosol layer,l1 : ',k
         endif
         if (xliq2.eq.1..and..not.cloud(2,k)) then
            cloud(2,k)=.true.
            print *,'new aerosol layer,l2 : ',k
         endif
         if (xliq3.eq.1..and..not.cloud(3,k)) then
            cloud(3,k)=.true.
            print *,'new cloud layer,l3 : ',k
         endif
         if (xliq4.eq.1..and..not.cloud(4,k)) then
            cloud(4,k)=.true.
            print *,'new cloud layer,l4 : ',k
         endif
         if (xliq1.eq.0..and.cloud(1,k))  cloud(1,k)=.false.
         if (xliq2.eq.0..and.cloud(2,k))  cloud(2,k)=.false.
         if (xliq3.eq.0..and.cloud(3,k))  cloud(3,k)=.false.
         if (xliq4.eq.0..and.cloud(4,k))  cloud(4,k)=.false.

c set heterogeneous rates to zero
         xhet1=1.
         xhet2=1.
         if (xliq1.eq.1.) xhet1=0.
         if (xliq2.eq.1.) xhet2=0.

c advection if eulerian view (neula=0), xadv in mol/(mol*day)

         if (neula.eq.0) then
            if (k.le.kinv) then
               do nisp=1,10
                  if (nspec(nisp).ne.0) 
     &   s1(nspec(nisp),k)=s1(nspec(nisp),k)+xadv(nisp)*dt_ch*air/86400.
               enddo
            endif
         endif

!         print('(i4,6f4.1)'),k,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2
c call kpp chemical integrator 
            if (xliq1.eq.1..or.xliq2.eq.1.) then
               if (xliq3.eq.1..or.xliq4.eq.1.) then
c gas + aerosol + cloud
c                  print *,k," tot ",tkpp,dt_ch
!                 call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4
!    &                 ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,air,co,h2o) ! jjb
!                 call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4
!    &     ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2,air,co,h2o) ! jjb working version
!                  call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4
!     &     ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2,air,co,h2o
!     &     ,ph_rat)                                                   ! jjb test version to handle ph_rat another way
                  call tot_drive (tkpp,dt_ch,k,cvv1,cvv2,cvv3,cvv4
     &     ,xhal,xiod,xliq1,xliq2,xliq3,xliq4,xhet1,xhet2,air,h2o
     &     ,ph_rat)                                                   ! jjb co removed
!         print *,'layer nb',k,'aqueous chem? TOT mech'
!         print *,'xliq1-4=',xliq1,xliq2,xliq3,xliq4
!         print *,'xhet1-2=',xhet1,xhet2
!         print *,'cloud(1-4,k)=',cloud(1,k),cloud(2,k),cloud(3,k)
!     & ,cloud(4,k)

               else
c gas + aerosol               
c                  print *,k," aer ",tkpp,dt_ch
!                 call aer_drive (tkpp,dt_ch,k,cvv1,cvv2,xhal,
!    &                 xiod,xliq1,xliq2,xhet1,xhet2,air,co,h2o) ! jjb working version
!                  call aer_drive (tkpp,dt_ch,k,cvv1,cvv2,xhal,
!     &                 xiod,xliq1,xliq2,xhet1,xhet2,air,co,h2o
!     &     ,ph_rat)                                             ! jjb test version to handle ph_rat another way
                  call aer_drive (tkpp,dt_ch,k,cvv1,cvv2,xhal,
     &                 xiod,xliq1,xliq2,xhet1,xhet2,air,h2o
     &     ,ph_rat)                                             ! jjb co removed
!         print *,'layer nb',k,'aqueous chem? AER mech'
!         print *,'xliq1-4=',xliq1,xliq2,xliq3,xliq4
!         print *,'xhet1-2=',xhet1,xhet2
!         print *,'cloud(1-4,k)=',cloud(1,k),cloud(2,k),cloud(3,k)
!     & ,cloud(4,k)
               endif   
            else
c gas only            
c het reactions on dry aerosol always on!
               xhet1 = 1.
               xhet2 = 1.
c               print *,k," gas ",tkpp,dt_ch
!              call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2,conv1, ! jjb working version
!    &              air,co,h2o)
!               call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2,conv1, ! jjb test version to handle ph_rat another way
!     &              air,co,h2o,ph_rat)
!              call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2,conv1, ! jjb co removed
!    &              air,h2o,ph_rat)
!               write(1717,*)'Gas layer',time,k
               call gas_drive (tkpp,dt_ch,k,xhal,xiod,xhet1,xhet2,       ! jjb conv1 removed
     &              air,h2o,ph_rat)!         print *,'layer nb',k,'aqueous chem? GAS mech'
!         print *,'xliq1-4=',xliq1,xliq2,xliq3,xliq4
!         print *,'xhet1-2=',xhet1,xhet2
!         print *,'cloud(1-4,k)=',cloud(1,k),cloud(2,k),cloud(3,k)
!     & ,cloud(4,k)
            endif

      enddo ! k

c eliminate negative values
      where (s1 < 0.d0) s1=0.d0 ! jjb new using where construct for array s1
      where (s3 < 0.d0) s3=0.d0 ! jjb new using where construct for array s3

      sl1(:,:,:) = max(0.d0,sl1(:,:,:))
      sion1(:,:,:) = max(0.d0,sion1(:,:,:))

      if (lmin.eq.10) call ionbalance (box,n_bl)

      end subroutine kpp_driver

c
c----------------------------------------------------- 
c

      subroutine ionbalance (box,n_bl)
c check ion balance

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     n,
     &     nkc

      implicit none

      logical, intent(in) :: box
      integer, intent(in) :: n_bl

      integer :: k, kc
      integer :: n_min, n_max

      double precision :: xpos(nkc),xneg(nkc)

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1


      write (103,*) lday,lst,lmin,' aerosol'
      write (104,*) lday,lst,lmin,' droplet'

      n_min=1
      n_max=n
      if (box) then 
         n_min=n_bl
         n_max=n_bl
      endif

      do k=n_min,n_max
         do kc=1,nkc
            xpos(kc)=sion1(1,kc,k)+sion1(2,kc,k)+sion1(20,kc,k)
            xneg(kc)=sion1(3,kc,k)+sion1(4,kc,k)+sion1(5,kc,k)+
     &           2*sion1(6,kc,k)+sion1(7,kc,k)+2*sion1(8,kc,k)
     &           +sion1(9,kc,k)+sion1(10,kc,k)+sion1(11,kc,k)
     &           +sion1(12,kc,k)+sion1(13,kc,k)+sion1(14,kc,k)
     &           +sion1(15,kc,k)+sion1(16,kc,k)+sion1(19,kc,k)
     &           +sion1(21,kc,k)
     &           +sion1(22,kc,k)+sion1(23,kc,k)+sion1(24,kc,k)
     &           +sion1(25,kc,k)+sion1(26,kc,k)+sion1(27,kc,k)
     &           +sion1(28,kc,k)+sion1(29,kc,k)+sion1(30,kc,k)
     &           +sion1(31,kc,k)+sion1(32,kc,k)+sion1(33,kc,k)
     &           +sion1(34,kc,k)+sion1(35,kc,k)+sion1(36,kc,k)
     &           +sion1(37,kc,k)+sion1(38,kc,k)+sion1(39,kc,k)
         enddo
         write (103,101) k,xpos(1),xneg(1),xpos(1)-xneg(1),
     &                   xpos(2),xneg(2),xpos(2)-xneg(2)
         write (104,101) k,xpos(3),xneg(3),xpos(3)-xneg(3),
     &                   xpos(4),xneg(4),xpos(4)-xneg(4)
      enddo
 101  format (i3,6d16.8)

      end subroutine ionbalance
c
c------------------------------------------------------------
c

      subroutine dry_cw_rc (nmax)
c calculates LWC and mean radius for "dry" aerosol

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc


      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: nmax

! Local parameters:
      double precision, parameter :: xpi = 4./3.*pi

! Local scalars:
      double precision :: cwd1, cwd2
      double precision :: rcd1, rcd2
      double precision :: x0
      integer :: ia, ial, jt, k

! Common blocks:
      common /blck06/ kw(nka),ka
      integer kw, ka

      common /blck11/ rcd(nkc,n)
      double precision rcd

      common /blck12/ cwd(nkc,n),cm(nkc,n)
      double precision cwd, cm

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka), ! only rq is used
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw, ew, rn, rw, en, e, dew, rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n) ! only ff is used
      double precision ff, fsum
      integer nar

      common /nucfeed/ ifeed
      integer ifeed

!- End of header ---------------------------------------------------------------

      if (ifeed.eq.2) then
         ial = 2
      else
         ial = 1
      endif

c rcd(n,2): mean radius of aerosols in m -------
      do k=nf+1,nmax
c calculate rc for all rel. humidities drier than crystallization
         rcd1=0.
         rcd2=0.
         cwd1=0.
         cwd2=0.

c here TOTAL particle volume is used for calculating LWC
c this is correct only for completely soluble aerosol
c small aerosol

         do ia=ial,ka
            do jt=1,kw(ia)
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cwd1=cwd1+x0
               rcd1=rcd1+x0*rq(jt,ia)
            enddo
         enddo
c large aerosol
         do ia=ka+1,nka
            do jt=1,kw(ia)
               x0=ff(jt,ia,k)*xpi*rq(jt,ia)**3
               cwd2=cwd2+x0
               rcd2=rcd2+x0*rq(jt,ia)
            enddo
         enddo
c conversion: um^3/cm^3 --> m^3(aq)/m^3(air):10^-12
c           : um        --> m               :10^-6
         if (cwd1.gt.0.d0) then
            rcd(1,k)=rcd1/cwd1*1.d-6
         else
            rcd(1,k) = 0.d0
         end if

         if (cwd2.gt.0.d0) then
            rcd(2,k)=rcd2/cwd2*1.d-6
         else
            rcd(2,k) = 0.d0
         end if

         cwd(1,k)=cwd1*1.d-12
         cwd(2,k)=cwd2*1.d-12

      end do

      end subroutine dry_cw_rc


c
c------------------------------------------------------------
c

      subroutine dry_rates_g (tt,pp,nmax)
c calculates kmt for heterogeneous reactions on dry aerosol
c for ALL levels without real aerosol chemistry
c 1: sulfate aerosol, 2: seasalt aerosol

c to be included in gas-phase mechanism

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nkc

      implicit none

      include 'gas_Parameters.h'     !additional common blocks and other definitions
     
      double precision, intent(in) :: tt(n) ! Temperature
      double precision, intent(in) :: pp(n) ! Pressure
      integer, intent(in) :: nmax           ! Maximum layer for calculations

      integer,parameter :: ndr=4 ! number of species
      integer idr(ndr)
      integer k,kc,l
      double precision x1, FCT
      double precision xgamma(NSPEC,2),freep(n),vmean(NSPEC,n)
      common /blck11/ rcd(nkc,n) ! average particle radius in chemical bin
      double precision rcd       !  note that in the context of this subroutine, this
                                 !  radius is considered as a "dry" radius, thus labeled rcd

      common /kpp_dryg/ xkmtd(n,2,NSPEC),henry(n,NSPEC),xeq(n,NSPEC)
      double precision xkmtd, henry, xeq

      double precision func, funa, func3 ! local function
      double precision a,a0,b0           !   and their variables
      integer k0                         !   ...

c      data idr/ind_HNO3,ind_N2O5,ind_BrNO3,ind_ClNO3,ind_HOBr/
!     data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4,ind_HCl/     ! jjb HCl is not handled here (but it should be, probably !)
      data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4/           ! jjb removed at the moment

      func(a,k)=dsqrt(tt(k)/a)*4.60138
      funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))
      func3(a0,b0,k0)=a0*exp(b0*((1/tt(k0))-3.3557d-3))

c calculate freep for all heights:
      do k=2,nmax
         freep(k)=2.28d-5 * tt(k) / pp(k)
      enddo

c define gamma's for all species ----
      xgamma(ind_HNO3,1)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_HNO3,2)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_N2O5,1)  = 0.02       !estimated from JPL 2003
      xgamma(ind_N2O5,2)  = 0.02       !estimated from JPL 2003
      xgamma(ind_NH3,1)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_NH3,2)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_H2SO4,1) = 0.1        !estimated from alpha=0.65
      xgamma(ind_H2SO4,2) = 0.1        !estimated from alpha=0.65
c      xgamma(ind_HCl,1)   = 
c      xgamma(ind_HCl,2)   =
c      xgamma(ind_BrNO3,1) = 0.3
c      xgamma(ind_BrNO3,2) = 0.3
c      xgamma(ind_ClNO3,1) = 0.3
c      xgamma(ind_ClNO3,2) = 0.3
c      xgamma(ind_HOBr,1)  = 0.2
c      xgamma(ind_HOBr,2)  = 0.2
c      xgamma(ind_HOCl,1)  = 0.
c      xgamma(ind_HOCl,2)  = 0.
c      xgamma(ind_HOI,1)   = 0.
c      xgamma(ind_HOI,2)   = 0.

c the following are needed to calculate dry rates w/ assumption of Henry's
c law equilibrium

c define equilibrium constant
c obviously without Pitzer coefficients
      do k=2,nmax
         xeq(k,ind_HNO3) = funa(1.54d+10,8700.d0,k)
      enddo

c define inverse Henry's constant 
      do k=2,nmax
         henry(k,ind_HNO3)=func3(2.5d6/1.5d1,8694.d0,k) !RS_MOCCA_exakt
         FCT=0.0820577*tt(k)
         do l=1,ndr
            if (henry(k,idr(l)).ne.0.d0) henry(k,idr(l))=1./ ! jjb .d0
     &           (henry(k,idr(l))*FCT)
c "else": henry=0 <=> k_H^cc=infinity
         enddo
      enddo

c mean molecular speed (see SR v_mean_* for details)
c v_mean in m/s; sqrt(8*R_gas/pi)=4.60138
      do k=2,nmax
         vmean(ind_HNO3,k)  = func(6.3d-2,k)
         vmean(ind_N2O5,k)  = func(1.08d-1,k)
         vmean(ind_NH3,k)   = func(1.7d-2,k)
         vmean(ind_H2SO4,k) = func(9.8d-2,k)
!         vmean(ind_HCl,k)   = func(3.6d-2,k)
c         vmean(ind_BrNO3,k) = func(1.42d-1,k)
c         vmean(ind_ClNO3,k) = func(9.7d-2,k)
c         vmean(ind_HOBr,k)  = func(9.7d-2,k)
c         vmean(ind_HOCl,k)  = func(5.2d-2,k)
c         vmean(ind_HOI,k)   = func(1.44d-1,k)
      enddo

c calculate kmt ----

c k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha)) 
c if : D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
c k_mt=v_mean/(r**2/lambda+4*r/(3*alpha))

c lambda=freep : mean free path
!     x1=0. ! jjb useless
      do k=2,nmax
         do l=1,ndr
            do kc=1,2
               if (xgamma(idr(l),kc).gt.0..and.rcd(kc,k).gt.0.) then
                  x1=1./(rcd(kc,k)*(rcd(kc,k)/freep(k)+4./(3.*xgamma
     &                 (idr(l),kc)))) 
               else
                  x1=0.
               endif
               xkmtd(k,kc,idr(l))=vmean(idr(l),k)*x1
!                  print *,k,idr(l),kc
!                  print *,xkmtd(k,kc,idr(l)),x1,vmean(idr(l),k)
            enddo
         enddo
      enddo

c      do k=2,n
c            write (432, 1001) k,xkmtd(k,1,ind_HNO3),1./xkmtd(k,1,ind_HNO3),
c     &        xkmtd(k,2,ind_HNO3),1./xkmtd(k,2,ind_HNO3),
c     &        xkmtd(k,1,ind_HNO3)*cwd(k,1),xkmtd(k,2,ind_HNO3)*cwd(k,2)
c      enddo
! 1001 format(i4, 6d16.8)

      end subroutine dry_rates_g


c
c------------------------------------------------------------
c

      subroutine dry_rates_a (freep,nmaxf)
c calculates kmt for heterogeneous reactions on dry aerosol
c 1: sulfate aerosol, 2: seasalt aerosol

c this SR calculates the gammas for the case where one aerosol class is
c allready above its crystallization point, but the other one is not yet
c so for the first class the complete aerosol chemistry is active, for
c the second only the reactions on dry aerosol (regulated via xliq1/2 and
c xhet1/2 in SR kpp_driver)

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h'     !additional common blocks and other definitions
     
      parameter (ndr=4)
      common /blck11/ rcd(nkc,n)
!     common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),      ! jjb  none of the objects of the common block is used
!    &              e(nkt),dew(nkt),rq(nkt,nka)         !   (after commenting rqm = rq line below)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_drya/ xkmtd(nf,2,NSPEC),xeq(nf,NSPEC)
!      common /kpp_dryp/ rcd(n,2),cwd(n,2)
!     dimension xgamma(NSPEC,2),freep(nf),idr(ndr),rqm(nkt,nka) ! jjb rqm is unreferenced
      dimension xgamma(NSPEC,2),freep(nf),idr(ndr)              ! jjb thus removed


c      data idr/ind_HNO3,ind_N2O5,ind_BrNO3,ind_ClNO3,ind_HOBr/
!     data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4,ind_HCl/     ! jjb HCl is not handled here (but it should be, probably !)
      data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4/           ! jjb removed at the moment

!     funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3)) ! jjb another copy-paste mistake from SR dry_rates_g
      funa(a0,b0,k)=a0*exp(b0*(1/t(k)-3.354d-3))  ! jjb here the temperature is passed through a CB, not as a parameter

c change rq in um to rqm in m            ! jjb rqm is not used in this SR
!     do ia=1,nka                        !   ( see commented block at the end)
!        do jt=1,nkt
!           rqm(jt,ia)=rq(jt,ia)*1.d-6
!        enddo
!     enddo

c define gamma's for all species ----
      xgamma(ind_HNO3,1)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_HNO3,2)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_N2O5,1)  = 0.02       !estimated from JPL 2003
      xgamma(ind_N2O5,2)  = 0.02       !estimated from JPL 2003
      xgamma(ind_NH3,1)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_NH3,2)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_H2SO4,1) = 0.1        !estimated from alpha=0.65
      xgamma(ind_H2SO4,2) = 0.1        !estimated from alpha=0.65
c      xgamma(ind_BrNO3,1)=0.3
c      xgamma(ind_BrNO3,2)=0.3
c      xgamma(ind_ClNO3,1)=0.3
c      xgamma(ind_ClNO3,2)=0.3
c      xgamma(ind_HOBr,1)=0.2
c      xgamma(ind_HOBr,2)=0.2
c      xgamma(ind_HOCl,1)=0.
c      xgamma(ind_HOCl,2)=0.
c      xgamma(ind_HOI,1)=0.
c      xgamma(ind_HOI,2)=0.


c the following are needed to calculate dry rates w/ assumption of Henry's
c law equilibrium

c define equilibrium constant
c obviously without Pitzer coefficients
!     do k=2,nmax  ! jjb nmax is not defined (copy-paste mistake from SR dry_rates_g)
      do k=2,nmaxf ! jjb nmaxf is the correct index here
         xeq(k,ind_HNO3) = funa(1.54d+10,8700.d0,k)
      enddo

c calculate kmt ----

c k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha)) 
c if : D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
c k_mt=v_mean/(r**2/lambda+4*r/(3*alpha))

c lambda=freep : mean free path
!     x1=0. ! jjb useless
      do k=2,nmaxf
         do l=1,ndr
            do kc=1,2
               if (xgamma(idr(l),kc).gt.0..and.rcd(kc,k).gt.0.) then
                  x1=1./(rcd(kc,k)*(rcd(kc,k)/freep(k)+4./(3.*xgamma
     &                 (idr(l),kc)))) 
               else
                  x1=0.
               endif
               xkmtd(k,kc,idr(l))=vmean(idr(l),k)*x1
!                  print *,k,idr(l),kc
!                  print *,xkmtd(k,kc,idr(l)),x1,vmean(idr(l),k)
            enddo
         enddo
      enddo

c      write (434,1100)
c      do k=2,n
c            write (434, 1101) k,xkmtd(k,1,ind_HNO3),1./xkmtd(k,1,ind_HNO3),
c     &        xkmtd(k,2,ind_HNO3),1./xkmtd(k,2,ind_HNO3),
c     &        xkmtd(k,1,ind_HNO3)*cwd(k,1),xkmtd(k,2,ind_HNO3)*cwd(k,2)
c      enddo
! 1100 format ('bulk values')
! 1101 format(i4, 6d16.8)
! 1102 format ('integrated values')


c now the same but with integrated values - - - - - - - - - - - - - - - - 
c if working: include smart way of saving CPU time
c this leads to kmt that are about 10x greater than the bulk approach,
c nevertheless use bulk as this is what gamma's are determined for, plus
c it's faster (CPU-wise) 


cc loop over vertical grid
c      do  k=2,nmaxf
cc loop over the nkc different chemical bins
c         do kc=1,2!nkc
cc loop over the species to be exchanged between gas and aqueous phase---
c            do l=1,ndr
c               xkmtd(k,kc,idr(l))=0.
cc define summation limits (1) ---
c               if (kc.eq.1) then
c                  iia_0=1
c                  iia_e=ka
c               endif
c               if (kc.eq.2) then
c                  iia_0=ka+1
c                  iia_e=nka
c               endif
c               ! kc eq. 3/4
cc fast version without logarithmic integration
c               x1=0.
c               xk1=0.
c               if (l.eq.1) xx1=0.
c               if (alpha(k,idr(l)).gt.0.) x1=4./(3.*alpha(k,idr(l)))
c               do ia=iia_0,iia_e
cc define summation limits (2)
c                  if (kc.eq.1) then
c                     jjt_0=1
c                     jjt_e=kw(ia)
c                  endif
c                  if (kc.eq.2) then
c                     jjt_0=1
c                     jjt_e=kw(ia)
c                  endif
c                  ! kc eq. 3/4
c                  do jt=jjt_0,jjt_e
cc conversion: um      --> m               :10^-6
c                     rqq=rqm(jt,ia)
cc kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
cc     with D_g=lambda*vmean/3.
cc here a volume weighted value is calculated, therefore weighting with r^3:
cc kmt=4/3*pi/L*sum(a)*sum(r){r^3*N*kt}
c                     x2=vmean(idr(l),k)/(rqq/freep(k)+x1) ![1/s]
cc conversion: 1/cm^3 --> 1/m^3(air):10^6
c                     xk1=xk1+x2*rqq*rqq*ff(jt,ia,k)*1.d6
c                  enddo !jt
c               enddo !ia
cc k_mt=4*pi/(3*LWC)*sum
c!               xkmtd(k,kc,idr(l))=4.*3.1415927/(3.*cwd(k,kc))*xk1 ![1/s]
c               xkmtd(k,kc,idr(l))=4.*pi/(3.*cwd(k,kc))*xk1 ![1/s]
c            enddo !l 
c         enddo                  ! kc
c      enddo                     !k

c      write (434,1102)
c      do k=2,n
c            write (434, 1101) k,xkmtd(k,1,ind_HNO3),1./xkmtd(k,1,ind_HNO3),
c     &        xkmtd(k,2,ind_HNO3),1./xkmtd(k,2,ind_HNO3),
c     &        xkmtd(k,1,ind_HNO3)*cwd(k,1),xkmtd(k,2,ind_HNO3)*cwd(k,2)
c      enddo

      end subroutine dry_rates_a

c
c------------------------------------------------------------
c

      subroutine dry_rates_t (freep,nmaxf)
c calculates kmt for heterogeneous reactions on dry aerosol
c 1: sulfate aerosol, 2: seasalt aerosol

c this SR calculates the gammas for the case where one aerosol class is
c allready above its crystallization point, but the other one is not yet
c so for the first class the complete aerosol chemistry is active, for
c the second only the reactions on dry aerosol (regulated via xliq1/2 and
c xhet1/2 in SR kpp_driver)

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h'     !additional common blocks and other definitions
     
      parameter (ndr=4)
      common /blck11/ rcd(nkc,n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_dryt/ xkmtd(nf,2,NSPEC),xeq(nf,NSPEC)
!     common /kpp_dryp/ rcd(n,2),cwd(n,2)
!     dimension xgamma(NSPEC,2),freep(nf),idr(ndr),rqm(nkt,nka) ! jjb rqm is unreferenced
      dimension xgamma(NSPEC,2),freep(nf),idr(ndr)              ! jjb thus removed


c      data idr/ind_HNO3,ind_N2O5,ind_BrNO3,ind_ClNO3,ind_HOBr/
!     data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4,ind_HCl/     ! jjb HCl is not handled here (but it should be, probably !)
      data idr/ind_HNO3,ind_N2O5,ind_NH3,ind_H2SO4/           ! jjb removed at the moment

!     funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3)) ! jjb another copy-paste mistake from SR dry_rates_g
      funa(a0,b0,k)=a0*exp(b0*(1/t(k)-3.354d-3))  ! jjb here the temperature is passed through a CB, not as a parameter

c change rq in um to rqm in m            ! jjb rqm is not used in this SR
!     do ia=1,nka                        !   ( see commented block at the end)
!        do jt=1,nkt
!           rqm(jt,ia)=rq(jt,ia)*1.d-6
!        enddo
!     enddo

c define gamma's for all species ----
      xgamma(ind_HNO3,1)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_HNO3,2)  = 0.02       !JPL 2003: uptake on halide salts
      xgamma(ind_N2O5,1)  = 0.02       !estimated from JPL 2003
      xgamma(ind_N2O5,2)  = 0.02       !estimated from JPL 2003
      xgamma(ind_NH3,1)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_NH3,2)   = 0.05       !Dentener and Crutzen, 1994, #744
      xgamma(ind_H2SO4,1) = 0.1        !estimated from alpha=0.65
      xgamma(ind_H2SO4,2) = 0.1        !estimated from alpha=0.65
c      xgamma(ind_BrNO3,1)=0.3
c      xgamma(ind_BrNO3,2)=0.3
c      xgamma(ind_ClNO3,1)=0.3
c      xgamma(ind_ClNO3,2)=0.3
c      xgamma(ind_HOBr,1)=0.2
c      xgamma(ind_HOBr,2)=0.2
c      xgamma(ind_HOCl,1)=0.
c      xgamma(ind_HOCl,2)=0.
c      xgamma(ind_HOI,1)=0.
c      xgamma(ind_HOI,2)=0.


c the following are needed to calculate dry rates w/ assumption of Henry's
c law equilibrium

c define equilibrium constant
c obviously without Pitzer coefficients
!     do k=2,nmax  ! jjb nmax is not defined (copy-paste mistake from SR dry_rates_g)
      do k=2,nmaxf ! jjb nmaxf is the correct index here
         xeq(k,ind_HNO3) = funa(1.54d+10,8700.d0,k)
      enddo

c calculate kmt ----

c k_mt=1/(r**2/(3*D_gas)+4*r/(3*v_mean*alpha)) 
c if : D_gas=lambda*v_mean/3. then we can rewrite k_mt as:
c k_mt=v_mean/(r**2/lambda+4*r/(3*alpha))

c lambda=freep : mean free path
!     x1=0. ! jjb useless
      do k=2,nmaxf
         do l=1,ndr
            do kc=1,2
               if (xgamma(idr(l),kc).gt.0..and.rcd(kc,k).gt.0.) then
                  x1=1./(rcd(kc,k)*(rcd(kc,k)/freep(k)+4./(3.*xgamma
     &                 (idr(l),kc)))) 
               else
                  x1=0.
               endif
               xkmtd(k,kc,idr(l))=vmean(idr(l),k)*x1
            enddo
         enddo
      enddo

      end subroutine dry_rates_t

c
c------------------------------------------------------------
c

      subroutine activ (box,n_bl)
! front end for Beiping Luo's model to calculate activity coefficients

! jjb work done:
!     cleaned, declared all variables, implicit none, removed old, commented stuff
!     cleaned final write (+ avoid useless write when cm < cm_min)
!     use parameter for min cm value
!     (probable) bugfix to be done: loop kc = 1,nkc: use nkc_l instead of nkc, if nkc_l is also
!               used (future development) to restrict tot mechanism
!     use xip (not xit) to check is pitzer will be called or not

! 03-Mar-2017  <J. Bock> When Pitzer is not called, set gamma=1 (default value) to avoid previous
!                        values to be still used

! 04-Mar-2017  <J. Bock> Introduce "xip2" to check the lower ionic strengh value allowed in Piter.
!                        See SR gammann: 1./xo4 where xo4=(omega1*I2)**4 and I2=sqrt(I) (I=xip here)
!                        min(omega(i,j)) = 0.43 thus xip2 = 0.01*xip**2

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nkc

      implicit none

      logical, intent(in) :: box
      integer, intent(in) :: n_bl

! Local parameters:
      integer, parameter :: nk = 10  ! Number of vertical levels outputed
      integer, parameter :: kk(nk) = (/2,12,22,32,42,52,62,72,82,92/) ! Indexes of levels in output

! Local scalars:
      integer :: jg, k, kc   ! loop indexes for liq. species, vertical layers, liq. bins
      integer :: nmin, nmax  ! lower and upper bounds for calculations

      double precision :: wact
      double precision :: xip2

! Local arrays:
      double precision :: wa(nkc,nf)
      double precision :: xip(nkc,nf), xit(nkc,nf)

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      double precision xm1, xm2, feu, dfddt, xm1a, xm2a

      common /blck12/ cw(nkc,n),cm(nkc,n)
      double precision cw, cm

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /kpp_mol/ xgamma(nf,j6,nkc)
      double precision xgamma


! Initialise lower and upper bounds (in the vertical grid) to do the calculations
      nmin=2
      nmax=nf
      if (box) then
         nmin=n_bl
         nmax=n_bl
      endif

! Initialise gammas
      xgamma(:,:,:) = 1.d0

      do k=nmin,nmax
         do kc=1,nkc

            ! Skip if cm too small
            if (cm(kc,k) <= tiny(0.d0)) cycle

! calculate ionic strength of species accounted for in pitzer module
!       ionic strength: I=0.5*sum(molality*charge^2)
! concentrations --> molality:
!  mol/m^3_air * m^3_air/m^3_solvent * 10^-3 * 1 dm^3/kg (density of water) = mol/kg_solvent
!   sion1      * cm^-1               * 1.d-3
            xip(kc,k)=0.5*(sion1(1,kc,k)+sion1(2,kc,k)+sion1(20,kc,k)+
     &           sion1(19,kc,k)+4.*sion1(8,kc,k)+sion1(13,kc,k)+
     &           sion1(14,kc,k))
            xip(kc,k)=xip(kc,k)*1.d-3/cm(kc,k)
! calculate ionic strength of all species
            xit(kc,k)=0.5*(sion1(1,kc,k)+sion1(2,kc,k)+sion1(3,kc,k)+
     &           sion1(4,kc,k)+sion1(5,kc,k)+4.*sion1(6,kc,k)+
     &           sion1(7,kc,k)+4.*sion1(8,kc,k)+sion1(9,kc,k)+
     &           sion1(10,kc,k)+sion1(11,kc,k)+sion1(12,kc,k)+
     &           sion1(13,kc,k)+sion1(14,kc,k)+sion1(15,kc,k)+
     &           sion1(16,kc,k)+sion1(19,kc,k)+sion1(20,kc,k)+
     &           sion1(21,kc,k)+sion1(22,kc,k)+sion1(23,kc,k)+
     &           sion1(24,kc,k)+sion1(25,kc,k)+sion1(26,kc,k)+
     &           sion1(27,kc,k)+sion1(28,kc,k)+sion1(29,kc,k)+
     &           sion1(30,kc,k)+sion1(31,kc,k)+sion1(32,kc,k)+
     &           sion1(33,kc,k)+sion1(34,kc,k)+sion1(35,kc,k)+
     &           sion1(36,kc,k)+sion1(37,kc,k)+sion1(38,kc,k)+
     &           sion1(39,kc,k))
     &           *1.d-3/cm(kc,k)


! check that ionic strength remains in a "reasonable" range; if too high there
! will also be a "floating overflow" in SR pitzer
            xip2 = 0.01*xip(kc,k)**2 ! <jjb> see explanations in header. 0.01 < 0.43**4
            !if (xit(kc,k).le.80..and.xit(kc,k).gt.0.) then
            if (xip(kc,k).le.80..and.xip2.gt.tiny(0.)) then
! calculate pitzer coefficients for layer k and liquid size bin kc
               call pitzer (k,kc,wact) ! for sulfate data can be used for ionic strengths
                                       ! up to 40 M, for sea salt only up to 6 M
               wa(kc,k)=wact
            else
               print *,'ionic strength > 80. or <=0.',k,kc
               print *,'I=',xip(kc,k),xit(kc,k)
               cycle
            endif

! gamma's calculated by SR pitzer are for molalities, used here are molarities ==> conversion factor needed:
! mol/kg(solvent) --> mol/l(solution): cm/cw (assuming unity density for the solvent water)


! don't apply to all xgamma's only to those that are <> 1:
            if (cw(kc,k).gt.0.d0) then
               do jg=1,j6
                  if (xgamma(k,jg,kc).ne.1.) xgamma(k,jg,kc) =
     &                 xgamma(k,jg,kc) * cm(kc,k)/cw(kc,k)
               enddo
            end if

! define gamma's for species that are not included in pitzer module
! L+J: Liang and Jacobson, 1999, JGR, 104, 13749, #554
! C+S: Chameides and Stelson, 1992, JGR, 97,20565,#470
            xgamma(k,3,kc)=xgamma(k,13,kc) !OH-   = NO3- (assumed, Luo, pers comm 2000)
            xgamma(k,5,kc)=xgamma(k,19,kc) !HSO3- = HSO4- (L+J)
            xgamma(k,6,kc)=xgamma(k,8,kc)  !SO3=  = SO4=  (L+J)
            xgamma(k,7,kc)=xgamma(k,19,kc) !SO4-  = HSO4- (L+J)
            xgamma(k,9,kc)=xgamma(k,5,kc)  !HCO3- = HSO3- (C+S)
            xgamma(k,11,kc)=xgamma(k,5,kc) !O2-   = Cl2- = HSO3- (C+S)
            xgamma(k,12,kc)=xgamma(k,13,kc)!NO2-  = NO3-  (L+J)
            xgamma(k,15,kc)=xgamma(k,5,kc) !Cl2-  = HSO3- (C+S)
            xgamma(k,16,kc)=xgamma(k,5,kc)  !HCOO- = HSO3- (assumed)????
!            xgamma(k,21,kc)=xgamma(k,,kc) !NO4-  = ??   (assumed)
            xgamma(k,22,kc)=xgamma(k,14,kc)!ClO-  = Cl-   (assumed)
            xgamma(k,24,kc)=xgamma(k,14,kc)!Br-   = Cl-   (assumed)
            xgamma(k,25,kc)=xgamma(k,5,kc) !Br2-  = HSO3- (C+S)
            xgamma(k,26,kc)=xgamma(k,24,kc)!BrO-  = Br-   (assumed)
!            xgamma(k,28,kc)=xgamma(k,,kc)  !BrCl2-= ? not used in SR equil_co*
!            xgamma(k,29,kc)=xgamma(k,,kc)  !Br2Cl-= ?       -"-
            xgamma(k,37,kc)=xgamma(k,5,kc) !ICl2- = HSO3- (assumed)
            xgamma(k,38,kc)=xgamma(k,5,kc) !IBr2- = HSO3- (assumed)

         enddo  !kc
      enddo     !k

! output
      if (lmin/30*30.eq.lmin) then
 3000    continue
         open (109,file='gam.out',status='unknown',position='append',
     &        err=3000)
         do k=1,nk
            do kc=1,nkc
               write (109,20) kk(k),kc,lday,lst,lmin,feu(kk(k)),
     &                        cm(kc,kk(k)),cw(kc,kk(k)),xm2(kk(k))*1.d-3

               ! Skip next 'write' instructions if cm too small
               if (cm(kc,kk(k)) <= tiny(0.d0)) then
                  write(109,*)' -> cm too small, xgamma set equal to 1.'
                  cycle
               end if

               write (109,13) xip(kc,kk(k)),xit(kc,kk(k)),wa(kc,kk(k))
               write (109,10) (xgamma(kk(k),jg,kc),jg=1,j6)
               write (109,10) (sion1(jg,kc,kk(k))*1.d-3/cm(kc,kk(k)),
     &                         jg=1,j6)
            enddo
         enddo
         close (109)
      endif
 10   format (8d14.6)
 13   format (3d14.6)
 20   format ('layer ',i3,' bin ',i1,' day ',2i3,':',i2,' rH ',d12.4,
     &        ' LWC ',3d12.4)

      end subroutine activ

c
c------------------------------------------------------------
c

      subroutine activ_init
! initialise activity coefficients (gammas)

      USE global_params, ONLY :
! Imported Parameters:
     &     j6,
     &     nf,
     &     nkc

      implicit none

      common /kpp_mol/ xgamma(nf,j6,nkc) 
      double precision xgamma

      xgamma(:,:,:) = 1.d0

      open (109,file='gam.out',status='unknown',form='formatted')
      write(109,*)'File opened by SR activ_init, written by SR activ'
      write(109,*)'  For each outputed layer:'
      write(109,*)'  partial ionic strength, total ionic str, water act'
      write(109,*)'  gamma (all ionic species)'
      write(109,*)'  molalities (all ionic species)'
      write(109,*)'  --------------------------------------------------'
      close (109)

      end subroutine activ_init

c
c-------------------------------------------------------------
c

      subroutine gasdrydep (xra,tt,rho,freep)
c calculate dry deposition velocities for gas phase after Sehmel cited in
c Seinfeld and Pandis, 1999


! jjb IMPORTANT NOTE
!
! Some points need to be addressed in this subroutine.

! 1) use mk tables to convert mistra <--> KPP indexes
! 2) change the if test in the middle (after hs(i) and "FCT" definitions) : if hs .gt. 0 (instead of .ne.)
!    and replace the test at the end : if hs == -1 instead of if hs == -1/FCT
! 3) correction of one formula, double check


      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nkc

      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
! Imported Array Variables with intent(in):
     &     ind_gas,
     &     vmean,
     &     vg

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h'     !additional common blocks and other definitions

      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
!      common /gas_vdd/ vg(j1)
      common /cb46/ ustern,gclu,gclt
!     dimension tt(n),freep(nf),rho(n),rb(j1),rc(j1),vm(j1),hs(j1), ! jjb rb & rc not used
!    &     f0(j1)
      dimension tt(n),freep(nf),rho(n),vm(j1),hs(ind_gas(j1)),
     &          f0(ind_gas(j1))
c function used in calculation of Hstar
      funa(a0,b0,k)=a0*exp(b0*(1/tt(k)-3.354d-3))

c gas phase dry deposition velocity:v_d=1/(ra + rb + rc)
c    ra=1/(kappa ustar) (ln (z/z0) +Phi) ;where z=height of surface (constant flux) layer
c                                         Phi takes stratification into account
c      see SR partdep, ra=xra
c    rb=5.*Sc^(2/3)/ustar ;Sc=nu/D, nu=kinematic viscosity of air, D diffusivity of gas
c      ustar=ustern
c      gas phase diffusivity is calculated for each species using the mean molecular speed vmean:
c        D=lambda*vmean/3.  (lambda=freep)
c rc after Sehmel, 1980
c    rc=2.54d+4/(Hstar * T * ustar) ; Hstar=effective Henry constant
c rc after Wesely, 1989
c     rc=1./(H*/(10^5*r_gsS)+f_0/r_gsO)  ; r_gsS=1., r_gsO=2000., f_0: standard 0.1
c                                                                      low reactivity:0., high:1.

c      Hstar calculated from henry and sea water pH=8.1 (Riley and Skirrow, 1965)

      k=2 ! only in lowest model layer
      xeta=1.8325d-5*(416.16/(tt(k)+120.))*((tt(k)/296.16)**1.5) !dynamic viscosity of air, Jacobsen p. 92
      xnu=xeta/rho(k)  !kinematic viscosity of air

c get vmean for j1-list (no deposition or radicals or fixed)
      !do jspec = 1,j1
         vm(:) = vmean(:j1,k)
      !end do

c     vm(radical 29)=vmean(ind_OIO,k)
c get henry constant for j1-order from KPP names and calculate Hstar
      sac=10.**(-8.1d0) ![H+]=10^(- pH), pH=8.1
      hs(1)=henry(ind_NO,k)
      hs(2)=henry(ind_NO2,k)
      hs(3)=henry(ind_HNO3,k)
      hs(4)=henry(ind_NH3,k)
      hs(5)=henry(ind_SO2,k)
      hs(6)=henry(ind_H2SO4,k)
      hs(7)=henry(ind_O3,k)
      hs(8)=henry(ind_CH4,k)
      hs(9)=henry(ind_C2H6,k)
c      hs(10)=henry(ind_C3H8,k)
c      hs(11)=henry(ind_ALKA,k)
      hs(12)=henry(ind_ETHE,k)
c      hs(13)=henry(ind_ALKE,k)
c      hs(14)=henry(ind_AROM,k)
      hs(15)=henry(ind_ACO2,k)
      hs(16)=henry(ind_ACTA,k)
      hs(17)=henry(ind_HCHO,k)
      hs(18)=henry(ind_ALD2,k)
      hs(19)=henry(ind_H2O2,k)
      hs(20)=henry(ind_ROOH,k)
      hs(21)=henry(ind_HONO,k)
      hs(22)=henry(ind_PAN,k)
c      hs(23)=henry(ind_TPAN,k)
c      hs(24)=henry(ind_KET,k)
c      hs(25)=henry(ind_CRES,k)
c      hs(26)=henry(ind_DIAL,k)
c      hs(27)=henry(ind_GLYX,k)
c      hs(28)=henry(ind_MGLY,k)
c      hs(29)=henry(ind_NH4NO3,k)
      hs(30)=henry(ind_HCl,k)
c      hs(31)=henry(ind_R3N2,k)
c      hs(32)=henry(ind_RAN2,k)
c      hs(33)=henry(ind_RAN1,k)
      hs(34)=-1.!henry(ind_N2O5,k)
      hs(35)=henry(ind_HNO4,k)
      hs(36)=henry(ind_NO3,k)
      hs(37)=henry(ind_DMS,k)
      hs(38)=henry(ind_HOCl,k)
      hs(39)=henry(ind_ClNO2,k)
      hs(40)=-1.!henry(ind_ClNO3,k)
      hs(41)=henry(ind_Cl2,k)
      hs(42)=henry(ind_HBr,k)
      hs(43)=henry(ind_HOBr,k)
      hs(44)=henry(ind_BrNO2,k)
      hs(45)=-1.!henry(ind_BrNO3,k)
      hs(46)=henry(ind_Br2,k)
      hs(47)=henry(ind_BrCl,k)
      hs(48)=-1.!henry(ind_HI,k)!*(1. + /sac)
      hs(49)=henry(ind_HOI,k)!*(1. + /sac)
      hs(50)=henry(ind_I2O2,k)
      hs(51)=henry(ind_INO2,k)
      hs(52)=-1.!henry(ind_INO3,k)
      hs(53)=henry(ind_I2,k)
      hs(54)=henry(ind_ICl,k) 
      hs(55)=henry(ind_IBr,k)
      hs(56)=henry(ind_CH3I,k)
      hs(57)=henry(ind_CH2I2,k)
      hs(58)=henry(ind_CH2ClI,k)
      hs(59)=henry(ind_C3H7I,k)
      hs(60)=henry(ind_DMSO,k)
      hs(61)=henry(ind_CH3SO2,k)
      hs(62)=henry(ind_CH3SO3,k)
      hs(63)=henry(ind_CH3SO3H,k)! MSA=CH3S(OO)OH
      hs(64)=henry(ind_CO,k)
      hs(65)=henry(ind_Cl2O2,k)
      hs(66)=henry(ind_DMOO,k) ! CH3SCH2OO
      hs(67)=henry(ind_CH3S,k)
      hs(68)=henry(ind_CH3SO,k)
      hs(69)=henry(ind_CH3SO2H,k) ! MSIA=CH3S(O)OH
      hs(70)=henry(ind_DMSO2,k)
      hs(71)=henry(ind_CH2BrI,k)
!      hs(72)=henry(ind_CHBr2I,k)
      hs(73)=henry(ind_C2H5I,k)
      hs(74)=henry(ind_HIO3,k)
c      hs(75)=henry(ind_NUCV,k)
      hs(76)=henry(ind_SO3,k)
      hs(77)=henry(ind_HOSO2,k)
      hs(78)=henry(ind_CO2,k)
!      hs(79)=henry(ind_I2O,k)
!      hs(80)=henry(ind_I2O3,k)
!      hs(81)=henry(ind_I2O4,k)
!      hs(82)=henry(ind_I2O5,k)
!      hs(83)=henry(ind_INO,k)
      hs(84)=henry(ind_Br2O,k)
      hs(85)=henry(ind_ClONO,k)
      hs(86)=henry(ind_ClO3,k)
      hs(87)=henry(ind_Cl2O3,k)
      hs(88)=henry(ind_CH3OH,k)
      hs(89)=henry(ind_C2H5OH,k)
      hs(90)=henry(ind_H2,k)
      hs(91)=henry(ind_NHS,k)
      hs(92)=henry(ind_RCl,k)
      hs(93)=henry(ind_RBr,k)
      hs(94)=henry(ind_XOR,k)
      hs(95)=henry(ind_SOR,k)
      hs(96)=henry(ind_SPAN,k)
c      hs(97)=henry(ind_Hg,k)
c      hs(98)=henry(ind_HgO,k)
c      hs(99)=henry(ind_HgCl,k)
c      hs(100)=henry(ind_HgCl2,k)
c      hs(101)=henry(ind_HgBr,k)
c      hs(102)=henry(ind_HgBr2,k)

c henry constant were already transformed into dimensionless values
      FCT=0.0820577*tt(k)
      do i=1,j1
         if (hs(ind_gas(i)).ne.0.)hs(ind_gas(i))=1./(hs(ind_gas(i))*FCT)
         f0(ind_gas(i))=0.1
      enddo
c some species dissociate:
      hs(3)=hs(3)*(1. +funa(1.54d+1,8700.d0,k)/sac)
      hs(4)=hs(4)*(1. + funa(1.7d-5,-4325.d0,k)* 
     &     sac/funa(1.d-14,-6710.d0,k)) 
      hs(5)=hs(5)*(1.+ funa(1.7d-2,2090.d0,k)/sac + 
     &     funa(1.7d-2,2090.d0,k)*funa(6.0d-8,1120.d0,k) /sac**2)
      hs(6)=hs(6)*(1. + 1.0d+3/sac +
     &     1.0d+3*funa(1.02d-2,2720.d0,k)/sac**2) 
      hs(30)=hs(30)*(1. + funa(1.7d6,6896.d0,k)/sac) 
      hs(38)=hs(38)*(1. + 3.2d-8/sac)
      hs(42)=hs(42)*(1. + 1.d9/sac)
      hs(43)=hs(43)*(1. + funa(2.3d-9,-3091.d0,k)/sac)
c define some f_0 that deviate from standard (Pandis and Seinfeld, Table 19.3):
      f0(1) =0.  !NO
      f0(3) =0.  !HNO3
      f0(4) =0.  !NH3
      f0(5) =0.  !SO2
      f0(7) =1.  !O3
      f0(8) =0.  !CH4
      f0(9) =0.  !C2H6
      f0(10)=0.  !C3H8
      f0(11)=0.  !ALKA
      f0(14)=0.  !AROM
      f0(15)=0.  !HCOOH, formic acid
      f0(16)=0.  !CH3COOH, acetic acid
      f0(17)=0.  !HCHO
      f0(19)=1.  !H2O2
      f0(20)=1.  !ROOH
      f0(30)=0.  !HCl
      f0(35)=0.  !HNO4
      f0(36)=1.  !NO3
      f0(42)=0.  !HBr
c      f0(43)=1. !HOBr - didn't find a reference for that, so I commented it again
c      f0(49)=1. !HOI
c      f0() =0.  !CH3CHO, acetaldehyde
c      f0() =0.  !CH3CH2CHO, propionaldehyde
c      f0() =0.3 !CH3OOH, methylhydroperoxide

c calculate vg from rb and rc for each species
!     rb_fact=5./ustern*(xnu*freep(k)/3.)**(2./3.) ! jjb mistake ?
      rb_fact=5./ustern*(3.*xnu/freep(k))**(2./3.)
c      rc_fact=2.54d+4/(tt(k)*ustern)  !Sehmel, 1980
c      print *,rb_fact,rc_fact
      do i=1,j1
       if (vm(i).eq.0.) then 
         vg(i)=0.
       else
         if (hs(ind_gas(i)).ne.0.) then
c           vg(i)=1./(xra+(rb_fact/(vm(i)**(2./3.)))+(rc_fact/hs(i))) !Sehmel, 1980
            vg(i)=1./(xra+(rb_fact/(vm(i)**(2./3.)))+1./(hs(ind_gas(i))*
     &           1.d-5+f0(ind_gas(i))/2000.))  !Wesely, 1989
         else
c            vg(i)= 1./(xra+(rb_fact/(vm(i)**(2./3.)))+(rc_fact/1.d-7)) !hs=0 => 1/hs -> infinity, 
                                                              !here just a value for hstar is chosen
                                                              !to get a small v_dd
            if (f0(ind_gas(i)) > 0.) then
!               print*,xra,rb_fact,vm(i),f0(ind_gas(i))
               vg(i)=1./(xra+(rb_fact/(vm(i)**(2./3.)))
     &               +1./(f0(ind_gas(i))/2000.))
            else
               vg(i)=0.
            end if
         endif
         if (hs(ind_gas(i)).eq.(-1./FCT)) 
     &        vg(i)=1./(xra+(rb_fact/(vm(i)**(2./3.)))+.1)     ! "infinite solubility"
      endif
      enddo

      end subroutine gasdrydep


c
c-----------------------------------------------------------------------------
c

      subroutine mass_ch
c calculates total number of Br and Cl atoms in all phases [mol/m2]
c including deposited atoms
c emitted atoms (see SR aer_source) are subtracted

      USE config, ONLY:
     &     nkc_l

      USE gas_common, ONLY:
     &     s1,
     &     j1_br, ind_gas_br,
     &     j1_cl, ind_gas_cl,
     &     s3,
     &     j5_br, ind_rad_br,
     &     j5_cl, ind_rad_cl

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j3,
     &     j6,
     &     n,
     &     nkc


      implicit double precision (a-h,o-z) 
      logical cl_out

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /sss/ brsss,clsss,xnasss

      if ((lday.eq.0.and.lst.eq.0.and.lmin.le.1).or.lmin/30*30.eq.lmin)
     &     then       
      cl_out=.false.
      if (nkc_l.gt.2) cl_out=.true.
c gas phase----------------------------
c atmosphere
      brg=0.
      clg=0.
      brgd=0.
      clgd=0.
      do k=2,n
         ! Brominated species
         brgp=0.
         do j=1,j1_br
            brgp=brgp+ind_gas_br(1,j)*s1(ind_gas_br(2,j),k)
         end do
         do j=1,j5_br
            brgp=brgp+ind_rad_br(1,j)*s3(ind_rad_br(2,j),k)
         end do
         brg=brg+brgp*detw(k) ! to get mol/m^2

         ! ditto for chlorinated species
         clgp=0.
         do j=1,j1_cl
            clgp=clgp+ind_gas_cl(1,j)*s1(ind_gas_cl(2,j),k)
         end do
         do j=1,j5_cl
            clgp=clgp+ind_rad_cl(1,j)*s3(ind_rad_cl(2,j),k)
         end do
         clg=clg+clgp*detw(k)
      enddo

c  + deposited ! already in mol/m^2
      do j=1,j1_br
         brgd=brgd+ind_gas_br(1,j)*s1(ind_gas_br(2,j),1)
      end do
      do j=1,j5_br
         brgd=brgd+ind_rad_br(1,j)*s3(ind_rad_br(2,j),1)
      end do

      do j=1,j1_cl
         clgd=clgd+ind_gas_cl(1,j)*s1(ind_gas_cl(2,j),1)
      end do
      do j=1,j5_cl
         clgd=clgd+ind_rad_cl(1,j)*s3(ind_rad_cl(2,j),1)
      end do

c aqueous phase-------------------------
c atmosphere
      bra1=0.
      bra2=0.
      bra3=0.
      bra4=0.
      cla1=0.
      cla2=0.
      cla3=0.
      cla4=0.
      xnaa2=0.
      bra1d=0.
      bra2d=0.
      bra3d=0.
      bra4d=0.
      cla1d=0.
      cla2d=0.
      cla3d=0.
      cla4d=0.
      xnaa2d=0.
      do k=2,n
         bra1=bra1+(sl1(42,1,k)+sl1(43,1,k)+sl1(44,1,k)+sl1(45,1,k)+
     &        2.*sl1(46,1,k)+sl1(47,1,k)+sl1(55,1,k)+sl1(j2-j3+9,1,k)
     &        +sion1(24,1,k)+2.*sion1(25,1,k)+sion1(26,1,k)+
     &    sion1(27,1,k)+sion1(28,1,k)+2.*sion1(29,1,k)+2.*sion1(38,1,k))
     &        *detw(k) 
         bra2=bra2+(sl1(42,2,k)+sl1(43,2,k)+sl1(44,2,k)+sl1(45,2,k)+
     &        2.*sl1(46,2,k)+sl1(47,2,k)+sl1(55,2,k)+sl1(j2-j3+9,2,k)
     &        +sion1(24,2,k)+2.*sion1(25,2,k)+sion1(26,2,k)+
     &    sion1(27,2,k)+sion1(28,2,k)+2.*sion1(29,2,k)+2.*sion1(38,2,k))
     &        *detw(k) 
         if (cl_out) then
         bra3=bra3+(sl1(42,3,k)+sl1(43,3,k)+sl1(44,3,k)+sl1(45,3,k)+
     &        2.*sl1(46,3,k)+sl1(47,3,k)+sl1(55,3,k)+sl1(j2-j3+9,3,k)
     &        +sion1(24,3,k)+2.*sion1(25,3,k)+sion1(26,3,k)+
     &    sion1(27,3,k)+sion1(28,3,k)+2.*sion1(29,3,k)+2.*sion1(38,3,k))
     &        *detw(k) 
         bra4=bra4+(sl1(42,4,k)+sl1(43,4,k)+sl1(44,4,k)+sl1(45,4,k)+
     &        2.*sl1(46,4,k)+sl1(47,4,k)+sl1(55,4,k)+sl1(j2-j3+9,4,k)
     &        +sion1(24,4,k)+2.*sion1(25,4,k)+sion1(26,4,k)+
     &    sion1(27,4,k)+sion1(28,4,k)+2.*sion1(29,4,k)+2.*sion1(38,4,k))
     &        *detw(k) 
         endif

         cla1=cla1+(sl1(30,1,k)+sl1(38,1,k)+sl1(39,1,k)+sl1(40,1,k)+
     &       2.*sl1(41,1,k)+sl1(47,1,k)+sl1(54,1,k)+sl1(58,1,k)
     &       +sl1(j2-j3+8,1,k)
     &       +sion1(14,1,k)+2.*sion1(15,1,k)+sion1(22,1,k)+sion1(23,1,k)
     &       +2.*sion1(28,1,k)+sion1(29,1,k)+2.*sion1(37,1,k))*detw(k)
         cla2=cla2+(sl1(30,2,k)+sl1(38,2,k)+sl1(39,2,k)+sl1(40,2,k)+
     &       2.*sl1(41,2,k)+sl1(47,2,k)+sl1(54,2,k)+sl1(58,2,k)
     &       +sl1(j2-j3+8,2,k)
     &       +sion1(14,2,k)+2.*sion1(15,2,k)+sion1(22,2,k)+sion1(23,2,k)
     &       +2.*sion1(28,2,k)+sion1(29,2,k)+2.*sion1(37,2,k))*detw(k)
        if (cl_out) then
         cla3=cla3+(sl1(30,3,k)+sl1(38,3,k)+sl1(39,3,k)+sl1(40,3,k)+
     &       2.*sl1(41,3,k)+sl1(47,3,k)+sl1(54,3,k)+sl1(58,3,k)
     &       +sl1(j2-j3+8,3,k)
     &       +sion1(14,3,k)+2.*sion1(15,3,k)+sion1(22,3,k)+sion1(23,3,k)
     &       +2.*sion1(28,3,k)+sion1(29,3,k)+2.*sion1(37,3,k))*detw(k)
         cla4=cla4+(sl1(30,4,k)+sl1(38,4,k)+sl1(39,4,k)+sl1(40,4,k)+
     &       2.*sl1(41,4,k)+sl1(47,4,k)+sl1(54,4,k)+sl1(58,4,k)
     &       +sl1(j2-j3+8,4,k)
     &       +sion1(14,4,k)+2.*sion1(15,4,k)+sion1(22,4,k)+sion1(23,4,k)
     &       +2.*sion1(28,4,k)+sion1(29,4,k)+2.*sion1(37,4,k))*detw(k)
        endif

         xnaa2=xnaa2+sion1(20,2,k)*detw(k)
      enddo

c + deposited
      bra1d=sl1(42,1,1)+sl1(43,1,1)+sl1(44,1,1)+sl1(45,1,1)+
     &     2.*sl1(46,1,1)+sl1(47,1,1)+sl1(55,1,1)+sl1(j2-j3+9,1,1)
     &     +sion1(24,1,1)+2.*sion1(25,1,1)+sion1(26,1,1)+
     &     sion1(27,1,1)+sion1(28,1,1)+2.*sion1(29,1,1)+
     &     2.*sion1(38,1,1)
      bra2d=sl1(42,2,1)+sl1(43,2,1)+sl1(44,2,1)+sl1(45,2,1)+
     &     2.*sl1(46,2,1)+sl1(47,2,1)+sl1(55,2,1)+sl1(j2-j3+9,2,1)
     &     +sion1(24,2,1)+2.*sion1(25,2,1)+sion1(26,2,1)+
     &     sion1(27,2,1)+sion1(28,2,1)+2.*sion1(29,2,1)+
     &     2.*sion1(38,2,1)
       if (cl_out) then
      bra3d=sl1(42,3,1)+sl1(43,3,1)+sl1(44,3,1)+sl1(45,3,1)+
     &     2.*sl1(46,3,1)+sl1(47,3,1)+sl1(55,3,1)+sl1(j2-j3+9,3,1)
     &     +sion1(24,3,1)+2.*sion1(25,3,1)+sion1(26,3,1)+
     &     sion1(27,3,1)+sion1(28,3,1)+2.*sion1(29,3,1)+
     &     2.*sion1(38,3,1)
      bra4d=sl1(42,4,1)+sl1(43,4,1)+sl1(44,4,1)+sl1(45,4,1)+
     &     2.*sl1(46,4,1)+sl1(47,4,1)+sl1(55,4,1)+sl1(j2-j3+9,4,1)
     &     +sion1(24,4,1)+2.*sion1(25,4,1)+sion1(26,4,1)+
     &     sion1(27,4,1)+sion1(28,4,1)+2.*sion1(29,4,1)+
     &     2.*sion1(38,4,1)
        endif

      cla1d=sl1(30,1,1)+sl1(38,1,1)+sl1(39,1,1)+sl1(40,1,1)+
     &     2.*sl1(41,1,1)+sl1(47,1,1)+sl1(54,1,1)+sl1(58,1,1)
     &     +sl1(j2-j3+8,1,1)
     &     +sion1(14,1,1)+2.*sion1(15,1,1)+sion1(22,1,1)+sion1(23,1,1)
     &     +2.*sion1(28,1,1)+sion1(29,1,1)+2.*sion1(37,1,1)
      cla2d=sl1(30,2,1)+sl1(38,2,1)+sl1(39,2,1)+sl1(40,2,1)+
     &     2.*sl1(41,2,1)+sl1(47,2,1)+sl1(54,2,1)+sl1(58,2,1)
     &     +sl1(j2-j3+8,2,1)
     &     +sion1(14,2,1)+2.*sion1(15,2,1)+sion1(22,2,1)+sion1(23,2,1)
     &     +2.*sion1(28,2,1)+sion1(29,2,1)+2.*sion1(37,2,1)
       if (cl_out) then
      cla3d=sl1(30,3,1)+sl1(38,3,1)+sl1(39,3,1)+sl1(40,3,1)+
     &     2.*sl1(41,3,1)+sl1(47,3,1)+sl1(54,3,1)+sl1(58,3,1)
     &     +sl1(j2-j3+8,3,1)
     &     +sion1(14,3,1)+2.*sion1(15,3,1)+sion1(22,3,1)+sion1(23,3,1)
     &     +2.*sion1(28,3,1)+sion1(29,3,1)+2.*sion1(37,3,1)
      cla4d=sl1(30,4,1)+sl1(38,4,1)+sl1(39,4,1)+sl1(40,4,1)+
     &     2.*sl1(41,4,1)+sl1(47,4,1)+sl1(54,4,1)+sl1(58,4,1)
     &     +sl1(j2-j3+8,4,1)
     &     +sion1(14,4,1)+2.*sion1(15,4,1)+sion1(22,4,1)+sion1(23,4,1)
     &     +2.*sion1(28,4,1)+sion1(29,4,1)+2.*sion1(37,4,1)
      endif

      xnaa2d=sion1(20,1,1)+sion1(20,2,1)+sion1(20,3,1)+sion1(20,4,1)

c output


 100  continue
      open (74,file='mass.out',status='unknown', position='append'
     &     ,err=100)
      write (74,10) lday,lst,lmin
       if (cl_out) then
        write (74,20) brg, brgd, bra1, bra2, bra1d,bra2d, brsss*detw(1), 
     &     brg+bra1+bra2+brgd+bra1d+bra2d-brsss*detw(1)
        write (74,20) brg, brgd, bra3, bra4,bra3d, bra4d, brsss*detw(1), 
     &     brg+bra3+bra4+brgd+bra3d+bra4d-brsss*detw(1),brg+bra1+bra2+
     &     brgd+bra1d+bra2d+bra3+bra4+brgd+bra3d+bra4d-brsss*detw(1)
        write (74,25) clg, clgd, cla1, cla2,cla1d, cla2d, clsss*detw(1), 
     &     clg+cla1+cla2+clgd+cla1d+cla2d-clsss*detw(1)
        write (74,25) clg, clgd, cla3, cla4,cla3d, cla4d, clsss*detw(1), 
     &     clg+cla3+cla4+clgd+cla3d+cla4d-clsss*detw(1),clg+cla1+cla2+
     &     clgd+cla1d+cla2d+cla3+cla4+clgd+cla3d+cla4d-clsss*detw(1)
        write (74,30) xnaa2, xnaa2d,xnasss*detw(1),xnaa2+xnaa2d-xnasss*
     &     detw(1)
       else
        write (74,20) brg, brgd, bra1, bra2, bra1d,bra2d, brsss*detw(1), 
     &     brg+bra1+bra2+brgd+bra1d+bra2d-brsss*detw(1)
        write (74,25) clg, clgd, cla1, cla2, cla1d,cla2d, clsss*detw(1), 
     &     clg+cla1+cla2+clgd+cla1d+cla2d-clsss*detw(1)
        write (74,30) xnaa2, xnaa2d,  xnasss*detw(1), xnaa2+xnaa2d
     &     -xnasss*detw(1)
       endif
      close (74)

 10   format(3i3)
 20   format('br: ',8d15.7)
 25   format('cl: ',8d15.7)
 30   format('na: ',45x,d15.7,15x,3d15.7)

      endif

      end subroutine mass_ch



c
c----------------------------------------------------------------
c


      subroutine ave_parms (n_bl,nz_box) 
c average cw and rc over BL in box model if BL_box=.true.
c      implicit double precision (a-h,o-z)

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'gas_Parameters.h' !additional common blocks and other definitions 

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_dryg/ xkmtd(n,2,NSPEC),henry(n,NSPEC),xeq(n,NSPEC)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb unused now

c note: activity coeff. are being correctly calculated with the averaged 
c concentrations in layer n_bl

c start averaging one level above "working level" 
c - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3

      do k=nstart,nz_box
         cwsum  = 0.
         rcsum = 0.
         cmsum  = 0.
!         cwmsum = 0.
!         acmsum = 0. ! jjb removed
!         convsum= 0. ! jjb unused below, now removed
         do kc=1,nkc
            cwsum  = cwsum  + cw(kc,k)
            rcsum = rcsum + rc(kc,k)
            cmsum  = cmsum  + cm(kc,k)
!            acmsum = acmsum + acm(k,kc) ! jjb removed
!            convsum= convsum+ conv2(k,kc) ! jjb unused
         enddo
         cw(kc,n_bl)  = cwsum / (nz_box-nstart+1)
         rc(kc,n_bl) = rcsum / (nz_box-nstart+1)
         cm(kc,n_bl)  = cmsum / (nz_box-nstart+1)
!         cwm(kc,n_bl) = cwmsum / (nz_box-nstart+1)
!         acm(n_bl,kc) = acmsum / (nz_box-nstart+1) ! jjb removed
      enddo


      am3s = 0.
!      am32s = 0.
      cm3s = 0.
!      cm32s = 0.
      ps    = 0.
      rhos  = 0.
      ts    = 0.
      
      do k=nstart,nz_box
         am3s = am3s + am3(k)
!         am32s = am32s + am3(k,2)
         cm3s = cm3s + cm3(k)
!         cm32s = cm32s + cm3(k,2)
         ps    = ps    + p(k)
         rhos  = rhos  + rho(k)
         ts    = ts    + t(k)
      enddo
      am3(n_bl) = am3s / (nz_box-nstart+1)
!      am3(n_bl,2) = am32s / (nz_box-nstart+1)
      cm3(n_bl) = cm3s / (nz_box-nstart+1)
!      cm3(n_bl,2) = cm32s / (nz_box-nstart+1)

      p(n_bl)     = ps   / (nz_box-nstart+1)
      rho(n_bl)   = rhos / (nz_box-nstart+1)
      t(n_bl)     = ts   / (nz_box-nstart+1)


      do j=1,NSPEC
         do kc=1,2
            xkmtds = 0.
            do k=nstart,nz_box
               xkmtds = xkmtds + xkmtd(k,kc,j)
            enddo
            xkmtd(n_bl,kc,j) = xkmtds / (nz_box-nstart+1)
         enddo
      enddo

      end subroutine ave_parms


c
c-------------------------------------------------------
c

      subroutine ave_j (nz_box,n_bl)
c calculate average photolysis rates for BL and put these values on level n_bl

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nphrxn

      implicit double precision (a-h,o-z)

      common /band_rat/ photol_j(nphrxn,n)
      double precision photol_j

c start averaging one level above "working level" 
c - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3

      do j=1,nphrxn ! jjb
         xjsum = 0.
         do k=nstart,nz_box
            xjsum = xjsum + photol_j(j,k)
         enddo
         photol_j(j,n_bl) = xjsum/(nz_box-nstart+1)   
      enddo

      end subroutine ave_j

c
c----------------------------------------------------------------
c

      subroutine ave_aer (n_bl,nz_box) 
c average k_mt over BL in box model if BL_box=.true.

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions          

      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
!     dimension fs(n,nka),ffsum(nka) ! jjb both not used here

c start averaging one level above "working level" 
c - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3

 
      do j=1,NSPEC
         do kc=1,2
            xkmsum = 0.
            xkfsum = 0.
            xkbsum = 0.
            do k=nstart,nz_box
               xkmsum = xkmsum + xkmt(k,kc,j)
               xkfsum = xkfsum + xkef(k,kc,j)
               xkbsum = xkbsum + xkeb(k,kc,j)
            enddo
            xkmt(n_bl,kc,j) = xkmsum / (nz_box-nstart+1)
            xkef(n_bl,kc,j) = xkfsum / (nz_box-nstart+1)
            xkeb(n_bl,kc,j) = xkbsum / (nz_box-nstart+1)
         enddo
      enddo

      do j=1,NSPEC
         xhensum = 0.
         xalsum  = 0.
         xvmsum  = 0.
         do k=nstart,nz_box
            xhensum = xhensum + henry(j,k)
            xalsum  = xalsum  + alpha(j,k)
            xvmsum  = xvmsum  + vmean(j,k)
         enddo
         henry(j,n_bl)  = xhensum / (nz_box-nstart+1)
         alpha(j,n_bl)  = xalsum  / (nz_box-nstart+1)
         vmean(j,n_bl)  = xvmsum  / (nz_box-nstart+1)
      enddo

c particle size distribution
c some logical mistake - leeds to very low LWC after next call to SR cw_rc
c      do ia=1,nka
c         ffsum(ia) = 0.
c      enddo
cc "dry" the aerosol and sum up
c      do k=nstart,nz_box
c         do ia=1,nka
c            do jt=1,nkt
c               fs(k,ia)=fs(k,ia)+ff(jt,ia,k)
c            enddo
c            ffsum(ia)=ffsum(ia)+fs(k,ia)
c         enddo
c      enddo
cc average
c      do ia=1,nka
c         ff(1,ia,n_bl)=ffsum(ia)  / (nz_box-nstart+1)
c         do jt=2,nkt
c            ff(jt,ia,n_bl)=0.
c         enddo
c      enddo

      end subroutine ave_aer

c
c----------------------------------------------------------------
c

      subroutine ave_tot (n_bl,nz_box) 
c average k_mt over BL in box model if BL_box=.true.

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h' !additional common blocks and other definitions          

      common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)

c start averaging one level above "working level" 
c - semi hard coded for n_bl = 2
      nstart = 2
      if (n_bl.eq.2) nstart = 3


      do j=1,NSPEC
         do kc=1,nkc
            xkmsum = 0.
            xkfsum = 0.
            xkbsum = 0.
            do k=nstart,nz_box
               xkmsum = xkmsum + xkmt(k,kc,j)
               xkfsum = xkfsum + xkef(k,kc,j)
               xkbsum = xkbsum + xkeb(k,kc,j)
            enddo
            xkmt(n_bl,kc,j) = xkmsum / (nz_box-nstart+1)
            xkef(n_bl,kc,j) = xkfsum / (nz_box-nstart+1)
            xkeb(n_bl,kc,j) = xkbsum / (nz_box-nstart+1)
         enddo
      enddo

      do j=1,NSPEC
         xhensum = 0.
         xalsum  = 0.
         xvmsum  = 0.
         do k=nstart,nz_box
            xhensum = xhensum + henry(j,k)
            xalsum  = xalsum  + alpha(j,k)
            xvmsum  = xvmsum  + vmean(j,k)
         enddo
         henry(j,n_bl)  = xhensum / (nz_box-nstart+1)
         alpha(j,n_bl)  = xalsum  / (nz_box-nstart+1)
         vmean(j,n_bl)  = xvmsum  / (nz_box-nstart+1)
      enddo

      end subroutine ave_tot

c
c----------------------------------------------------------------
c

      subroutine print_k_mt_a


      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     nkc

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions          

      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)

      print *,'print aerosol kmt'
      do k=1,nf
         print *,k,xkmt(k,1,ind_O3),xkmt(k,2,ind_O3)
      enddo

      end subroutine print_k_mt_a



c
c-------------------------------------------------------
c

      subroutine set_box_gas (nlevbox,n_bl) 
c     pick the values from the designated level: nlevbox

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit double precision (a-h,o-z)

      include 'gas_Parameters.h' !additional common blocks and other definitions

      common /kpp_dryg/ xkmtd(n,2,NSPEC),henry(n,NSPEC),xeq(n,NSPEC)

      do kc=1,2!nkc
         do j=1,NSPEC
            xkmtd(n_bl,kc,j) = xkmtd(nlevbox,kc,j)
         enddo
      enddo

      end subroutine set_box_gas

c
c-------------------------------------------------------
c

      subroutine set_box_lev_a (nlevbox,n_bl) 
c     pick the values from the designated level: nlevbox

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nrlay,
     &     nka,
     &     nkt,
     &     nkc,
     &     mb

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions

      common /cb11/ totrad (mb,nrlay)
      double precision totrad

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck13/ conv2(nkc,n)
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)
      common /kpp_drya/ xkmtd(nf,2,NSPEC),xeq(nf,NSPEC)

      do kc=1,2!nkc
!         acm(n_bl,kc)  = acm(nlevbox,kc) ! jjb removed
!        conv2(n_bl,kc)= conv2(nlevbox,kc)! jjb updated
         conv2(kc,n_bl) = conv2(kc,nlevbox)
         cw(kc,n_bl)   = cw(kc,nlevbox)
         rc(kc,n_bl)  = rc(kc,nlevbox)
         do j=1,NSPEC
            xkmt(n_bl,kc,j) = xkmt(nlevbox,kc,j)
            xkef(n_bl,kc,j) = xkef(nlevbox,kc,j)
            xkeb(n_bl,kc,j) = xkeb(nlevbox,kc,j)
            xkmtd(n_bl,kc,j) = xkmtd(nlevbox,kc,j)
         enddo
      enddo

      am3(n_bl) = am3(nlevbox)
!      am3(n_bl,2) = am3(nlevbox,2)
      cm3(n_bl) = cm3(nlevbox)
!      cm3(n_bl,2) = cm3(nlevbox,2)

      p(n_bl)     = p(nlevbox)
      rho(n_bl)   = rho(nlevbox)
      t(n_bl)     = t(nlevbox)
      totrad(1,n_bl) = totrad(1,nlevbox) ! jjb check this. Why is only the first wavelength band picked up?

      do j=1,NSPEC
         alpha(j,n_bl) = alpha(j,nlevbox)
         henry(j,n_bl) = henry(j,nlevbox)
         vmean(j,n_bl) = vmean(j,nlevbox)
      enddo

c SR cw_rc is still being called, therefore also init f
      do ia=1,nka
         do jt=1,nkt
            ff(jt,ia,n_bl) = ff(jt,ia,nlevbox) 
         enddo
      enddo

      end subroutine set_box_lev_a

c
c-------------------------------------------------------
c


      subroutine set_box_lev_t (nlevbox,n_bl) 
c     pick the values from the designated level: nlevbox    

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nrlay,
     &     nkc,
     &     mb

      implicit double precision (a-h,o-z)

      include 'tot_Parameters.h' !additional common blocks and other definitions          

      common /cb11/ totrad (mb,nrlay)
      double precision totrad

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck13/ conv2(nkc,n)
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_ltot/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2tot/ alpha(NSPEC,nf),vmean(NSPEC,nf)
!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated

      xph3=0.
      xph4=0.
      if (cm(n_bl,3).gt.0.) xph3 = 1.
      if (cm(n_bl,4).gt.0.) xph4 = 1.

      if (xph3.eq.1..and.xph4.eq.1.) return

      do kc=1,nkc
!         acm(n_bl,kc)  = acm(nlevbox,kc) ! jjb removed
!        conv2(n_bl,kc)= conv2(nlevbox,kc) ! jjb updated
         conv2(kc,n_bl) = conv2(kc,nlevbox)
         cw(kc,n_bl) = cw(kc,nlevbox)
         rc(kc,n_bl) = rc(kc,nlevbox)
         do j=1,NSPEC
            xkmt(n_bl,kc,j) = xkmt(nlevbox,kc,j)
            xkef(n_bl,kc,j) = xkef(nlevbox,kc,j)
            xkeb(n_bl,kc,j) = xkeb(nlevbox,kc,j)
         enddo
      enddo

      am3(n_bl) = am3(nlevbox)
!      am3(n_bl,2) = am3(nlevbox,2)
      cm3(n_bl) = cm3(nlevbox)
!      cm3(n_bl,2) = cm3(nlevbox,2)

      p(n_bl)     = p(nlevbox)
      rho(n_bl)   = rho(nlevbox)
      t(n_bl)     = t(nlevbox)
      totrad(1,n_bl) = totrad(1,nlevbox) ! jjb check this. Why is only the first wavelength band picked up?

      do j=1,NSPEC
         alpha(j,n_bl) = alpha(j,nlevbox)
         henry(j,n_bl) = henry(j,nlevbox)
         vmean(j,n_bl) = vmean(j,nlevbox)
      enddo

      end subroutine set_box_lev_t


c
c-------------------------------------------------------
c

      subroutine print_vals (nlevbox,n_bl)
c     test output

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nrlay,
     &     nkc,
     &     mb

      implicit double precision (a-h,o-z)

      include 'aer_Parameters.h' !additional common blocks and other definitions          

      common /cb11/ totrad (mb,nrlay)
      double precision totrad

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck13/ conv2(nkc,n)
!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      common /kpp_2aer/ alpha(NSPEC,nf),vmean(NSPEC,nf)


      print *,'output info', lst,lmin

      k=n_bl

      print *,'k = ',k
      do kc=1,2!nkc
         print *,'kc = ',kc
!        print *,acm(k,kc),conv2(k,kc),cw(k,kc) ! jjb acm removed
         print *,conv2(kc,k),cw(kc,k)
         print *,rc(kc,k),cm(kc,k)
         do j=1,NSPEC
            print *,xkmt(k,kc,j),xkef(k,kc,j),xkeb(k,kc,j)
         enddo
      enddo

      print *,'parameters'
!     print *,am3(k,1),am3(k,2),p(k) ! jjb am3 CO removed
      print *,am3(k),p(k)
!     print *,cm3(k,1),cm3(k,2),rho(k) ! jjb cm3 CO removed
      print *,cm3(k),rho(k)
      print *,t(k), totrad(1,k)

      print *,'alpha,henry,vmean'
      do j=1,NSPEC
         print *,alpha(j,k), henry(j,k), vmean(j,k)
      enddo

      k=nlevbox

      print *,'crap'

      print *,'k = ',k
      do kc=1,2!nkc
         print *,'kc = ',kc
!        print *,acm(k,kc),conv2(k,kc),cw(k,kc) ! jjb acm removed
         print *,conv2(kc,k),cw(kc,k)
         print *,rc(kc,k),cm(kc,k)
         do j=1,NSPEC
            print *,xkmt(k,kc,j),xkef(k,kc,j),xkeb(k,kc,j)
         enddo
      enddo

      print *,'parameters'
!     print *,am3(k,1),am3(k,2),p(k) ! jjb am3 CO removed
      print *,am3(k),p(k)    
!     print *,cm3(k,1),cm3(k,2),rho(k) ! jjb cm3 CO removed
      print *,cm3(k),rho(k)  
      print *,t(k), totrad(1,k)

      print *,'alpha,henry,vmean'
      do j=1,NSPEC
         print *,alpha(j,k), henry(j,k), vmean(j,k)
      enddo

      print *,' used for k=n_bl'
      k=n_bl
      do kc=1,2                 !nkc
         print *,k,kc,cw(kc,k)
         print *,xkmt(k,kc,ind_O3),xkmt(k,kc,ind_HCl),
     &           xkmt(k,kc,ind_HNO3)
      enddo

      print *,'now column vals'
      do k=1,nf
         do kc=1,2!nkc
            print *,k,kc,cw(kc,k)
            print *,xkmt(k,kc,ind_O3),xkmt(k,kc,ind_HCl),
     &           xkmt(k,kc,ind_HNO3)
         enddo
      enddo

      end subroutine print_vals

c
c-------------------------------------------------------
c
! jjb commented as unused 24/03/2016

c$$$      subroutine gamma_surf (box,n_bl)
c$$$c     calculation of reaction rate coefficients for surface reactions
c$$$c     note that the value that is used here for the accommodation
c$$$c     coefficient is in fact a reaction probability (=gamma)
c$$$c     after Knipping and Dabdub 2002 (#1692) and von Glasow (2006)
c$$$
c$$$c transfer coefficient after Schwarz, 1986 (see Sander & Crutzen '96, JGR, 9127)
c$$$c but no mean values used (like in SR k_mt_a/t) but integrated values
c$$$
c$$$      USE global_params, ONLY :
c$$$! Imported Parameters:
c$$$     &     j2,
c$$$     &     j6,
c$$$     &     nf,
c$$$     &     n,
c$$$     &     nka,
c$$$     &     nkt,
c$$$     &     nkc
c$$$
c$$$      implicit double precision (a-h,o-z)
c$$$      logical box
c$$$
c$$$      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
c$$$     &              e(nkt),dew(nkt),rq(nkt,nka)
c$$$      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
c$$$      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
c$$$      double precision theta, thetl, t, talt, p, rho
c$$$      common /blck06/ kw(nka),ka
c$$$      common /blck12/ cw(nkc,n),cm(nkc,n)
c$$$      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
c$$$!     common /kpp_1/ am3(n,2), cm3(n,2),cw(nf,nkc),conv2(nf,nkc),xconv1 ! jjb old CB, updated
c$$$!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
c$$$!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
c$$$!     common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc) ! jjb none of the objects of the common block is used
c$$$!     common /kpp_mol/ cm(nf,nkc),xgamma(nf,j6,nkc) ! jjb updated
c$$$      common /k_surf/ xkmt_OHClm(nf,nkc)
c$$$
c$$$      dimension freep(nf),rqm(nkt,nka),xkmt_surf(nf,nkc)
c$$$
c$$$      func(a,k)=dsqrt(t(k)/a)*4.60138
c$$$
c$$$c change rq in um to rqm in m
c$$$      do ia=1,nka
c$$$         do jt=1,nkt
c$$$            rqm(jt,ia)=rq(jt,ia)*1.d-6
c$$$         enddo
c$$$      enddo
c$$$
c$$$      nmin=2
c$$$      nmax=nf
c$$$      if (box) then
c$$$        nmin=n_bl
c$$$        nmax=n_bl
c$$$      endif
c$$$
c$$$c free path length (lambda=freep):
c$$$      do k=nmin,nmax
c$$$         freep(k)=2.28d-5 * t(k) / p(k)
c$$$      enddo
c$$$c loop over vertical grid
c$$$      do 1000 k=nmin,nmax
c$$$         vmean_OH= func(1.7d-2,k)
c$$$c loop over the nkc different chemical bins
c$$$         do kc=1,4
c$$$            if (cm(kc,k).eq.0.) goto 1001 ! switch changed from cw
c$$$c define summation limits (1) ---
c$$$            if (kc.eq.1) then
c$$$               iia_0=1
c$$$               iia_e=ka
c$$$            endif
c$$$            if (kc.eq.2) then
c$$$               iia_0=ka+1
c$$$               iia_e=nka
c$$$            endif
c$$$            if (kc.eq.3) then
c$$$               iia_0=1
c$$$               iia_e=ka
c$$$            endif
c$$$            if (kc.eq.4) then
c$$$               iia_0=ka+1
c$$$               iia_e=nka
c$$$            endif
c$$$c fast version without logarithmic integration
c$$$            x1=0.
c$$$            xk1=0.
c$$$c --- OH + Cl- --> 0.5 Cl2 + OH- --- 
c$$$c "best guess": use alpha from Knipping and Dabdub and gas phase limitation:
c$$$c alpha = 0.02 * gamma_s * [Cl-], where [Cl-] is in mol/l, gamma_s=2
c$$$c C(ind_Clmlx) / LWC *1.d-3:  in mol/m3_air * m3_air/m3_aq * m3_aq/l_aq
c$$$            if (cw(kc,k).gt.0.) then
c$$$               gamma = min(1.,4.d-5*sion1(14,kc,k)/cw(kc,k))
c$$$            else
c$$$               gamma = 0.d0
c$$$            end if
c$$$            if (gamma.gt.0.) x1= 4./(3. * gamma)
c$$$c -------------------------------------------
c$$$            do ia=iia_0,iia_e
c$$$c define summation limits (2)
c$$$               if (kc.eq.1) then
c$$$                  jjt_0=1
c$$$                  jjt_e=kw(ia)
c$$$               endif
c$$$               if (kc.eq.2) then
c$$$                  jjt_0=1
c$$$                  jjt_e=kw(ia)
c$$$               endif
c$$$               if (kc.eq.3) then
c$$$                  jjt_0=kw(ia)+1
c$$$                  jjt_e=nkt
c$$$               endif
c$$$               if (kc.eq.4) then
c$$$                  jjt_0=kw(ia)+1
c$$$                  jjt_e=nkt
c$$$               endif
c$$$               do jt=jjt_0,jjt_e
c$$$c conversion: um      --> m               :10^-6
c$$$                  rqq=rqm(jt,ia)
c$$$c kt=1./(r^2/(3*D_g)+4*r/(3*vmean*alpha))=vmean/r*1/(r/lambda+4/(3*alpha))
c$$$c     with D_g=lambda*vmean/3.
c$$$c here a volume weighted value is calculated, therefore weighting with r^3:
c$$$c kmt=4/3*piL*/sum(a)*sum(r){r^3*N*kt}
c$$$c --- OH + Cl- --> 0.5 Cl2 + OH- ----
c$$$                  x2=vmean_OH/(rqq/freep(k)+x1) ![1/s]
c$$$c conversion: 1/cm^3 --> 1/m^3(air):10^6
c$$$                  xk1=xk1+x2*rqq*rqq*ff(jt,ia,k)*1.d6
c$$$               enddo            !jt
c$$$            enddo               !ia
c$$$c k_mt=4*pi/(3*LWC)*sum
c$$$            if (cw(kc,k).gt.0.d0) then
c$$$!               xkmt_surf(k,kc)=4.*3.1415927/(3.*cw(kc,k))*xk1 ![1/s]
c$$$               xkmt_surf(k,kc)=4.*pi/(3.*cw(kc,k))*xk1 ![1/s]
c$$$            end if
c$$$
c$$$c kmt is supposed to be a first-order rate coefficient but in KPP kmt will 
c$$$c be multiplied by [Cl-]/[Br-], so divide by [X-] here:
c$$$
c$$$            if (sion1(14,kc,k).gt.0.) 
c$$$     &           xkmt_OHClm(k,kc)=xkmt_surf(k,kc)/sion1(14,kc,k)  
c$$$
c$$$c            print *,'SR gamma_surf'
c$$$c            print *,lday,lst,lmin
c$$$c            print *,k,kc,cw(k,kc)!,sion1(14,kc,k)
c$$$c            print *,xkmt_OHClm(k,kc)
c$$$c multiply with conc of X- to get "real" kmt:
c$$$c            print *,xkmt_OHClm(k,kc)*sion1(14,kc,k),xkmt_OHClm(k,kc)*
c$$$c     &           sion1(24,kc,k)
c$$$
c$$$ 1001       continue
c$$$         enddo                  ! kc
c$$$         
c$$$ 1000 continue                  !k
c$$$
c$$$      end subroutine gamma_surf



c
c-------------------------------------------------------------
c

c functions for KPP

c third body reaction rate constants are formulated as pseudo 
c second order: cm^3/s, so [aircc]=cm^3/s has to be used

      double precision function farr (a,b)
c Arrhenius relation
      implicit none

      double precision a
      integer b

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      farr=a*exp(b/te)
      end function farr

c----------------------------------------------------------------

      double precision function farr_sp (a,b,c,d)
c Arrhenius relation (with complex temperature dependence)
      implicit none

      double precision a,c
      integer b,d

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      farr_sp=a*((te/b)**c)*exp(d/te)
      end function farr_sp

c----------------------------------------------------------------

      double precision function ATK_3 (a1,a2,b1,b2,fc)
c calculate third body reactions according to Atkinson '92
      implicit none

      double precision a0,b0,a1,a2,b1,b2,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2
      x2=fc

      atk_3=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      end function ATK_3

c----------------------------------------------------------------

      double precision function ATK_3a (a1,a2,b1,b2,tfc)
c calculate third body reactions according to Atkinson '92
      double precision a0,b0,a1,a2,b1,b2,tfc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2
      x2=exp(-te/tfc)
      atk_3a=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      end function ATK_3a

c----------------------------------------------------------------

      double precision function ATK_3c (a1,b1,fc)
c calculate third body reactions according to Atkinson '92
      implicit none

      double precision a0,b0,a1,b1,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*exp(-10000./te)*aircc
      b0=b1*exp(-10900./te)
      x2=fc
      if (fc.eq.0.) x2=exp(-te/250.)+exp(-1050./te)
      atk_3c=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      end function ATK_3c

c----------------------------------------------------------------

      double precision function ATK_3d (a1,b1,fc)
c calculate third body reactions according to IUPAC 2004
      implicit none

      double precision a0,b0,a1,b1,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*exp(-8000./te)*aircc
      b0=b1*exp(-8820./te)
      x2=fc
      atk_3d=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      end function ATK_3d

c----------------------------------------------------------------

      double precision function ATK_3e (a1,a2,b1,b2,fc)
c calculate third body reactions according to Atkinson '92 (used for OH + OIO)
c  NOT USED: Plane et al., 2006 reported a pressure dependent rate coefficient
c  from a theoretical study, however they suggested an Arrhenius relationship
c  for use in atmospheric modelling
      implicit none

      double precision a0,b0,a1,a2,b1,b2,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2*exp(46./te)
      x2=fc

      atk_3e=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      end function ATK_3e

c----------------------------------------------------------------

      double precision function ATK_3f (a1,a2,b1,b2,fc)
c calculate third body reactions according to Atkinson '92 (used for OCLO + O3P)
      implicit none

      double precision a0,b0,a1,a2,b1,b2,fc,x2

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*aircc*(te/298.)**a2
      b0=b1*(te/298.)**b2
      x2=fc

      atk_3f=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      end function ATK_3f

c----------------------------------------------------------------

      double precision function sHNO3 (a1,b1,a2,b2,a3,b3)
c calculate special rate function for OH + HNO3 (JPL 2003)
      implicit none

      double precision a0,a1,a2,a3,func,tte
      integer b0,b1,b2,b3

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      func(a0,b0)=a0*exp(b0*tte)
      tte=1./te
      sHNO3=func(a1,b1)+(func(a3,b3)*aircc/
     &     (1+func(a3,b3)*aircc/func(a2,b2)))
      end function sHNO3

c----------------------------------------------------------------

      double precision function fbck (a1,a2,b1,b2,fc,ak,bk)
c calculate thermal decomposition rate from forward and
c equilibrium rate
      implicit none

      double precision a0,b0,a1,a2,b1,b2,fc,x1,x2,ak,bk

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.d0)**a2
      b0=b1*(te/300.d0)**b2
      x2=fc

      x1=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      fbck=x1/(ak*dexp(bk/te))
      end function fbck

c----------------------------------------------------------------

      double precision function fbckJ (a1,a2,b1,b2,ak,bk)
c calculate thermal decomposition rate from forward and
c equilibrium rate
      implicit none

      double precision a0,b0,a1,a2,b1,b2,x1,x2,ak,bk

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      a0=a1*aircc*(te/300.d0)**a2
      b0=b1*(te/300.d0)**b2
      x2= 0.6d0

      x1=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      fbckJ=x1/(ak*dexp(bk/te))
      end function fbckJ

c----------------------------------------------------------------

      double precision function fbck2 (a1,a2,b1,b2,fc,ck)
c calculate thermal decomposition rate from forward and
c equilibrium rate (used for BrNO3 decomposition)
      implicit none

      double precision a0,b0,a1,a2,b1,b2,fc,x1,x2,ak,bk,ck

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

c parameters to calculate K_eq in atm-1 (Orlando and Tyndall, 1996)
      ak=5.44d-9
      bk=14192.d0
      a0=a1*aircc*(te/300.)**a2
      b0=b1*(te/300.)**b2
      x2=fc

      x1=(a0/(1+a0/b0))*(x2**(1/(1+dlog10(a0/b0)*
     &     dlog10(a0/b0))))
      fbck2 = 0.d0
      if (ck.ne.0.d0) fbck2=x1/(ak*dexp(bk/te)*8.314/101325.*te/ck)
      end function fbck2


c----------------------------------------------------------------

      double precision function sp_17_old (a1)
c calculate special rate function for rxn 17 (OH + CO)
      implicit none

      double precision a1
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      sp_17_old=a1*(1+0.6*pk/101325.) !(pressure in atm)
      end function sp_17_old

c----------------------------------------------------------------

      double precision function sp_17 (a,b)
c calculate special rate function for rxn 17 (OH + CO)
      implicit none

      double precision a,b
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      sp_17=a*(1.d0+aircc/b) ! simpler IUPAC parametrization
c$$$      sp_17=a1*(1+0.6*pk/101325.) !(pressure in atm)
      end function sp_17

c----------------------------------------------------------------

      double precision function sp_23 (a1,b1,a2,b2,a3,b3)
c calculate special rate function for rxn 23 (HO2+HO2 - including H2O correction)
      implicit none

      double precision a0,a1,a2,a3,func,tte
      integer b0,b1,b2,b3
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      func(a0,b0)=a0*exp(b0*tte)
      tte=1./te
      sp_23=(func(a1,b1)+func(a2*aircc,b2))*
     &     (1+(func(a3*aircc*h2oppm*1.0d-6,b3)))
      end function sp_23

c----------------------------------------------------------------

      double precision function sp_29 (a1,b1,a2,b2,c)
c calculate special rate function for rxn 29
      implicit none

      double precision a1,b1,a2,b2,c,fun1,fun2,fun3,num,den,z,
     & a0,b0,c0,d0
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      fun1(a0,b0,c0)=a0*(b0**c0)
      fun2(a0,b0)=1./(1.+(dlog10(a0/b0))*(dlog10(a0/b0)))
      fun3(a0,b0,c0,d0)=a0/(1.+a0/b0)*(c0**d0)

      num=aircc*fun1(a1,te,b1)
      den=    fun1(a2,te,b2)
      z=fun2(num,den)
      sp_29=fun3(num,den,c,z)    
      end function sp_29

c----------------------------------------------------------------

      double precision function fcn (x1)
c rate constant for thermal decomposition of ClNO3 (#2730)
      implicit none

      double precision x1,x2,xmg
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      x2=8.314*te
      xmg=pk/x2
      fcn=10**(-6.16)*dexp(-90.7d3/x2)*xmg*x1
      end function fcn

c----------------------------------------------------------------

      double precision function farr2 (a0,b0)
c Arrhenius function but with b0 as the value for T=298K
      implicit none

      double precision a0
      integer b0
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk
c 1/298.=3.3557d-3
      farr2=a0*exp(dble(b0)*(1.d0/te-3.3557d-3))
      end function farr2

c----------------------------------------------------------------

      double precision function fhet_t (a0,b0,c0)
c heterogeneous rate function 
c ClFCT     = 5.0D2                ; factor for H02/H01, i.e Cl-/H2O
c BrFCT     = 3.0D5                ; factor for H03/H01, i.e Br-/H2O
c a0=1..4  liquid size class
c b0=1..3  branch of het reaction: H2O, Cl-, Br-
c c0=1..3  gas phase reactant:     N2O5, ClNO3, BrNO3

      implicit none

      INCLUDE 'tot_Parameters.h'
      INCLUDE 'tot_Global.h'
      integer, intent(in) :: a0,b0,c0

      double precision :: h2oa, hetT
      double precision :: xbr, xtr

      if (a0.eq.1) then
         h2oa=FIX(indf_H2Ol1)
         hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)

      else if (a0.eq.2) then
         h2oa=FIX(indf_H2Ol2)
         hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)

      else if (a0.eq.3) then
         h2oa=FIX(indf_H2Ol3)
         hetT=h2oa + 5.0D2*C(ind_Clml3) + 3.0D5*C(ind_Brml3)

      else if (a0.eq.4) then
         h2oa=FIX(indf_H2Ol4)
         hetT=h2oa + 5.0D2*C(ind_Clml4) + 3.0D5*C(ind_Brml4)

      else ! undefined case, shouldn't happend
         stop 'Wrong a0 index in function fhet_t'
      endif


      if (b0.eq.1) then
         xbr=h2oa

      else if (b0.eq.2) then
         xbr=5.0D2

      else if (b0.eq.3) then
         xbr=3.0D5

      else ! undefined case, shouldn't happend
         stop 'Wrong b0 index in function fhet_t'
      endif


      if (c0.eq.1) then
         xtr=yxkmt(a0,ind_N2O5)
      else if (c0.eq.2) then
         xtr=yxkmt(a0,ind_ClNO3)
      else if (c0.eq.3) then
         xtr=yxkmt(a0,ind_BrNO3)
      else ! undefined case, shouldn't happend
         stop 'Wrong c0 index in function fhet_t'
      endif


      if (hetT.gt.0.d0) then
         fhet_t=xtr * ycw(a0) * xbr /hetT
      else
         fhet_t=0.
      endif

      end function fhet_t

c----------------------------------------------------------------

      double precision function fliq_60 (a1,b1,c,d)
c calculate special rate function for rxn 60

      implicit none

      double precision a1,c,d
      integer b1
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      if (d.gt.0.d0) then
c         fliq_60=farr2(a1,b1)*c/(c+0.1/d)
         fliq_60=a1*dexp(dble(b1)*(1.d0/te-3.3557d-3))*c/(c+0.1d0/d)
      else
         fliq_60=0.
      endif
      end

c----------------------------------------------------------------

      double precision function dmin2 (a)
c confine rate constant to upper limit (diffusion control)
c a=k; dclim=upper limit due to diffusion-control 

      implicit none

      double precision a,dclim

      dclim = 1.d10

      dmin2 = dmin1(a,dclim )

      end

c----------------------------------------------------------------

      double precision function dmin3 (a)
c confine rate constant to upper limit (diffusion control)
c a=k; dclim=upper limit due to diffusion-control
c factor 2.d0 is to account for larger upper limit for
c    2nd order reactions between differtly-charged ions

      implicit none

      double precision a,dclim

      dclim = 1.d10

      dmin3 = dmin1(a,dclim*2.d0 )

      end

c----------------------------------------------------------------

      double precision function flsc (a,b,c,d)
c calculate special rate function  
c after #s_364, Schmitz (1999), eq.(4) / #s_650, Schmitz (2000)
c: dio3/dt = k1*[IO3-][H+]^2[I-]^2 + k2*[IO3-][H+]^2[I-]
c a=k1, b=H+, c=I-, d=cvvz

      implicit none

      double precision a,b,c,d

      if (d.gt.0.d0) then
         flsc=( a*b**2*d**4 + 1.2d3*b**2/c*d**3 )
      else
         flsc=0.
      endif
      end

c----------------------------------------------------------------

      double precision function flsc4 (a,b,c)
c calculate special rate function 
c after #s_650, Schmitz (2000)
c a=k4, b=H+, c=cvvz

      implicit none

      double precision a,b,c

      if (c.gt.0.d0) then
         flsc4=( a*b*c**3 )
      else
         flsc4=0.
      endif
      end

c----------------------------------------------------------------

      double precision function flsc5 (a,b,c)
c calculate special rate function 
c after #s_650, Schmitz (2000)
c a=k5, b=H+, c=cvvz

      implicit none

      double precision a,b,c

      if (c.gt.0.d0) then
        flsc5=( a*b**2*c**4 )
      else
        flsc5=0.
      endif
      end

c----------------------------------------------------------------

      function flsc6 (a,b)

! Description :
! -----------
!    calculate special rate function
!    after G. Schmitz (pers.comm.)
!    a=k6, b=H+

! Author :
! ------
!    Roland von Glasow

! Modifications :
! -------------
!  17-Oct-2016  Josue Bock  implicit none (u8.5)
!
!  05-Mar-2017  Josue Bock  introduced a lower limit for [H+] to prevent flsc6
!                           to reach very high values if [H+] is very small
!                           The chosen value might need to be further adjusted
!
!  19-Oct-2017  Josue Bock  - added a test to warn only if negative values (avoid warning
!                             messages to be displayed just after initialisation)
!                           - cosmetic improvements: header, precision from module, ...
! End modifications
!-----------------------------------------------------------------------------------------------------------------------


! Declarations:
! ------------
! Modules used:
      USE precision, ONLY:
! Imported Type Definitions:
     &     dp                   ! kind double precision real

      implicit none

! Function result
      real(kind=dp) :: flsc6

! Function arguments
! Scalar arguments with intent(in):
      real(kind=dp), intent(in) :: a
      real(kind=dp), intent(in) :: b

! Local parameters:
      real(kind=dp), parameter :: min_Hp = 1.e-15_dp ! <jjb> introduce a lower limit for [H+] to avoid overflow flsc6

!- End of header ------------------------------------------------------------

      if (b.gt.min_Hp) then
         flsc6= a/b
      else
         flsc6=0._dp
         if (b.lt.0._dp) print*,"Warning: flsc6 encountered [H+]<0",b
      endif

      end function flsc6

c----------------------------------------------------------------

      double precision function uplim (a,b,c,d)

      implicit none

      double precision a,b,c,d,dclim
c dclim = upper limit for diffusion-controlled reactions
c a=k-; b=k+; alpha=b/dclim; c=H+; d=cvvz;

      dclim = 1.d10
      if (d.gt.0.d0) then
         if (c<0.) print*,"Warning uplim encountered [H+]<0" ! Track unexpected case
        !uplim = ( a/(1.+b/dclim*c*d) )
        uplim = ( a/(1.+b/dclim*max(c,0.d0)*d) ) ! <jjb> avoid unexpected values if [H+]<0
      else
        uplim = 0.
      endif
      end

c----------------------------------------------------------------

      double precision function uparm (a0,b0,c,d,e)
c Arrhenius function but with b0 as the value for T=298K
c dclim = upper limit for diffusion-controlled reactions
c c=k+; alpha=c/dclim; d=[H+]; e=cvvz

      implicit none

      double precision a0,c,d,e,dclim
      integer b0
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk
c 1/298.=3.3557d-3

      dclim = 1.d10
      if (d.gt.0.d0) then
        uparm=a0*exp(dble(b0)*(1.d0/te-3.3557d-3))/(1.+c/dclim*d*e)
      else
        uparm=0.
      endif
      end

c----------------------------------------------------------------

      double precision function uplip (a,b,c)

      implicit none

      double precision a,b,c,dclim
c dclim = upper limit for diffusion-controlled 3rd order reactions
c with H+ as reactant
c a=k+; b=[H+]; c=cvvz; alpha=a/dclim

      dclim = 1.d10
      if (c.gt.0.d0) then
         if (b<0.) print*,"Warning uplip encountered [H+]<0" ! Track unexpected case
        !uplip = ( a/(1.+a/dclim*b*c)*c**2 )
         uplip = ( a/(1.+a/dclim*max(b,0.d0)*c)*c**2 ) ! <jjb> avoid unexpected values if [H+]<0
      else
        uplip = 0.
      endif
      end

c----------------------------------------------------------------

      double precision function uparp (a0,b0,c,d)
c Arrhenius function but with b0 as the value for T=298K
c dclim = upper limit for diffusion-controlled reactions
c c=H+; d=cvvz; alpha=f(a0,b0)/dclim

      implicit none

      double precision a0,c,d,dclim
      integer b0
      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

c 1/298.=3.3557d-3

      dclim = 1.d10
      if (d.gt.0.d0) then
        uparp=a0*exp(dble(b0)*(1.d0/te-3.3557d-3))/
     &       (1.+(a0*exp(dble(b0)*(1.d0/te-3.3557d-3)))/dclim*c*d)*d**2
      else
        uparp=0.
      endif
      end

c----------------------------------------------------------------

      double precision function fhet_da (xliq,xhet,a0,b0,c0)
c heterogeneous rate function
c ClFCT     = 5.0D2                ; factor for H02/H01, i.e Cl-/H2O
c BrFCT     = 3.0D5                ; factor for H03/H01, i.e Br-/H2O
c a0=1..2  liquid size class
c b0=1..3  branch of het reaction: H2O, Cl-, Br-
c c0=1..3  gas phase reactant:     N2O5, ClNO3, BrNO3

      implicit double precision (a-h,o-z)

      INCLUDE 'aer_Parameters.h'
      INCLUDE 'aer_Global.h'
      integer a0,b0,c0


      if (xhet.eq.0.) then
         if (c0.eq.1) xtr=yxkmt(a0,ind_N2O5)
         if (c0.eq.2) xtr=yxkmt(a0,ind_ClNO3)
         if (c0.eq.3) xtr=yxkmt(a0,ind_BrNO3)
         if (a0.eq.1) then
            h2oa=FIX(indf_H2Ol1)
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycw(a0)
         endif
         if (a0.eq.2) then
            h2oa=FIX(indf_H2Ol2)
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycw(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=FIX(indf_H2Ol1)
            if (a0.eq.2) hetT=FIX(indf_H2Ol2)
         endif
      else
!        if (c0.eq.2) xtr=yxkmtd(a0,ind_N2O5)  ! jjb wrong index
!        if (c0.eq.3) xtr=yxkmtd(a0,ind_BrNO3) ! jjb wrong index
!        if (c0.eq.4) xtr=yxkmtd(a0,ind_ClNO3) ! jjb wrong index
         if (c0.eq.1) xtr=yxkmtd(a0,ind_N2O5)  ! jjb corrected
         if (c0.eq.2) xtr=yxkmtd(a0,ind_BrNO3) ! jjb corrected
         if (c0.eq.3) xtr=yxkmtd(a0,ind_ClNO3) ! jjb corrected
!         print*,xliq,a0,c0,xtr
         if (a0.eq.1) then
            h2oa=55.55*ycwd(1)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycwd(a0)
         endif
         if (a0.eq.2) then
            h2oa=55.55*ycwd(2)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycwd(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=55.55*ycwd(1)*1.d+3
            if (a0.eq.2) hetT=55.55*ycwd(2)*1.d+3
         endif
      endif

      if (b0.eq.1) xbr=h2oa
      if (b0.eq.2) xbr=5.0D2
      if (b0.eq.3) xbr=3.0D5

      if (hetT.gt.0.d0) then
         fhet_da=xtr * yw * xbr /hetT
      else
         fhet_da=0.
      endif
      if ((c0.eq.2.or.c0.eq.3.or.b0.eq.2.or.b0.eq.3).and.xhal.eq.0.) 
     &     fhet_da=0.
      if (xliq.eq.0.) fhet_da=0.
!      print*,xliq,a0,c0,fhet_da
      end function fhet_da

c----------------------------------------------------------

      double precision function fhet_dt (xliq,xhet,a0,b0,c0)
c heterogeneous rate function 
c ClFCT     = 5.0D2                ; factor for H02/H01, i.e Cl-/H2O
c BrFCT     = 3.0D5                ; factor for H03/H01, i.e Br-/H2O
c a0=1..2  liquid size class
c b0=1..3  branch of het reaction: H2O, Cl-, Br-
c c0=1..3  gas phase reactant:     N2O5, ClNO3, BrNO3

      implicit double precision (a-h,o-z)

      INCLUDE 'tot_Parameters.h'
      INCLUDE 'tot_Global.h'
      integer a0,b0,c0


      if (xhet.eq.0.) then
         if (c0.eq.1) xtr=yxkmt(a0,ind_N2O5)
         if (c0.eq.2) xtr=yxkmt(a0,ind_ClNO3)
         if (c0.eq.3) xtr=yxkmt(a0,ind_BrNO3)
         if (a0.eq.1) then
            h2oa=FIX(indf_H2Ol1)
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycw(a0)
         endif
         if (a0.eq.2) then
            h2oa=FIX(indf_H2Ol2)
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycw(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=FIX(indf_H2Ol1)
            if (a0.eq.2) hetT=FIX(indf_H2Ol2)
         endif
      else
!        if (c0.eq.2) xtr=yxkmtd(a0,ind_N2O5)  ! jjb wrong index
!        if (c0.eq.3) xtr=yxkmtd(a0,ind_BrNO3) ! jjb wrong index
!        if (c0.eq.4) xtr=yxkmtd(a0,ind_ClNO3) ! jjb wrong index
         if (c0.eq.1) xtr=yxkmtd(a0,ind_N2O5)  ! jjb corrected
         if (c0.eq.2) xtr=yxkmtd(a0,ind_BrNO3) ! jjb corrected
         if (c0.eq.3) xtr=yxkmtd(a0,ind_ClNO3) ! jjb corrected
         if (a0.eq.1) then
            h2oa=55.55*ycwd(1)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml1) + 3.0D5*C(ind_Brml1)
            yw=ycwd(a0)
         endif
         if (a0.eq.2) then
            h2oa=55.55*ycwd(2)*1.d+3
            hetT=h2oa + 5.0D2*C(ind_Clml2) + 3.0D5*C(ind_Brml2)
            yw=ycwd(a0)
         endif
         if (xhal.eq.0.) then
            if (c0.eq.2) xtr=0.
            if (c0.eq.3) xtr=0.
            if (a0.eq.1) hetT=55.55*ycwd(1)*1.d+3
            if (a0.eq.2) hetT=55.55*ycwd(2)*1.d+3
         endif
      endif

      if (b0.eq.1) xbr=h2oa
      if (b0.eq.2) xbr=5.0D2
      if (b0.eq.3) xbr=3.0D5

      if (hetT.gt.0.d0) then
         fhet_dt=xtr * yw * xbr /hetT
      else
         fhet_dt=0.
      endif
      if ((c0.eq.2.or.c0.eq.3.or.b0.eq.2.or.b0.eq.3).and.xhal.eq.0.) 
     &     fhet_dt=0.
      if (xliq.eq.0.) fhet_dt=0.

      end function fhet_dt

c----------------------------------------------------------

      double precision function fdhetg (na,nb)
c heterogeneous rate function 
c a0=1..2  liquid size class

! na is bin number
! nb is reaction:
!   1 HNO3
!   2 N2O5
!   3 NH3
!   4 H2SO4

      implicit double precision (a-h,o-z)

      INCLUDE 'gas_Parameters.h'
      INCLUDE 'gas_Global.h'
      integer na,nb
c net mass transfer coefficient including Henry's law equilibrium for HNO3

      if (nb.eq.1) then
c not limited by Henry's law:
c         xkt=yxkmtd(na,ind_HNO3) * ycwd(na) 

c limited by Henry's law:
c see Diss RvG (3.10): dcg/dt = ... + kmt(LWC*Cg - Ca/H), this term is implemented here
c note that Cg is multiplied to rate in KPP, therefore the heterogeneous reaction rate is:
c     kmt(LWC - Ca/(Cg*H)) = x1 + x2
c     in x2 the "aqueous" concentration of HNO3 on the dry aerosol is calculated via 
c     Henry's law and a HARDCODED particle pH = 2
     
         x1 = yxkmtd(na,ind_HNO3) * ycwd(na) 
c index out of bounds in C(ind_NO3mlz) as this is not known in gas_Parameters.h
c         if (na.eq.1) caq=((C(ind_HNO3l1)+C(ind_NO3ml1))*1.d-2)/
c     &        (yxeq(ind_HNO3) + 1.d-2)
c         if (na.eq.2) caq=((C(ind_HNO3l2)+C(ind_NO3ml2))*1.d-2)/
c     &        (yxeq(ind_HNO3) + 1.d-2)
c if pH=2, the fractionation between HNO3 and NO3- can be calculated with equil. const:
c [NO3-]=Kq [HNO3] 1/[H+] = 1500 [HNO3] at pH=2
c this is all VERY rough and should be replaced!!
         if (na.eq.1) caq=((C(ind_HNO3l1)*1.5d3)*1.d-2)/
     &        (yxeq(ind_HNO3) + 1.d-2)
         if (na.eq.2) caq=((C(ind_HNO3l2)*1.5d3)*1.d-2)/
     &        (yxeq(ind_HNO3) + 1.d-2)
         x2 = 0.d0
          if (C(ind_HNO3).ne.0.d0.and.yhenry(ind_HNO3).ne.0.d0) 
     &        x2=-yxkmtd(na,ind_HNO3)/(C(ind_HNO3)*yhenry(ind_HNO3))*caq
         xkt = max(0.d0,(x1 + x2))
      endif
c kmt in (m^3_air/(m^3_aq*s)) therefore multiplication with LWC (m^3_aq/m^3_air)
c to get k in 1/s:
      if (nb.eq.2) xkt=yxkmtd(na,ind_N2O5) * ycwd(na) 
      if (nb.eq.3) xkt=yxkmtd(na,ind_NH3) * ycwd(na) 
      if (nb.eq.4) xkt=yxkmtd(na,ind_H2SO4) * ycwd(na) 


      fdhetg=xkt

c      print *,k,na,nb
c      print *,fdhetg,xkt,ycwd(na)
c      if (nb.eq.1) write (440, 1001) k,na,nb,fdhetg,yxkmtd(na,ind_HNO3) 
c     &     * ycwd(na)
! 1001 format(3i4, 6d16.8)

      end function fdhetg

c----------------------------------------------------------------

      double precision function fdheta (na,nb)
c heterogeneous rate function 
c a0=1..2  liquid size class

      implicit double precision (a-h,o-z)

      INCLUDE 'aer_Parameters.h'
      INCLUDE 'aer_Global.h'
      integer na,nb
c see explanation in FCN fdhetg

      if (nb.eq.1) then
         x1 = yxkmtd(na,ind_HNO3) * ycwd(na) 
         caq = 0.d0
         if ((yxeq(ind_HNO3)+1.d-2).ne.0.d0) then
            if (na.eq.1) caq=((C(ind_HNO3l1)+C(ind_NO3ml1))*1.d-2)/
     &           (yxeq(ind_HNO3) + 1.d-2)
            if (na.eq.2) caq=((C(ind_HNO3l2)+C(ind_NO3ml2))*1.d-2)/
     &           (yxeq(ind_HNO3) + 1.d-2)
         endif
         x2 = 0.d0
          if (C(ind_HNO3).ne.0.d0.and.yhenry(ind_HNO3).ne.0.d0) 
     &        x2=-yxkmtd(na,ind_HNO3)/(C(ind_HNO3)*yhenry(ind_HNO3))*caq
         xkt = max(0.d0,(x1 + x2))
      endif
      if (nb.eq.2) xkt=yxkmtd(na,ind_N2O5) * ycwd(na) 
      if (nb.eq.3) xkt=yxkmtd(na,ind_NH3) * ycwd(na) 
      if (nb.eq.4) xkt=yxkmtd(na,ind_H2SO4) * ycwd(na) 

      fdheta=xkt 

      end function fdheta

c-----------------------------------------------------------------------------


      double precision function fdhett (na,nb)
c heterogeneous rate function 
c a0=1..2  liquid size class

      implicit double precision (a-h,o-z)

      INCLUDE 'tot_Parameters.h'
      INCLUDE 'tot_Global.h'
      integer na,nb
c see explanation in FCN fdhetg

      if (nb.eq.1) then
         x1 = yxkmtd(na,ind_HNO3) * ycwd(na) 
         caq = 0.d0
         if ((yxeq(ind_HNO3)+1.d-2).ne.0.d0) then
            if (na.eq.1) caq=((C(ind_HNO3l1)+C(ind_NO3ml1))*1.d-2)/
     &           (yxeq(ind_HNO3) + 1.d-2)
            if (na.eq.2) caq=((C(ind_HNO3l2)+C(ind_NO3ml2))*1.d-2)/
     &           (yxeq(ind_HNO3) + 1.d-2)
         endif
         x2 = 0.d0
          if (C(ind_HNO3).ne.0.d0.and.yhenry(ind_HNO3).ne.0.d0) 
     &        x2=-yxkmtd(na,ind_HNO3)/(C(ind_HNO3)*yhenry(ind_HNO3))*caq
         xkt = max(0.d0,(x1 + x2))
      endif
      if (nb.eq.2) xkt=yxkmtd(na,ind_N2O5) * ycwd(na) 
      if (nb.eq.3) xkt=yxkmtd(na,ind_NH3) * ycwd(na) 
      if (nb.eq.4) xkt=yxkmtd(na,ind_H2SO4) * ycwd(na) 

      fdhett=xkt 

      end function fdhett

c-----------------------------------------------------------------------------
!     double precision function DMS_add (c) ! jjb argument not used
      double precision function DMS_add ()  ! jjb removed (also in mech/master_gas.eqn)
c calculate special rate function for DMS + OH addition; IUPAC 10/06
c k(298K)=2.2d-12 cm3/(mlc s)
      implicit none

      double precision o2,tte

      common /cb_1/ aircc,te,h2oppm,pk
      double precision aircc,te,h2oppm,pk

      o2=0.21*aircc
      tte=1./te
      DMS_add=9.5d-39*exp(5270.*tte)*o2/(1.+7.5d-29*exp(5610.*tte)*o2)
      end function DMS_add

c-----------------------------------------------------------------------------

c      double precision function xkHgBr (x1)
cc rate coefficient for recmonbination Hg+Br --> HgBr (Donohoue et al., 2006, #4161
c      double precision aircc,te,h2oppm,pk,x1,x2
c      common /cb_1/ aircc,te,h2oppm,pk
c
cc reaction is 3rd order, multiply with conversion factors here (instead of in 
cc master_gas.eqn) in order to have an argument for the function
cc note that the fit that they give does not exactly reproduce the measured values in their Tables 1 and 2)
c      xkHgBr = 1.46d-32 * (te/298.)**(-1.86) * x1 * x1
c      end function xkHgBr

c-----------------------------------------------------------------------------

c      double precision function xkHgBrBr (x1)
cc rate coefficient for recombination HgBr+Br --> HgBr2 (Goodsite et al., 2004, #3244
c      double precision aircc,te,h2oppm,pk,x1
c      common /cb_1/ aircc,te,h2oppm,pk
c
cc reaction is 2nd order, multiply with conversion factor here (instead of in 
cc master_gas.eqn) in order to have an argument for the function
c      xkHgBrBr = 2.5d-10*(te/298.d0)**(-.57d0) * x1
c      end function xkHgBrBr

c-----------------------------------------------------------------------------

c      double precision function xkGood (x1)
cc rate coefficient for HgBr dissociation, use Goodsite et al., 2004
cc (#3244) but scaled with ratio of Donohoue and Goodsite as
cc suggested in Seigneur and Lohmann, 2008 (#4136)
c      double precision aircc,te,h2oppm,pk,x1,x2,xkHgBr
c      common /cb_1/ aircc,te,h2oppm,pk
c      
c      if (te.eq.0.d0) then
c         xkGood = 0.d0
c      else
cc argument for kHgBr is unity to avoid scaling with conversion factors
c         x2=xkHgBr(1.d0) / (1.1d-12 * (te/298.d0) **(-2.37))
cc multiply x2 with air density as Donohoue is 3rd order and Goodsite 2nd order
cc [M in mol/m3] = p / RT; convert to molec/cm3 : * conv1
c         x2 = x2 * pk / (8.314 * te) * x1
c         xkGood = 1.2d10 * exp(-8357.d0 / te) * x2
c      endif
c
c      end function xkGood

c-----------------------------------------------------------------------------



