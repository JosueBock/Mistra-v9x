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


! MISTRA with chemistry
!
! str.f : main, meteo-init, microphysics, turbulence

! box = .true. is option to run model in one level only
!     so far without dynamics, microphysics and radiation.
!     Only the photolysis rates are updated and T, rh, .. 
!     are changed with sinus function.
!     Init is the same as in 1D run, to make sure that all
!     variables are defined and to avoid too much differences
!     between the box and the 1D version.

!
! INITIALISATION
!
! ... no restart
!       |
!       |_____SR vgleich
!       |      |_____FN rgl

! ... if (mic)
!       |_____SR difp
!       |_____SR kon
!       |      |    _SR equil (case 1)   ! if rH < 70%
!       |      |   /  |_____FN rgl
!       |      |__/
!       |      |  \
!       |      |   \_SR subkon           ! if rH >= 70%
!       |      |      |_____SR advec
!       |      |
!       |     ... if (chem)
!       |      |_____SR konc
!       |
!       |_____SR sedp
!       |_____SR equil (case 2)
!
!       


      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nm,
     &     nrlay,
     &     nka,
     &     nkt,
     &     nphrxn,
     &     mbs

      implicit double precision (a-h,o-z)

      logical chem,mic,rst,halo,iod,box,netCDF,BL_box,nuc,binout
      logical Napari, Lovejoy, both

      common /cb16/ u0,albedo(mbs),thk(nrlay)
      double precision u0, albedo, thk

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /hall/ halo,iod
      common /band_rat/ photol_j(nphrxn,n)
      common /band_o3/ scaleo3_m
      common /liq_pl/ nkc_l
      common /kpp_eul/ xadv(10),nspec(10),neula
      common /nucfeed/ ifeed


      dimension aer(n,nka)
      character *1 fogtype
      character *10 fname
! options for program evaluation
! chem     : chemistry included
! mic      : microphysics included
! rst      : restart or initialization of program run
! netCDF   : output in netCDF format
! binout   : output in binary format
! box      : box model run
! halo     : halogen chemistry on/off
! iod      : iodine chemistry on/off
! nuc      : nucleation on/off
! fogtype  : suffix for output filenames
! iaertyp  : type of aerosol; 1=urban, 2=rural, 3=ocean, 4=background
! lstmax   : integration time in hours
! scaleo3_m: total O3 in DU (for photolysis only)
! nkc_l    : number of output classes for aq. chem.
! neula    : eulerian (0) or lagrangian (1) view
! z_box    : height of MBL (if box run)
! BL_box   : box only, average init cond over BL and/or mean of J-values over BL
! nlevbox  : box only, level to be used for init cond of box if  BL_box=false
      open (11,file='istart', status='old')
      read (11,5000) chem
      read (11,5000) mic
      read (11,5000) rst
      read (11,5000) netCDF
      read (11,5000) binout
      read (11,5000) box
      read (11,5000) halo
      read (11,5000) iod
      read (11,5000) nuc
      fogtype='a'
      read (11,5020) iaertyp
      read (11,5030) lstmax
      read (11,5040) scaleo3_m
      read (11,5020) nkc_l
      read (11,5020) neula
      read (11,5050) z_box
      read (11,5000) BL_box
      read (11,5030) nlevbox
      close (11)
      print *,chem,mic,rst
      print *,netCDF,binout,box
      print *,halo,iod,nuc
      print *,fogtype,iaertyp,lstmax
      print *,scaleo3_m,nkc_l,neula
      print *,z_box,BL_box,nlevbox
      ifeed = 1
      if (neula.eq.0) then
         open (12,file='euler_in.dat',status='old')
         do i=1,10
            read (12,5100) nspec(i),xadv(i)
         enddo
         close (12)
      endif
 5000 format(l1)
! 5010 format(a1)
 5020 format (i1)
 5030 format (i3)
 5040 format (f3.0)
 5050 format (f4.0)
 5100 format (i3,d8.2)

      call mk_interface

      Napari = .true. ; Lovejoy = .true.
      if (nuc) call nuc_init(Napari,Lovejoy,iod)

      if (box) print *,'box model run'
! it's important to keep this n_bl = 2 for box runs as loops are designed that way
! (especially output)
      if (box) then
         n_bl  = 2
         n_bln = 2
         n_bl8 = 1
      else
         n_bl  = nf
         n_bln = n
         n_bl8 = 15
      endif
! open input/output files
      call openm (fogtype)
      call openc (fogtype,nuc)
! netCDF output
      if (netCDF) call open_netcdf(n_bln,chem,mic,halo,iod,nuc)
! numerical gridpoints
      call grid
      nz_box = 0
      if (box) call get_n_box (z_box,nz_box)
      call write_grid ! writes information on grid that is not f(t)
! radiation parameters qabs, qext and asym
      call intrad
      dt = 60. ! jjb moved from below, was missing for restart case
      if (rst) go to 2000

! Continue the initialisation, no-restart case
! --------------------------------------------
! initial meteorological and chemical input
      call initm (iaertyp,fogtype)
      call initc(box,n_bl)
! number of iterations
      it0=0
      itmax=60*lstmax
! initial exchange coefficients and turbulent kinetic energy
      call atk0
! initial position of humidified aerosols
      call vgleich
! first radiation calculation
      call initstr
! output of meteorological and chemical constants of current run
      call constm (chem,mic,rst)
      if (chem) call constc
! output of initial vertical profiles
      call profm (0.d0)
      if (chem) call profc (0.d0,mic)
      go to 2010

! Continue the initialisation, restart case
! -----------------------------------------
! read meteorological and chemical input from output of previous run
 2000 call startm (fogtype)
! init some microphysical data that's not in SR startm
      call mic_init (iaertyp,fogtype)
!+      it0=it    ! use when time stamp from restart run is to be preserved
      it0=0
      it=0
      if (chem) call startc (fogtype)

! initialization of radiation code
      print*,'initialisation, call initstr (restart)'
      call initstr
! number of iterations
      itmax=it0+60*lstmax
! output of meteorological and chemical constants of current run
      call constm (chem,mic,rst)
! output of initial profiles of restart run
      call profm (dt)
      if (chem) call profc (dt,mic)
! allocate arrays and initialise vmean
      if (chem) call v_mean_init
      if (chem) call v_mean (t(:nf))


! Continue the initialisation, both cases
! ---------------------------------------
! initial photolysis rates
!      print*,'initialisation, call str'
 2010 call str (mic)
!      if (chem) call photol
      if (chem) then
         call photol_initialize
         call photol
      end if

! initial output for plotting
      if (binout) then
         call ploutm (fogtype,n_bln)
         if (mic.and..not.box) call ploutp (fogtype)
         call ploutr (fogtype,n_bln)
         call ploutt (fogtype,n_bln)
         if (chem) call ploutc (fogtype,mic,n_bl,n_bl8)
         if (chem) call ploutj (fogtype,n_bln)
      endif
      if (chem) call out_mass
      if (netCDF) call write_netcdf(n_bln,chem,mic,halo,iod,box,nuc)
      time=60.*float(it0)
! local time: day (lday), hours (lst), minutes (lmin)
      fname='tim .out'
      fname(4:4)=fogtype
      open (99, file=fname,status='unknown',err=2005)
      atmax=0.
      write (99,6000) lday,lst,lmin,atmax
      close (99)
 2005 continue
      if (box) call box_init (nlevbox,nz_box,n_bl,BL_box)
      if (box) box_switch=1.
      print*,'end initialisation str.f'
! ====================integration in time=====================
! outer time loop: minutes
      do 1000 it=it0+1,itmax                    
         if (lct.gt.nf) stop 'lct.gt.nf'
!         time=time+dt
         lmin=lmin+1
         if (lmin.lt.60) go to 2030
         lmin=lmin-60
         lst=lst+1
         if (lst.eq.24) then
            lst=0
            lday=lday+1
         endif
 2030    continue
! dry dep velocities
!         print*,'call partdep'
         call partdep (xra)
! dd: fractional timestep in sec
         dd=10.
! inner time loop: 10 sec
         do ij=1,6
            time=time+dd
! --------1D only start------------
! skip dynamics, microphysics and radiation for box model run
            if (.not.box) then
! if w-field variable in time call wfield
!         call wfield
! turbulent exchange of thermodynamic variables, particles and
! chemical species
!         print*,'call difm'
               call difm (dd)
               if (chem) call difc (dd)
! microphysics
               if (mic) then
                  call difp (dd)
! condensation/evaporation, update of chemical concentrations
!         print*,'call kon'
                  call kon (dd,chem)
! gravitational settling of particles
!         print*,'call sedp'
                  call sedp (dd)
! put aerosol into equilibrium with current rel hum for k>nf
!         print*,'call equil'
                  call equil (2,k)
               endif
! put aerosol into equilibrium with current rel hum 
               if (.not.mic) call equil (1,n_bl)
! radiative heating
               do k=2,nm
                  t(k)=t(k)+dtrad(k)*dd
               enddo
! temperature and humidity within the soil
! water surface: no call to soil
!         call soil (dd)
! flux balances at the earth's surface
! water surface: call surf0; else: call surf1
               call surf0 (dd)
!         call surf1 (dd)
! dry deposition and emission of chemical species
               if (chem) then
                  call sedc (dd)
! wet deposition of chemical species
!         if (lct.gt.1) call sedl (dd)
                  call sedl (dd)
! chemical reactions
!                 call stem_kpp (dd,xra,z_box,n_bl,box)     ! jjb
                  call stem_kpp (dd,xra,z_box,n_bl,box,nuc) ! jjb nuc is needed in this SR
                  if (nuc) then
! set switches for ternary nucleation: Napari: ternary H2SO4-H2O-NH3 nucleation
!                                      Lovejoy: homogeneous OIO nucleation
!                                      for further explanation see nuc.f
                    !Napari = .true.
                    !Lovejoy = .true.                 ! <jjb> defined previously
                    if ((Napari) .and. (Lovejoy)) then
                      both = .true.
                    else
                      both = .false.
                    endif
                    if ((.not.Napari) .and. (.not.Lovejoy)) 
     $                STOP 'Napari or Lovejoy must be true'
                    if (both) then
!         print*,'call appnucl2'
                      call appnucl2 (dd,both)
                    else
!         print*,'call appnucl'
                      call appnucl (dd,Napari,Lovejoy,both)
                    endif
!                   --- de-comment if .asc output for gnu-plotting is desired ---
!                    call nucout1
                  endif
               endif
            else                ! if .not.box
! --------1D only end------------
! --------box model version only start -------------------
! put aerosol into equilibrium with current rel hum 
!               if (.not.mic) call equil (1,n_bl)
! call u0, T, rh, J-values  .. update
               call box_update(box_switch,ij,nlevbox,nz_box,n_bl,
!     &              chem,halo,iod,BL_box) ! jjb 3 unused arguments
     &              BL_box)
! gas phase emissions and deposition
               call sedc_box (dd,z_box,n_bl)
! particle and aqueous phase deposition
               call box_partdep (dd,z_box,n_bl)
! aerosol emission and chemical reactions
               if (chem) then
!                 call stem_kpp (dd,xra,z_box,n_bl,box)     ! jjb
                  call stem_kpp (dd,xra,z_box,n_bl,box,nuc) ! jjb nuc is needed in this SR
               endif
            endif               ! if .not.box
! --------box model version only end -------------------
         enddo                  ! ij-loop : end of fractional timestep loop

! radiative fluxes and heating rates
!         print*,'call str'
         if (.not.box) call str (mic)
! new photolysis rates
         if (chem) then
            if (u0.gt.3.48e-2) then
! optimize this!!
               if (u0.gt.3.48e-2.and.u0.le.0.4.and.lmin/2*2.eq.lmin
     &             .or.u0.gt.0.4.and.lmin/2*2.eq.lmin) call photol
               if (box.and.BL_box) call ave_j (nz_box,n_bl)
            else
               do k=1,n
                  do i=1,nphrxn ! jjb
                     photol_j(i,k)=0.
                  enddo
               enddo
            endif
         endif

! output of meteorological and chemical variables ----------------------
         ilmin=15
!         ilmin=1 !output every minute
         if (lmin/ilmin*ilmin.eq.lmin) then
! calc 1D size distribution for output
!         print*,'call oneD_dist'
            call oneD_dist
! binary output
            if (binout) then 
               call ploutm (fogtype,n_bln)
               if (lmin/30*30.eq.lmin.and.mic.and..not.box) 
     &              call ploutp (fogtype)
               call ploutr (fogtype,n_bln)
               call ploutt (fogtype,n_bln)
               if (chem) call ploutc (fogtype,mic,n_bl,n_bl8)
            endif
! netCDF output
            if (netCDF) call write_netcdf(n_bln,chem,mic,halo,iod,
     &           box,nuc)
! output of data from nucleation
           if (chem.and.nuc) call nucout2
! output from mass balance
            if (chem) call out_mass
!         if (chem.and.lmin/60*60.eq.lmin) call ploutj(fogtype,n_bln)
         endif
! hourly output of profiles in ascii files
         if (lmin/60*60.eq.lmin) then 
            call profm (dt)
!           call profr
            if (chem) call profc (dt,mic)
         endif
! output for restart option
!     comment these calls to save disk space for production runs:
         if (lst/12*12.eq.lst.and..not.box.and.lmin.eq.0) then
            call outm
            if (chem) call outc
         endif
! output of "tima.out"
         atmax=0.
         tkemax=0.
         xm2max=0.
         do k=lcl,nf
            atmax=dmax1(atmax,atkh(k))
            tkemax=dmax1(tkemax,tke(k))
            xm2max=dmax1(xm2max,xm2(k)*1000./rho(k))
         enddo
         open (99, file=fname,status='unknown',err=1000)
         write (99,6010) lday,lst,lmin,
     &        tkemax,atmax,xm2max,eta(lcl),eta(lct)
         write (*,6010) lday,lst,lmin,
     &        tkemax,atmax,xm2max,eta(lcl),eta(lct)
 6000    format (' time: ',i2,':',i2,':',i2,3x,' iteration: ',f10.3,3x,
     &        'cloudy region: ',f7.1,' - ',f7.1)
 6010    format (1x,i2,':',i2,':',i2,3f10.3,3x,
     &        'cloudy region: ',f7.1,' - ',f7.1)
         close (99)
 1000 continue          
! =========================end of time integration=====================


! final output of restart files
      call outm
      if (chem) call outc
! final output of aerosol size distribution
      do k=1,n
         do ia=1,nka
            aer(k,ia)=0.
            do jt=1,nkt
               aer(k,ia)=aer(k,ia)+ff(jt,ia,k)
            enddo
         enddo
      enddo
      fname='ae .out'
      fname(3:3)=fogtype
      open (66, file=fname,status='unknown',form='unformatted')
      write (66) aer
      close (66)

      if (netCDF) call close_netcdf(mic,chem,nuc)

      stop 'main program'
      end

!
!-----------------------------------------------------------------------
!

      block data
! defines parameters that are accessible to all subroutines

      USE global_params, ONLY :
! Imported Parameters:
     &     nka

      implicit double precision (a-h,o-z)

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

! gravitational acceleration
      data g /9.8065d0/

! chose the water temperature and subsidence velocities depending on 
! what version of SR initm is used (see ./special_versions/SR_initm)
! water temperature
!      data tw /288.15d0/ !cloud and aer sub run
!      data tw /286.15d0/ !cloud and aer sub run
!      data tw /287.4d0/ !cloud no sub run
!      data tw /290.4d0/ !aerosol no sub run
!      data tw /288.4d0/ !aerosol no sub run
!      data tw /300.d0/ !INDOEX
      data tw /299.5d0/

! geostrophic wind, large scale subsidence
!      data ug,vg,wmin,wmax / 6.d0, 0.d0, 0.d0,-0.005d0/
!      data ug,vg,wmin,wmax / 6.d0, 0.d0, 0.d0,-0.004d0/
!      data ug,vg,wmin,wmax /10.d0, 0.d0, 0.d0,-0.004d0/
!      data ug,vg,wmin,wmax /10.d0, 0.d0, 0.d0, 0.d0/
!      data ug,vg,wmin,wmax /8.5d0, 0.d0, 0.d0, 0.d0/
!      data ug,vg,wmin,wmax /7.0d0, 0.d0, 0.d0, 0.d0/
!      data ug,vg,wmin,wmax / 7.d0, 0.d0, 0.d0,-0.006d0/   !cloud sub
!      data ug,vg,wmin,wmax /8.5d0, 0.d0, 0.d0,-0.0015d0/ !aerosol sub
!      data ug,vg,wmin,wmax / 6.d0, 0.d0, 0.d0, 0.d0/
!      data ug,vg,wmin,wmax / 6.d0, 0.d0, 0.d0, 0.d0/
!      data ug,vg,wmin,wmax / 6.d0, 0.d0, 0.d0,-0.01d0/
!      data ug,vg,wmin,wmax / 6.d0, 0.d0, 0.d0,-0.02d0/
!       data ug,vg,wmin,wmax /8.0d0, 0.d0, 0.d0, 0.d0/
!       data ug,vg,wmin,wmax /15.0d0, 0.d0, 0.d0, 0.d0/
!       data ug,vg,wmin,wmax /15.0d0, 0.d0, 0.d0, -0.0015d0/ !aerosol sub (value copied from above)
       data ug,vg,wmin,wmax /15.0d0, 0.d0, 0.d0, -0.006d0/ !cloud sub (value copied from above)
! surface roughness
!      data z0 /0.01d0/
      data z0 /0.00001d0/ 
! soil constants for sandy loam
      data ebs,psis,aks,bs,rhoc /.435d0,-.218d0,3.41d-05,4.9d0,1.34d+06/
      data rhocw,ebc,anu0,bs0
     &     /4.186d+06,.0742724d0,43.415524d0,2.128043d0/
      end block data

!
!-------------------------------------------------------------
!

      subroutine openm (fogtype)
! input/output files

      USE directories, ONLY : inpdir

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      character (len=1) fogtype

! Local scalars:
      character (len=10) fname    ! I/O files names
      integer i,k                 ! implied do loops indexes

! Common blocks:
      common /cb61/ fu(18,7),ft(18,7),xzpdl(18),xzpdz0(7) ! Clarks table data
      double precision fu, ft, xzpdl, xzpdz0
!- End of header ---------------------------------------------------------------


! input Clarke-tables
      open (50, file=trim(inpdir)//'clark.dat', status='old')
      read (50,5000) ((fu(i,k),i=1,9),(fu(i,k),i=10,18),k=1,7)
      read (50,5000) ((ft(i,k),i=1,9),(ft(i,k),i=10,18),k=1,7)
      read (50,5010) (xzpdl(i),i=1,9),(xzpdl(i),i=10,18)
      read (50,5020) (xzpdz0(i),i=1,7)
      close (50)
 5000 format (9f8.4)
 5010 format (9f5.2)
 5020 format (7f5.0)
! output vertical profiles of meteorological variables
      fname='profm .out'
      fname(6:6)=fogtype
      open (26, file=fname,status='unknown')

! output plot files for meteorological data
      fname='pm .out'
      fname(3:3)=fogtype
      open (17, file=fname,status='unknown',form='unformatted')
      close (17)
      fname(2:2)='t'
      open (18, file=fname,status='unknown',form='unformatted')
      close (18)
      fname(2:2)='b'
      open (19, file=fname,status='unknown',form='unformatted')
!      close (19)
      fname(2:2)='r'
      open (14, file=fname,status='unknown',form='unformatted')
      close (14)

! output radiative fluxes and heating rates
      fname='profr .out'
      fname(6:6)=fogtype
      open (40, file=fname,status='unknown')
      fname='f1 .out'
      fname(3:3)=fogtype
      open (41, file=fname,status='unknown',form='unformatted')
      close (41)
      fname(2:2)='2'
      open (42, file=fname,status='unknown',form='unformatted')
      close (42)
      fname(2:2)='3'
      open (43, file=fname,status='unknown',form='unformatted')
      close (43)

      end subroutine openm

!
!-------------------------------------------------------------
!

      subroutine openc (fogtype,nuc)
! input/output files of chemical species
      character *10 fname
      character *1 fogtype
      logical nuc
! nucleation output:
      if (nuc) then
! de-comment if .asc output for gnu-plotting is desired
!        open (21,file='nuc+part.asc',status='unknown',form='formatted')
!        open (22,file='ternuc.asc',status='unknown',form='formatted')
!        open (23,file='nh3_h2so4.asc',status='unknown',form='formatted')
!        open (24,file='NUCV.asc',status='unknown',form='formatted')

! un-comment if the output from SR nucout2 is desired
!        open (unit=25, file='nucout2.out', status='unknown',
!     &        form='unformatted')
      endif
! chemical concentrations for initialization
!      fname='initc .dat'
!      fname(6:6)=fogtype
!      open (10,file=fname,status='old')
! vertical profiles of chemical species
      fname='profc .out'
      fname(6:6)=fogtype
      open (60,file=fname,status='unknown')
! all plotfiles for chemical species
      fname='sg1 .out'
      fname(4:4)=fogtype
      open (61,file=fname,status='unknown',form='unformatted')
      fname='sl1 .out'
      fname(4:4)=fogtype
      open (62,file=fname,status='unknown',form='unformatted')
      close (62)
      fname='ion .out'
      fname(4:4)=fogtype
      open (63,file=fname,status='unknown',form='unformatted')
      close (63)
      fname='sr1 .out'
      fname(4:4)=fogtype
      open (64,file=fname,status='unknown',form='unformatted')
      fname='ara .out'
      fname(4:4)=fogtype
      open (65,file=fname,status='unknown',form='unformatted')
      close (65)
      fname='gr .out'
      fname(3:3)=fogtype
      open (66,file=fname,status='unknown',form='unformatted')
      close (66)
      fname='gs .out'
      fname(3:3)=fogtype
      open (67,file=fname,status='unknown',form='unformatted')
      close (67)
      fname='jra .out'
      fname(4:4)=fogtype
      open (69,file=fname,status='unknown',form='unformatted')
      close (69)
!      open (64,file=fname,status='unknown',form='unformatted')
      fname='sle .out'
      fname(4:4)=fogtype
      open (67,file=fname,status='unknown',form='unformatted')
      close (67)
      open (74,file='mass.out',status='unknown')
      write (74,101)
      write (74,102)
      write (74,103)
      close (74)
 101  format ('output of molecule burden/deposit/source; unit is',
     & ' [mol/m2]')
 102  format ('to get balance: divide last output by first; to get',
     & ' emitted salt mass (in [g/m2])')
 103  format ('multiply xnass with 68.108 (=23 g(Na)/mol(Na) / 0.3377',
     & ' g(Na)/g(seasalt))')

      end subroutine openc

!
!-------------------------------------------------------------
!

      subroutine initm (iaertyp,fogtype) !change also SR surf0 !_aerosol_nosub

      USE constants, ONLY :
! Imported Parameters:
!     &     pi,
     &     r0,                   ! Specific gas constant of dry air, in J/(kg.K)
     &     r1,                   ! Specific gas constant of water vapour, in J/(kg.K)
     &     rhow                  ! Water density [kg/m**3]

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nb,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)
! initial profiles of meteorological variables
      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb51/ dlgew,dlgenw,dlne
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb53a/ thet(n),theti(n)
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /cb63/ fcs(nka),xmol3(nka)
      common /kinv_i/ kinv
!     dimension wn(4,3),wr(4,3),ws(4,3),sr(nka,nkt),aer(nf,nka), ! jjb AER not used, removed
!    &          fnorm(n)
 !     dimension wn(4,3),wr(4,3),ws(4,3),sr(nka,nkt),fnorm(n) ! jjb removed ! jjb fnorm unused
      dimension wn(4,3),wr(4,3),ws(4,3),sr(nka,nkt)     ! jjb removed
      character *1 fogtype
      character *10 fname
      data xmol2 /18./
! p21 after magnus formula
      p21(tt)=610.7*dexp(17.15*(tt-273.15)/(tt-38.33))
! constants for aerosol distributions after jaenicke (1988)
! 3 modes j=1,2,3
! 4 aerosol types i=iaertyp: 1=urban; 2=rural; 3=ocean; 4=background
!      data ((wn(i,j),i=1,4),j=1,3)/1.6169d+05,1.1791d+04,80.76,79.788,
!     & 664.9,105.29,126.52,94.138,4.3091d+04,2.9846d+03,3.0827,0.0596/
!      data ((wr(i,j),i=1,4),j=1,3)/6.51d-03,7.39d-03,3.9d-03,3.6d-03,
!     & 7.14d-03,.0269,.133,.127,.0248,.0419,.29,.259/
!      data ((ws(i,j),i=1,4),j=1,3)/8.3299,9.8765,1.1583,1.2019,
!     & 1.1273,1.6116,11.338,7.8114,4.4026,7.0665,3.1885,2.7682/
! constants for aerosol distributions after jaenicke (1988)
! except constants for maritime aerosol distribution 
! after Hoppel et al. 1990 JGR 95, pp. 3659-3686
      data ((wn(i,j),i=1,4),j=1,3)
     & /1.6169d+05,1.1791d+04,159.576,79.788,
     & 664.9,105.29,427.438,94.138,
     & 4.3091d+04,2.9846d+03,5.322,0.0596/
      data ((wr(i,j),i=1,4),j=1,3)
     & /6.51d-03,7.39d-03,0.027,3.6d-03,
     & 7.14d-03,.0269,.105,.127,
     & .0248,.0419,.12,.259/
      data ((ws(i,j),i=1,4),j=1,3)
     & /8.3299,9.8765,8.,1.2019,
     & 1.1273,1.6116,39.86,7.8114,
     & 4.4026,7.0665,2.469,2.7682/
!c aerosol distribution; f=dfdlogr*dlogr=dfdlogr*dlgenw/3
      dfdlogr(rr,ka)=wn(ka,1)*dexp(-ws(ka,1)*dlog10(rr/wr(ka,1))**2)+
     &               wn(ka,2)*dexp(-ws(ka,2)*dlog10(rr/wr(ka,2))**2)+
     &               wn(ka,3)*dexp(-ws(ka,3)*dlog10(rr/wr(ka,3))**2)
      dfdlogr2(rr,ka)=wn(ka,1)*dexp(-ws(ka,1)*dlog10(rr/wr(ka,1))**2)+
     &                wn(ka,2)*dexp(-ws(ka,2)*dlog10(rr/wr(ka,2))**2)
! after Jaenicke/Sander/Kim:
!      dfdlogr(rr,ka)=2.8d2/(0.1106*sqrt(2*pi))*dexp(-dlog10(rr/8.8d-2)
!     & **2/(2*0.1106**2))+
!     &               6.6d-1/(0.1906*sqrt(2*pi))*dexp(-dlog10(rr/1.7d0)
!     & **2/(2*0.1906**2))
! see below: call adjust_f
      lcl=1
      lct=1
! declination of the sun
      declin=18.65
!      declin=0.
! geographical latitude
!      alat=-6.44
      alat=-5.8544
! starting time of calculations
      lday=0
      lst=0
      lmin=0
! initial inversion height
!      zinv=700.
!      zinv=500.
!      zinv=900.
!      zinv=600.
      zinv=800.
      do k=2,nf
         if (eta(k).lt.zinv.and.eta(k+1).gt.zinv) kinv=k
      enddo
! maritime size distribution after hoppel et al. 1994, jgr. 14,443
!      wr(3,1)=0.02
!      wr(3,2)=0.05
!      wr(3,3)=0.15
!      wn(3,1)=110.
!      wn(3,2)=72.
!      wn(3,3)=7.
!      ws(3,1)=0.14
!      ws(3,2)=0.16
!      ws(3,3)=0.18
!      x0=sqrt(2.*pi)
!      do k=1,3
!         wn(3,k)=wn(3,k)/(x0*ws(3,k))
!         ws(3,k)=1./(2.*ws(3,k)**2)
!      enddo
      do k=1,n
         do ia=1,nka
         do jt=2,nkt
            ff(jt,ia,k)=0.
         enddo
         enddo
      enddo
      do ia=1,nka
         ff(1,ia,1)=0.
      enddo
      do k=1,n
        fsum(k)=0.
        nar(k)=iaertyp
        x0=1.0
        if (iaertyp.lt.3.and.k.gt.nf) x0=0.2
        do ia=1,nka
          ff(1,ia,k)=dfdlogr(rn(ia),nar(k))*dlgenw/3.*x0
          if (k.gt.kinv) ff(1,ia,k)=dfdlogr2(rn(ia),nar(k))*dlgenw/3.*x0
          fsum(k)=fsum(k)+ff(1,ia,k)
        enddo
!        write (199,*)"k,fsum",k,fsum(k)
!        fnorm(k)=fsum(k) ! jjb variable unreferenced
      enddo
! read initial aerosol distribution from previous run
!#      fname='ae .out'
!#      fname(3:3)=fogtype
!#      open (66, file=fname,status='old',form='unformatted')
!#      read (66) aer
!#      close (66)
!#      do k=1,nf
!#         fsum(k)=0.
!#         do ia=1,nka
!#            ff(1,ia,k)=aer(k,ia)
!#            fsum(k)=fsum(k)+ff(1,ia,k)
!#         enddo
!#      enddo
! normalization of aerosol size distribution
!      do k=2,nf
!         x0=fnorm(k)/fsum(k)
!         fsum(k)=0.
!         do ia=1,nka
!            ff(1,ia,k)=ff(1,ia,k)*x0
!            fsum(k)=fsum(k)+ff(1,ia,k)
!         enddo
!      enddo
! parameters a0m, b0m of koehler curve of subroutine subkon: 
! sr(ia,jt)=exp(a0m/(r(ia,jt)*t)-b0m(ia)*en(ia)/ew(jt)))
! a0m see p. 141 pruppacher and klett a0m=2 sigma/(r1*t*rhow*a)
! 152200= 2 sigma*10**6 with sigma is surface tension = 76.1*10**-3
! see pruppacher and klett p. 104
      a0m=152200./(r1*rhow)
! aerosol types: 1=urban 2=rural 3=ocean 4=tropospheric
      k0=nar(2)
      do 1030 ia=1,nka
      if (k0-2) 2000,2010,2020
! b0m=fcs*xnue*xmol2/xmol3: fcs(ia) fraction of soluble species
! xnue number of ions; xmol2 (xmol3) mol masses of water (aerosol)
! NH4NO3 mole mass 80; (NH4)2SO4 mole mass 132
! soluble part of urban aerosol: 2 mole NH4NO3 and 1 mole (NH4)2SO4
 2000 continue
      fcs(ia)=.4-rn(ia)*(.4-.1)
      if (rn(ia).gt.1.) fcs(ia)=.1
      xmol3(ia)=(132.+80.*2.)/3.
      xnue=(3.+2.*2.)/3.
      go to 1030
! soluble part of rural aerosol: pure (NH4)2SO4
! 2010 fcs(ia)=.5-rn(ia)*(.5-.1)
 2010 continue
      fcs(ia)=.9-rn(ia)*(.9-.5)
      if (rn(ia).gt.1.) fcs(ia)=.5
!      fcs(ia)=0.1
      xnue=3.
      xmol3(ia)=132.
      go to 1030
!c soluble part of ocean aerosol: small pure (NH4)2SO4; large pure NaCl
! soluble part of ocean aerosol: pure (NH4)2SO4; 
 2020 continue
      fcs(ia)=1.
!      xnue=3.                   
!      xmol3(ia)=132.
! 32% (NH4)2SO4, 64% NH4HSO4, 4% NH4NO3
      xnue=0.32*3.+0.64*2+0.04*2                   
      xmol3(ia)=0.32*132.+0.64*115+0.04*80
! large are NaCl
      if (rn(ia).lt.0.5) go to 1030
      xnue=2.     !no change in microphysics due to halogen chemistry
      xmol3(ia)=58.4
 1030 b0m(ia)=fcs(ia)*xnue*xmol2/xmol3(ia)
! initial temperature
! in radiation code: background aerosol = rural aerosol
      t(1)=tw
      do k=2,kinv
         t(k)=t(1)-0.0098*eta(k)
!         t(k)=t(1)-0.0096*eta(k)
         if (nar(k).eq.4) nar(k)=2
      enddo
      x0=t(kinv)+3.
!      x0=t(kinv) ! no strong inversion
      do k=kinv+1,n !kinv+10
         t(k)=x0-0.006*(eta(k)-eta(kinv))
!;         t(k)=x0-0.008*(eta(k)-eta(kinv))
!;         t(k)=x0-0.0098*(eta(k)-eta(kinv))
         if (nar(k).eq.4) nar(k)=2
      enddo
!      do k=kinv+11,n
!         t(k)=x0-0.006*(eta(k)-eta(kinv))
!         if (nar(k).eq.4) nar(k)=2
!      enddo

!      call adjust_f !only if Kim aerosol is used!!

! large scale hydrostatic pressure
      poben=101400.
      cc=g/(2.*r0)
      do k=1,n
         punten=poben
         dd=detw(k)*cc/t(k)
         poben=punten*(1.-dd)/(1.+dd)
         p(k)=0.5*(poben+punten)
      enddo
! initial profiles of humidity and wind field
! xm1: = specific humidity in kg/kg
! xm2: = liquid water content in kg/m**3
      do k=1,n
         talt(k)=t(k)
         thet(k)=(p(1)/p(k))**0.286
         theti(k)=1./thet(k)
         theta(k)=t(k)*thet(k)
! r0/r1=287.05/461.51=0.62198; 1-r0/r1=0.37802;
         xm21s=0.62198*p21(t(k))/(p(k)-0.37802*p21(t(k)))
!         xm1(k)=dmin1(8.5d-03,xm21s)
!#         xm1(k)=dmin1(8.5d-03,0.95*xm21s)
!         xm1(k)=dmin1(8.5d-03,0.7*xm21s)
         xm1(k)=dmin1(8.5d-03,0.8*xm21s)
!         xm1(k)=dmin1(8.5d-03,0.4*xm21s) no cloud
!         if (eta(k).gt.zinv)  xm1(k)=7.6d-03-4.d-03/1000.*(eta(k)-zinv)
!#         if (eta(k).gt.zinv)  xm1(k)=4.d-03
         if (eta(k).gt.zinv)  xm1(k)=dmin1(4.d-03,0.4*xm21s) 
!;         if (eta(k).gt.zinv)  xm1(k)=dmin1(4.d-03,0.35*xm21s) 
!         if (eta(k).gt.zinv)  xm1(k)=dmin1(4.d-03,0.3*xm21s) 
!         if (eta(k).gt.zinv)  xm1(k)=dmin1(4.d-03,0.3*xm21s) no cloud
         feu(k)=xm1(k)*p(k)/((0.62198+0.37802*xm1(k))*p21(t(k)))
! r1/r0-1=0.61
         rho(k)=p(k)/(r0*(t(k)*(1.+0.61*xm1(k))))
         thetl(k)=theta(k)*(1.+0.61*xm1(k))
         xm1a(k)=xm1(k)
         dfddt(k)=0.
         xm1a(k)=xm1(k)
         xm2(k)=0.
         xm2a(k)=0.
         u(k)=ug
         v(k)=vg  
         w(k)=eta(k)/1000.*0.5*(wmin+wmax)
         tke(k)=1.d-05
         if (eta(k).lt.zinv) tke(k)=0.05
      enddo
      u(1)=0.
      v(1)=0.
      u(2)=0.25*ug
      v(2)=0.25*vg
      u(3)=0.75*ug
      v(3)=0.75*vg
      do k=n,1,-1
         w(k)=w(k)-w(1)
      enddo
      vbt=sqrt(u(2)*u(2)+v(2)*v(2))
      zp=deta(1)+z0
      zpdz0=dlog(zp/z0)
      zpdl=g*(theta(2)-t(1))*zp/(theta(2)*vbt)
      call claf (zpdl,zpdz0,cu,ctq)
      ustern=dmax1(0.01d0,vbt/cu)
      gclu=cu
      gclt=ctq
      ajs=0.
      ds1=0.
      ds2=0.
      trdep=0.
      tau=0.
      reif=0.
! temperature and volumetric moisture within soil
      x0=0.5*ebs
!      if (iaertyp.eq.1) x0=x0*.9
      do k=1,nb
         tb(k)=285.0
         eb(k)=x0
         if (zb(k).lt..1) tb(k)=(t(1)*(.1-zb(k))+285.*zb(k))/.1
      enddo
! initial output for plotting
      fname='pi .out'
      fname(3:3)=fogtype
      open (97, file=fname,status='unknown')
      write (97,6000) (eta(k),etw(k),rho(k),p(k),w(k),k=1,n)
 6000 format (5e16.8)
      close (97)
      write (19) zb
      close (19)
      do jt=1,nkt
!         jtp=min0(jt+1,nkt) ! jjb variable unreferenced
!         de0=dew(jt)  ! jjb variable unreferenced
!         dep=dew(jtp) ! jjb variable unreferenced
!         de0p=de0+dep ! jjb variable unreferenced
         do ia=1,nka
            rk=rw(jt,ia)
            sr(ia,jt)=dmax1(.1d0,dexp(a0m/(rk*t(2))
     &                -b0m(ia)*en(ia)/ew(jt)))
         enddo
      enddo
      fname='fi .out'
      fname(3:3)=fogtype
      open (44, file=fname,status='unknown')
      write (44,6010) rn,en,rq,e,sr
      close (44)
 6010 format (5e16.8)

      end subroutine initm

!
!-------------------------------------------------------------
!

      subroutine grid

      USE constants, ONLY :
! Imported Parameters:
     &     pi,
     &     rhow                  ! Water density [kg/m**3]

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nb,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)
! numerical grid for the atmosphere, aerosols and water droplets

! jjb test new grid
      double precision, parameter :: rho_aerosol = 2000  ! [kg/m**3]
      double precision, parameter :: kg_2_mg     = 1.e6  ! [mg/kg]
      double precision, parameter :: vol_fact    = 4./3.*pi ! [-]
      double precision, parameter :: r_min     =  0.005 ! [um]
      double precision, parameter :: r_max_dry = 15.    ! [um]
      double precision, parameter :: r_max_wet = 60.    ! [um]


      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb51/ dlgew,dlgenw,dlne
      common /blck06/ kw(nka),ka
! density and molecular weight of pure water
! minimum (rnw0) maximum (rnw1) radius of dry aerosol   [um]
! minimum (rw0) maximum (rw1) radius of total particle  [um]
      data rnw0,rnw1,rw0,rw1 /0.005,15.,0.005,150./
! minimum-maximum values of soil and atmospheric grid
      data dzbw0,zbw1,detamin,etaw1 /0.001,1.,10.,2000./

! vertical grid
! equidistant grid between earth's surface and eta(nf)
      etw(1)=0.
      do 1000 k=2,nf
 1000 etw(k)=float(k-1)*detamin
! nonequidistant grid above eta(nf)
      nf1=nf+1
      j=0
      x0=detamin
 3000 x0=x0+detamin
      j=j+1
      x3=detamin/x0+1.
      etw(nf1)=x0
      do k=nf+2,n
         etw(k)=etw(k-1)*x3
      enddo
      if (j.gt.10000) stop 4000
      if (etw(n)-etw(nf1).gt.etaw1-etw(nf)) go to 3000
      x0=nf*detamin-etw(nf1)
      do k=nf1,n
         etw(k)=etw(k)+x0
      enddo
      detw(1)=detamin
      eta(1)=0.
      do k=2,n
         detw(k)=etw(k)-etw(k-1)
         eta(k)=0.5*(etw(k)+etw(k-1))
         deta(k-1)=eta(k)-eta(k-1)
      enddo
      deta(n)=(1.+x3)*0.5*etw(n)-eta(n)

! grid within the soil
      zbw0=0.
 3010 zbw0=zbw0+0.0001
      dlgzbw=dlog10(zbw1/zbw0)/float(nb)
      x3=10.**dlgzbw
      zbw=zbw0*x3
      if (zbw-zbw0.lt.dzbw0) go to 3010
      zb(1)=zbw
      dzbw(1)=zbw-zbw0
      do k=2,nb
         zbw0=zbw
         zbw=zbw0*x3
         zb(k)=0.5*(zbw+zbw0)
         dzbw(k)=zbw-zbw0
         dzb(k-1)=zb(k)-zb(k-1)
      enddo
      dzb(nb)=(1.+x3)*0.5*zbw-zb(nb)
      x0=zb(1)
      do k=1,nb
         zb(k)=zb(k)-x0
      enddo

! aerosol grid
      x0=1./3.
      x1=4.*x0*pi*rhow
      rho3=2000.          ! jjb: dry aerosol density?
      x2=4.*x0*pi*rho3

! aerosol mass en(i) in mg
      enwmin=x2*rnw0**3*1.d-12
      enwmax=x2*rnw1**3*1.d-12
! begin prescribed lowest and largest aerosol class
      dlgenw=dlog10(enwmax/enwmin)/float(nka)
      x3=10.**dlgenw
      enw(1)=enwmin*x3
      en(1)=0.5*(enw(1)+enwmin)
      rn(1)=(en(1)/x2)**x0*1.d+04
      do ia=2,nka
         enw(ia)=enw(ia-1)*x3
         en(ia)=0.5*(enw(ia)+enw(ia-1))
         rn(ia)=(en(ia)/x2)**x0*1.d+04
      enddo
      print*,enwmax,enw(nka)
! ax: growth factor for consecutive masses
! old water grid
      ewmin=x1*rw0**3*1.d-12
      ewmax=x1*rw1**3*1.d-12
      dlgew=dlog10(ewmax/ewmin)/float(nkt) ! jjb nkt, not nkt-1
      ax=10.**dlgew
! new water grid 
! prescribed doubling of water mass after scal bins 
!      scal=1.36
! ln(10)=2.3025851  ln(x)=log(10)*log(x)
!      dlgew=dlog10(scal)
      dlne=2.3025851*dlgew
!      ax=2.d0**(1.0/scal)
      ew(1)=ewmin*ax
      e(1)=0.5*(ew(1)+ewmin)
      dew(1)=ew(1)-ewmin
      do jt=2,nkt
         ew(jt)=ew(jt-1)*ax
         e(jt)=0.5*(ew(jt)+ew(jt-1))
         dew(jt)=ew(jt)-ew(jt-1)
      enddo
      do ia=1,nka
      do jt=1,nkt
         rq(jt,ia)=(e(jt)*1.d-06/x1+(rn(ia)*1.d-06)**3)**x0*1.d+06
!         rq(jt,ia)=dmax1(rq(jt,ia),rn(ia))
         rw(jt,ia)=(ew(jt)*1.d-06/x1+(rn(ia)*1.d-06)**3)**x0*1.d+06
      enddo
      enddo
! partitioning of particle spectrum into aerosols, first and second
! liquid species region for chemical reactions
!      kg=-1
!      kd=-1
!      do jt=1,nkt
!         if (rq(jt,1).gt.3.d0.and.kg.lt.0) kg=jt-1
!         if (rq(jt,1).gt.15.d0.and.kd.lt.0) kd=jt-1
!      enddo
! aerosol + cloud chemistry: partition into small and large aerosol  : ka
!                       and small and large water content (aer-drop) : kw(nka)
      ka=-1
      do ia=1,nka
         if (rn(ia).gt..5d0.and.ka.lt.0) ka=ia-1
         kw(ia)=-1
      enddo
      if (ka.lt.0) ka=nka
! xfac: dilution as volume ratio: V_dry*x = V_water (r_dry*(x)^(1./3.)=r_water)
      xfac=10. ! volume ratio is 1000
      do ia=1,nka
         do jt=1,nkt
            if ((e(jt)*1.e-6/x1)**(1./3.)*1.e6.gt.xfac*rn(ia).and.
     &           kw(ia).lt.0) kw(ia)=jt-1
         enddo
         if (kw(ia).lt.0) kw(ia)=nkt
      enddo
! end aerosol + cloud chemistry



      end subroutine grid

!
!-------------------------------------------------------------
!

      subroutine startm (fogtype)

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nrlay,
     &     nb,
     &     nka,
     &     nkt,
     &     mb

      implicit double precision (a-h,o-z)
! vertical profiles of meteorological data if the program is restarted

      common /cb11/ totrad (mb,nrlay)
      double precision totrad

      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb53a/ thet(n),theti(n)
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /cb55/ dtrad0(n),dtrad1(n),sk0,sl0,sk1,sl1,time0,time2
      common /cb63/ fcs(nka),xmol3(nka)
      character *10 fname
      character *1 fogtype
      fname='rstm .dat'
      fname(5:5)=fogtype
      open (15,file=fname,status='unknown',form='unformatted')
! double precision arrays
      read (15) 
     &     atkm,atkh,b0m,dfddt,dtrad,dtrad0,dtrad1,eb,ff,fcs,feu,fsum,
     &     gh,p,rho,t,talt,tb,tke,tkep,theta,totrad,u,v,w,xl,xm1,xm1a,
     &     xm2,xmol3,
! double precision single vars
     &     a0m,alat,declin,ds1,ds2,reif,sk,sk0,sk1,sl,sl0,sl1,tau,time0,
     &     time2,trdep,
! integer arrays
     &     nar,
! integer single vars
     &     it,lcl,lct,lday,lmin,lst

! set/diagnose some parameters
      tw = t(1)
      do k=1,n
         thet(k)=theta(k)/t(k)
         theti(k)=1./thet(k)
         thetl(k)=theta(k)*(1.+0.61*xm1(k))
      enddo
      vbt=sqrt(u(2)*u(2)+v(2)*v(2))
      zp=deta(1)+z0
      zpdz0=dlog(zp/z0)
      zpdl=g*(theta(2)-t(1))*zp/(theta(2)*vbt)
      call claf (zpdl,zpdz0,cu,ctq)
      ustern=dmax1(0.01d0,vbt/cu)
      gclu=cu
      gclt=ctq
      close (15)

      print *,"restart file for meteo read, filename: ",fname
      print *,lday,lst,lmin

      end subroutine startm

!
!-------------------------------------------------------------
!

      subroutine startc (fogtype)

      USE gas_common, ONLY :
     &     s1,
     &     es1,
     &     s3

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc,
     &     nlev,
     &     nrxn,
     &     nphrxn

      implicit double precision (a-h,o-z)

! profiles of chemical data if the program is restarted

      common /band_rat/ photol_j(nphrxn,n)
      common /blck01/ am3(n),cm3(n)
      common /blck11/ rc(nkc,n)
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck13/ conv2(nkc,n) ! conversion factor = 1/(1000*cw)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /blck78/ sa1(nka,j2),sac1(nka,j2)

      common /budg/ bg(2,nrxn,nlev),il(nlev)
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /kinv_i/ kinv
      common /kpp_crys/ xcryssulf,xcrysss,xdelisulf,xdeliss

      common /kpp_l1/ cloudt(nkc,n)
      logical cloudt

      common /kpp_mol/ xgamma(nf,j6,nkc)
      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
      double precision is4(n,3)
      character *10 fname
      character *1 fogtype
      fname='rstc .dat'
      fname(5:5)=fogtype
      open (16,file=fname,status='unknown',form='unformatted')
! alpha, henry, vmean, xkmt are not read in because they depend on the
! chemical mechanism used (aer, tot) and are calculated every time step 
! anyways

! double precision arrays
      read (16) am3,cm,cm3,conv2,cw,es1,photol_j,rc,s1,s3,sa1,
     &     sac1,sl1,sion1,vd,vdm,vt,xgamma,
! double precision, single values
     &     xcryssulf,xcrysss,xdelisulf,xdeliss,
! logicals
     &     cloudt,
! integers
     &     il,kinv,lday,lmin,lst
      close (16)

! initial output for plotting - only needed for binary output
      do k=1,n
            is4(k,1)=am3(k) ! stay consistent with plot routine array size!
            is4(k,2)=0.
            is4(k,3)=0.
!            write (542,12) k,am3(k,1),am3(k,2)
!            write (543,12) k,cm3(k,1),cm3(k,2)
      enddo
! 12   format (i3,2d16.8)
      write (61) is4
      close (61)
      write (64) is4
      close (64) 
      
      print *,"restart file for chemistry read, filename: ",fname
      print *,lday,lst,lmin

      print *,'deliquescense',xdelisulf,xdeliss
      print *,'crystallization',xcryssulf,xcrysss
! get alpha's
      call st_coeff_a
      call st_coeff_t
      
      end subroutine startc

!
!-------------------------------------------------------------
!

      subroutine vgleich

      USE constants, ONLY :
! Imported Parameters:
     &     pi

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)
! equilibrium values for radius rg(i) and water mass eg(i)
! of humidified aerosol particles at given relative humidity
! new distribution of the particles on their equilibrium positions

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      dimension rg(nka),eg(nka)

      do k=2,n
         feu(k)=dmin1(feu(k),0.99999d0)
! equilibrium radius
         a0=a0m/t(k)
         do ia=1,nka
            b0=b0m(ia)*2.
! b0=b0m*rho3/rhow; rho3=2000; rhow=1000
            rg(ia)=rgl(rn(ia),a0,b0,feu(k))
            eg(ia)=4.d-09*pi/3.*(rg(ia)**3-rn(ia)**3)
         enddo
! new distribution
         do 1020 ia=1,nka
         do jt=1,nkt
            if (eg(ia).le.ew(jt)) then
               if (jt.eq.1) goto 1020
               ff(jt,ia,k)=ff(1,ia,k)
               ff(1,ia,k)=0.
               goto 1020
            endif
         enddo
 1020    continue
! total liquid water
         xm2(k)=0.
         do ia=1,nka
         do jt=1,nkt
            xm2(k)=xm2(k)+ff(jt,ia,k)*e(jt)
         enddo
         enddo
      enddo

      end subroutine vgleich

!
!-------------------------------------------------------------
!

      function rgl (r_dry,a,b,feu)
      implicit double precision (a-h,o-z)
! equilibrium radius of aerosol particle at given relative humidity
! r_dry dry particle radius; a scaling radius for surface effect
      f(x)=(x**3-1.)*(x*dlogf-alpha)+b*x
      fstr(x)=(4.*x**3-1.)*dlogf-3.*x*x*alpha+b
      if (feu.ge.1.) then
         write (6,6000)
 6000    format (10x,'initial value of relative humidity exceeding one')
         rgl=r_dry
         return
      endif
      dlogf=dlog(feu)
      alpha=a/r_dry
! initial value
      xalt=dexp(feu)
! Newton iteration
      do ij=1,100
         xneu=xalt-f(xalt)/fstr(xalt)
         if (abs(xneu-xalt).lt.1.d-7*xalt) goto 2000
         xalt=xneu
      enddo
 2000 continue
      rgl=r_dry*xneu

      end function rgl

!
!-------------------------------------------------------------
!

      subroutine sedp (dt)
! gravitational settling of particles with terminal velocity w in m/s
! for further details on the determination of w see function vterm

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nb,
     &     nka,
     &     nkt,
     &     nkc

      implicit double precision (a-h,o-z)

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb58/ c(nf),psi(nf)
      common /blck06/ kw(nka),ka
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)

      ajs=0.
      c(nf)=0.
      x3=-deta(2)

      do ia=1,nka
         do jt=1,nkt
!            x4=rq(jt,ia)
!            ww=-1.25d-4*x4*x4*(1.+8.6d-02/x4)
            ww=-1.*vterm(rq(jt,ia)*1.d-6,t(nf),p(nf)) !"first guess" for determination of ww
            dt0=dt     
            xsum=0.
            do k=2,nf
               psi(k)=ff(jt,ia,k)
               xsum=xsum+psi(k)
            enddo
            if (xsum.gt.1.d-06) then
               x0=0.
! 3000          dtmax=dmin1(dt0,x3/(ww+w(nf)))
! see also SR difp: subsidence treated now consistently (i.e. like for other 
! tracers): w df/dz instead of d(wf)/dz
 3000          dtmax=dmin1(dt0,x3/(ww))  
               do k=2,nf
!                  c(k)=dtmax/deta(k)*(ww+w(k))
                  c(k)=dtmax/deta(k)*(-1.*vterm(rq(jt,ia)*1.d-6,t(k)
     &                 ,p(k)))
!     &                 ,p(k))+w(k))
               enddo
! particle dry deposition velocity in lowest model layer:
               c(2)=dmin1(c(2),dtmax/deta(k)*vd(jt,ia)*(-1.))
               c(1)=c(2)
               dt0=dt0-dtmax
               x1=psi(2)
               psi(1)=x1
               if (rq(jt,ia).lt.1.) then
                  call advsed0
               else
                  call advsed1
               endif
               x0=x0+psi(1)-x1
               if (dt0.gt.0.1) go to 3000
               do k=2,nf-1
                  ff(jt,ia,k)=psi(k)
               enddo
               ff(jt,ia,nf)=ff(jt,ia,nf-1)
            endif
!         enddo
! droplet sedimentation
            x2=x0*e(jt)*detw(2)
            ajs=ajs+x2/dt
            trdep=trdep+x2
c Droplet sedimentation has been evaluated at ground
c trdep :  cumulative deposition [kg liquid water/m^2]
c          being equivalent to [mm] precipitation
c ajs   :  sedimentation rate [kg liquid water/m^2/sec]
c          being equivalent to precip. rate of [mm/sec]
c update total liquid water [kg/m^3]

!            if (jt.ge.kgp) then
!               if (jt.lt.kdp) then
!                  ds1=ds1+x2
!               else
!                  ds2=ds2+x2
!               endif
!            endif
            if (jt.gt.kw(ia)) then
               ds1=ds1+x2
            else
               ds2=ds2+x2
            endif
         enddo
      enddo

      end subroutine sedp

!
!-------------------------------------------------------------
!

      subroutine sedc (dt)
! dry deposition and emission of gaseous species

! jjb work done = implicit none, missing declarations, little cleaning, modules including constants

      USE constants, ONLY :
! Imported Parameters:
     &     Avogadro

      USE gas_common, ONLY:
! Imported Parameters:
     &     j1,
! Imported Array Variables with intent (in):
     &     es1,
     &     ind_gas_rev,
! Imported Array Variables with intent (inout):
     &     s1,
     &     vg

      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) ::  dt

! Local scalars:
      integer j
      double precision s12old
      double precision x4, w

! Common blocks:
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

!- End of header ---------------------------------------------------------------


! constant dry deposition velocities in m/sec, not for all species defined
! old data set:
!      data vg/0.10e-2,0.10e-2,0.10e-1,0.1e-1,0.8e-02,
!c     &      0.01e-2,0.50e-2,7*0.10e-2,2*0.0,0.30e-2,
!     &      0.01e-2,0.50e-3,7*0.10e-2,2*0.0,0.30e-2,
!     &      3*0.2e-2,0.50e-2,2*0.20e-02,2*0.10e-2,
!c     &      3*0.20e-2,0.10e-1,0.,3*0.01e-2,3*0.10e-1,
!c     &      4*0.0/
!     &      3*0.20e-2,0.10e-1,0.2e-1,3*0.01e-2,3*0.10e-1,
!     &      0.1e-2,0.2e-3,3*0.0,0.2e-1,0.2e-2,4*0.,
!     &      6*0.1e-1,2*0.0,0.1e-1,9*0.0/

! Laurens data set:
!x      data vg/0.0,0.10e-4,0.050e-1,0.27e-2,0.5e-02,
!x     &      1.0e-2,0.40e-3,7*0.0,0.5e-2,0.0,0.30e-2,
!x     &      0.2e-2,0.4e-2,0.2e-2,0.0,2*0.2e-02,2*0.10e-2,
!x     &      3*0.20e-2,0.10e-1,0.2e-1,3*0.01e-2,0.10e-1,0.0,0.0,
!x     &      0.0,0.2e-2,3*0.0,0.2e-1,0.2e-2,4*0.,
!x     &      5*0.1e-1,4*0.0,3*0.0,1.e-2,2*0.0,1.e-2,2*0.0/

!      vg(4)=0.27e-2 ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or 
!      vg(34)=vg(30) ! N2O5=HCl
!      vg(37)=0.     ! DMS           emission is net flux 
!      vg(38)=vg(30) ! HOCl = HCl
!      vg(43)=vg(30) ! HOBr = HCl
!      vg(50)=vg(49) ! I2O2=HOI
!      vg(51)=vg(49) ! INO2=HOI
!      vg(56)=0.     ! CH3I          emission is net flux
!      vg(57)=0.     ! CH2I2         emission is net flux
!      vg(58)=0.     ! CH2ClI        emission is net flux
!      vg(59)=0.     ! C3H7I         emission is net flux
!      vg(63)=vg(30) ! CH3SO3H = HCl
!      vg(71)=0.     ! CH2BrI         emission is net flux
!      vg(72)=0.     ! CHBr2I         emission is net flux
!      vg(73)=0.     ! C2H5I          emission is net flux

      if(ind_gas_rev(4) /= 0)
     & vg(ind_gas_rev(4))=0.27e-2              ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or 
      if(ind_gas_rev(34) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(34))=vg(ind_gas_rev(30)) ! N2O5=HCl
      if(ind_gas_rev(37) /= 0)
     &  vg(ind_gas_rev(37))=0.                  ! DMS           emission is net flux 
      if(ind_gas_rev(38) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(38))=vg(ind_gas_rev(30)) ! HOCl = HCl
      if(ind_gas_rev(43) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(43))=vg(ind_gas_rev(30)) ! HOBr = HCl
      if(ind_gas_rev(50) /= 0 .and. ind_gas_rev(49) /= 0)
     &  vg(ind_gas_rev(50))=vg(ind_gas_rev(49)) ! I2O2=HOI
      if(ind_gas_rev(51) /= 0 .and. ind_gas_rev(49) /= 0)
     &  vg(ind_gas_rev(51))=vg(ind_gas_rev(49)) ! INO2=HOI
      if(ind_gas_rev(56) /= 0)
     &  vg(ind_gas_rev(56))=0.                  ! CH3I          emission is net flux
      if(ind_gas_rev(57) /= 0)
     &  vg(ind_gas_rev(57))=0.                  ! CH2I2         emission is net flux
      if(ind_gas_rev(58) /= 0)
     &  vg(ind_gas_rev(58))=0.                  ! CH2ClI        emission is net flux
      if(ind_gas_rev(59) /= 0)
     &  vg(ind_gas_rev(59))=0.                  ! C3H7I         emission is net flux
      if(ind_gas_rev(63) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(63))=vg(ind_gas_rev(30)) ! CH3SO3H = HCl
      if(ind_gas_rev(71) /= 0)
     &  vg(ind_gas_rev(71))=0.                  ! CH2BrI         emission is net flux
      if(ind_gas_rev(72) /= 0)
     &  vg(ind_gas_rev(72))=0.                  ! CHBr2I         emission is net flux
      if(ind_gas_rev(73) /= 0)
     &  vg(ind_gas_rev(73))=0.                  ! C2H5I          emission is net flux

      if (lst/4*4.eq.lst.and.lmin.eq.1) then
         print *,lday,lst,lmin
         print *,' dry deposition velocities'
         do j=1,j1
            print *,j,vg(j)
         enddo
      endif


!      x3=deta(2)
! emission rates: 80% of total emission during day 20% during night
! 1.6e-2=1.6/100: factor 100 due to detw(2) with units in cm
!      x0=float(lmin)/60.
! emission rates of chemical species variable in time
!      x4=1.6e-02*dt
!      if (lst.ge.19.or.lst.lt.6) x4=.4e-02*dt
!      if (lst.eq.6) x4=(1.6e-02*x0+0.4e-02*(1.-x0))*dt
!      if (lst.eq.18) x4=(0.4e-02*x0+1.6e-02*(1.-x0))*dt
      x4=1.
!      dt0=dt
      do j=1,j1
         w=vg(j)
         if (w.ge.1.e-05) then
!            psi2=s1(j,2)*detw(2)
!            psi3=s1(j,3)*detw(3)
!            dtmax=dmin1(dt0,x3/w)
!            cl=dtmax*w/x3
!            a0=(25.*psi2-psi3)/24.
!            a1=(psi3-psi2)/16.
!            a2=(psi3-psi2)/48.
!            x1=1.-2.*cl
!            x2=x1*x1
!            xxx2=a0*cl-a1*(1.0-x2)+a2*(1.0-x1*x2)
!            fm=dmax1(0.d0,xxx2)
!            flux1=fm*psi2/dmax1(fm+1.e-15,psi2)
!            s1(j,2)=(psi2-flux1)/detw(2)
!c total dry deposition of gas phase phase species in mol/m^2        : depfac=1.
!            depfac=1.
!c            s1(j,1)=s1(j,1)+flux1*100.
!            s1(j,1)=s1(j,1)+flux1*depfac
            s12old=s1(j,2)
            s1(j,2)=s1(j,2)*exp(-dt/deta(2)*vg(j))
            s1(j,1)=s1(j,1)+(s12old-s1(j,2))*deta(2)
         endif
! es1: emission rates in molec./cm**2/s, s1 in mol/m**3            
         s1(j,2)=s1(j,2)+es1(j)*x4*dt*1.e+4/(detw(2)*Avogadro)
      enddo

      end subroutine sedc

!
!----------------------------------------------------------------
!

      subroutine sedl (dt)
! new aqueous phase concentrations due to 
! gravitational settling of droplets

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc

      implicit double precision (a-h,o-z)
!      double precision dt,f,fsum,c,psi,detw,deta,eta,etw

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb58/ c(nf),psi(nf)
      common /blck11/ rc(nkc,n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
      common /liq_pl/ nkc_l
      dimension cc(nf)
      c(nf)=0.
! changes have to be made BOTH here and below for ions
! rc in m 
      xfac=1.e6
      do kc=1,nkc_l
         do k=2,nf
            xxx=0.01 ! jjb apparently, in um
            x4=dmax1(xxx,xfac*rc(kc,k)) ! jjb here as well
! subsidence see SR difl
!            cc(k)=(-1.25e-4*x4*x4*(1.+8.6e-02/x4))/deta(k) 
            cc(k)=(-1.*vterm(x4*1.d-6,t(k),p(k)))/deta(k) ! jjb: in vterm, radius in m
! mass weighted terminal velocity
            cc(k)=dmin1(cc(k),-1.*vt(kc,k)/deta(k))
         enddo
! particle dry deposition velocity in lowest model layer:
         cc(2)=dmin1(cc(2),-1./deta(k)*vdm(kc))
         do l=1,j2
            do k=2,nf
               psi(k)=sl1(l,kc,k)
            enddo
            dt0=dt
            x0=0.
            xxxt=-.999/cc(2)
 3000       continue
            dtmax=dmin1(dt0,xxxt)
            dt0=dt0-dtmax
            do k=2,nf
               c(k)=cc(k)*dtmax
            enddo
            c(1)=c(2)
            x1=psi(2)
            psi(1)=x1
            call advsed1
            x0=x0+psi(1)-x1
            if (dt0.gt.0.1) go to 3000
            do k=2,nf-1
               sl1(l,kc,k)=psi(k)
            enddo
! wet deposition to the ground in mole/m**2
            sl1(l,kc,1)=sl1(l,kc,1)+x0*deta(2)
         enddo
      enddo
! dito for ions
      c(nf)=0.
      do kc=1,nkc_l
         do k=2,nf
            xxx=0.01
            x4=dmax1(xxx,xfac*rc(kc,k))
! subsidence see SR difl
!            cc(k)=(-1.25e-4*x4*x4*(1.+8.6e-02/x4))/deta(k)
            cc(k)=(-1.*vterm(x4*1.d-6,t(k),p(k)))/deta(k)
! mass weighted terminal velocity
            cc(k)=dmin1(cc(k),-1.*vt(kc,k)/deta(k))
         enddo
! particle dry deposition velocity in lowest model layer:
         cc(2)=dmin1(cc(2),-1./deta(k)*vdm(kc))
         do l=1,j6
            do k=2,nf
               psi(k)=sion1(l,kc,k)
            enddo
            dt0=dt
            x0=0.
            xxxt=-.999/cc(2)
 3010       continue
            dtmax=dmin1(dt0,xxxt)
            dt0=dt0-dtmax
            do k=2,nf
               c(k)=cc(k)*dtmax
            enddo
            c(1)=c(2)
            x1=psi(2)
            psi(1)=x1
            call advsed1
            x0=x0+psi(1)-x1
            if (dt0.gt.0.1) go to 3010
            do k=2,nf-1
               sion1(l,kc,k)=psi(k)
            enddo
! wet deposition to the ground in mole/m**2
            sion1(l,kc,1)=sion1(l,kc,1)+x0*deta(2)
         enddo
      enddo !ions

      end subroutine sedl

!
!------------------------------------------------------------
!

      function vterm(a,t,p)

! terminal velocity for droplets
! after Stokes with Cunningham correction for small droplets (regime 1)
! and after Beard for large droplets (r > 10 micron, regime 2)
! all formulas after Pruppacher and Klett Chapter 10.

      USE constants, ONLY :
! Imported Parameters:
     &     g,
     &     r0,                   ! Specific gas constant of dry air, in J/(kg.K)
     &     rhow                  ! Water density [kg/m**3]

      implicit none

      double precision :: vterm

      double precision, intent(in) :: a ! radius       in [m]
      double precision, intent(in) :: t ! temperature  in [K]
      double precision, intent(in) :: p ! pressure     in [Pa]

      ! Polynomial coefficients for Beard approximation
      ! Pruppacher & Klett, p. 417, equation (10-145)
      double precision, parameter :: b0=-.318657d+1
      double precision, parameter :: b1= .992696d+0 
      double precision, parameter :: b2=-.153193d-2
      double precision, parameter :: b3=-.987059d-3
      double precision, parameter :: b4=-.578878d-3
      double precision, parameter :: b5=+.855176d-4
      double precision, parameter :: b6=-.327815d-5

      double precision, parameter :: c1 = 2.d0 * g / 9.d0       ! constant factor in equation (10-138)
      double precision, parameter :: c2 = 1.26d0                ! constant in equation (10-139)
      double precision, parameter :: P0 = 101325                ! Standard pressure    [Pa]
      double precision, parameter :: T0 = 293.15                ! Standard temperature [K]
      double precision, parameter :: lambda0 = 6.6d-8           ! mean free path STP   [m]
      double precision, parameter :: c3 = c2 * lambda0 * P0/T0  ! constant factor in equation (10-139) and (10-140)
      double precision, parameter :: c4 = 32.d0 * g / 3.d0      ! constant factor in equation (10-142)

      double precision :: best  ! Davies or Best number, see equation (10-142) [-]
      double precision :: x     ! ln(best)                                     [-]
      double precision :: rho_a ! air density                                  [kg/m3]
      double precision :: eta   ! dynamic viscosity                            [kg/(m.s)]
      double precision :: y     ! Reynolds number = exp(y) in Regime 2         [-]

      rho_a=p/(r0*t)
      eta=3.7957d-06+4.9d-08*t

      if (a.le.1.d-5) then
         ! Regime 1: see P & K pp. 415-417, equations (10-138) to (10-140)
         vterm=c1*a*a*(rhow-rho_a)/eta*(1.+c3*t/(a*p))
         ! 2.17... = 2*g/9, see (10-138)
         ! 3.0849d-5 = 1.26 * lamda_0 * P_0 / T_0 with a mistake over T_0: 273.15 instead of 293.15 in P & K
      else
         ! Regime 2
         best=c4*a**3*(rhow-rho_a)*rho_a/(eta*eta) ! equation (10-142) with 104.60267 = 32*g/3
         x=log(best)
         ! evaluation of BEARD-polynomial with Horner-scheme
         y=b6*x+b5
         y=y*x+b4
         y=y*x+b3
         y=y*x+b2
         y=y*x+b1
         y=y*x+b0
         vterm=eta*exp(y)/(2.*rho_a*a)
      end if

      end function vterm

!
!----------------------------------------------------------------
!

      subroutine wfield
! calculation of subsidence if chosen to be time dependent

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka

      implicit double precision (a-h,o-z)

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      zeit=lst*3600.+lmin*60.
      u0=dcos(2.*pi*zeit/86400.)
      do k=1,n
         w(k)=eta(k)/1000.*0.5*((wmax+wmin)+(wmin-wmax)*u0)
      enddo
      do k=n,1,-1
         w(k)=w(k)-w(1)
      enddo

      end subroutine wfield

!
!-------------------------------------------------------------------
!

      subroutine difm (dt)
! fully implicit procedure for the solution of the diffusion equations
! after Roache, 1972: Computational fluid dynamics, Appendix A.
! all quantities are similarly defined with a(k) --> xa(k), etc.
! except xd(k) which has another meaning than d(k) of Roache's program.
! dirichlet conditions at the surface and at the top of the atmosphere

      USE constants, ONLY :
! Imported Parameters:
     &     r0               ! Specific gas constant of dry air, in J/(kg.K)

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nm,
     &     nka

      implicit double precision (a-h,o-z)

      parameter (fcor=1.d-4)
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb53a/ thet(n),theti(n)
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /cb57/ xa(n),xb(n),xc(n),xd(n),xe(n),xf(n),oldu(n)
      dimension c(n)
      p21(tt)=610.7*dexp(17.15*(tt-273.15)/(tt-38.33))
      tke(1)=dmax1(1.d-06,3.2537*ustern**2)
      do k=1,n
         rho(k)=p(k)/(r0*(t(k)*(1.+0.61*xm1(k))))
         theta(k)=t(k)*thet(k)
         tke(k)=dmax1(1.d-05,tke(k)+tkep(k)*dt)
         c(k)=w(k)*dt/deta(k)
      enddo
! exchange coefficients
      call atk1 
! turbulent exchange with k_m:
      xa(1)=atkm(1)*dt/(detw(1)*deta(1))
      xe(1)=0.
      do k=2,nm
         xa(k)=atkm(k)*dt/(detw(k)*deta(k))
         xc(k)=xa(k-1)*detw(k-1)/detw(k)
         xb(k)=1.+xa(k)+xc(k)
         xd(k)=xb(k)-xc(k)*xe(k-1)
         xe(k)=xa(k)/xd(k)
      enddo
! u-component of horizontal wind field
      fdt=fcor*dt
      xf(1)=0.
      do k=2,nm
         oldu(k)=u(k)
         xf(k)=(u(k)+fdt*(v(k)-vg)+xc(k)*xf(k-1))/xd(k)
      enddo
      do k=nm,2,-1
         u(k)=xe(k)*u(k+1)+xf(k)
      enddo
! v-component of the horizontal wind field
      do k=2,nm
         xf(k)=(v(k)-fdt*(oldu(k)-ug)+xc(k)*xf(k-1))/xd(k)
      enddo
      do k=nm,2,-1
         v(k)=xe(k)*v(k+1)+xf(k)
      enddo
! turbulent kinetic energy
      xa(1)=atke(1)*dt/(detw(1)*deta(1))
      xe(1)=0.
      do k=2,nm
         xa(k)=atke(k)*dt/(detw(k)*deta(k))
         xc(k)=xa(k-1)*detw(k-1)/detw(k)
         xb(k)=1.+xa(k)+xc(k)
         xd(k)=xb(k)-xc(k)*xe(k-1)
         xe(k)=xa(k)/xd(k)
      enddo
      xf(1)=tke(1)
      do k=2,nm
         xf(k)=(tke(k)+xc(k)*xf(k-1))/xd(k)
      enddo
      do k=nm,2,-1
         tke(k)=xe(k)*tke(k+1)+xf(k)
      enddo
! turbulent exchange with k_h:
      xa(1)=atkh(1)*dt/(detw(1)*deta(1))
      do k=2,nm
         xa(k)=atkh(k)*dt/(detw(k)*deta(k))
         xc(k)=xa(k-1)*detw(k-1)/detw(k)
         xb(k)=1.+xa(k)+xc(k)
         xd(k)=xb(k)-xc(k)*xe(k-1)
         xe(k)=xa(k)/xd(k)
      enddo
! specific humidity
      xf(1)=xm1(1)
      do k=2,nm
         xf(k)=(xm1(k)+xc(k)*xf(k-1))/xd(k)
      enddo
      do k=nm,2,-1
         xm1(k)=xe(k)*xm1(k+1)+xf(k)
      enddo
! new temperature obtained from potential temperature
      xf(1)=theta(1)
      do k=2,nm
         xf(k)=(theta(k)+xc(k)*xf(k-1))/xd(k)
      enddo
      do k=nm,2,-1
         theta(k)=xe(k)*theta(k+1)+xf(k)
      enddo
! large scale subsidence
      do k=2,nm
         kp=k+1
         theta(k)=theta(k)-c(k)*(theta(kp)-theta(k))
         u(k)=u(k)-c(k)*(u(kp)-u(k))
         v(k)=v(k)-c(k)*(v(kp)-v(k))
         xm1(k)=xm1(k)-c(k)*(xm1(kp)-xm1(k))
         tke(k)=tke(k)-0.5*(c(k)+c(k+1))*(tke(kp)-tke(k))
      enddo
      do k=2,nm
         t(k)=theta(k)*theti(k)
         feu(k)=xm1(k)*p(k)/((0.62198+0.37802*xm1(k))*p21(t(k)))
      enddo

      end subroutine difm

!
!----------------------------------------------------------------------
!

      subroutine difp (dt)
! fully implicit procedure for the solution of the turbulent transport
! of aerosols and cloud droplets
! for further details see subroutine difm
! for diffusion mixing ratio is needed --> factor 1./am3(k,1)
! (#/cm^3 --> #/mol)

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nm,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)

      common /blck01/ am3(n),cm3(n)
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb45/ u(n),v(n),w(n)
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb57/ xa(n),xb(n),xc(n),xd(n),xe(n),xf(n),oldf(n)

      dimension c(n)

      xe(1)=0.
      xa(1)=atkh(1)*dt/(detw(1)*deta(1))
!      do k=2,nf
      do k=2,nm
         xa(k)=atkh(k)*dt/(detw(k)*deta(k))
         xc(k)=xa(k-1)*detw(k-1)/detw(k)*am3(k-1)/am3(k)
         xb(k)=1.+xa(k)+xc(k)
         xd(k)=xb(k)-xc(k)*xe(k-1)
         xe(k)=xa(k)/xd(k)
         c(k)=w(k)*dt/deta(k)
      enddo
!      if (lct.gt.1) then
         do ia=1,nka
         do jt=1,nkt
            xf(1)=ff(jt,ia,2)/am3(2)*1.d6
!            do k=2,nf
            do k=2,nm
               xf(k)=(ff(jt,ia,k)/am3(k)*1.d6+xc(k)*xf(k-1))/xd(k)
            enddo
!            do k=nf,2,-1
            do k=nm,2,-1
               ff(jt,ia,k)=(xe(k)*ff(jt,ia,k+1)/am3(k+1)*1.d6+xf(k))
     &              *am3(k)/1.d6
            enddo
!           large scale subsidence
            do k=2,nm
               kp=k+1
               ff(jt,ia,k)=ff(jt,ia,k)-c(k)*(ff(jt,ia,kp)-ff(jt,ia,k))
            enddo
         enddo
         enddo
!         do k=2,nf
         do k=2,n
            fsum(k)=0.
            do ia=1,nka
            do jt=1,nkt
               fsum(k)=fsum(k)+ff(jt,ia,k)
            enddo
            enddo
         enddo
!      else
!         xf(1)=fsum(2)
!         do k=2,nf
!            oldf(k)=fsum(k)
!            xf(k)=(fsum(k)+xc(k)*xf(k-1))/xd(k)
!         enddo
!         do k=nf,2,-1
!            fsum(k)=xe(k)*fsum(k+1)+xf(k)
!         enddo
!         do k=2,nf
!            x0=fsum(k)/oldf(k)
!            do ia=1,nka
!            do jt=1,nkt
!               ff(jt,ia,k)=ff(jt,ia,k)*x0
!            enddo
!            enddo
!         enddo
!      endif

      end subroutine difp

!
!-----------------------------------------------------------------
!

      subroutine difc (dt)
! fully implicit procedure for the solution of the turbulent transport
! of chemical species
! for further details see subroutine difm

! jjb work done: removal of unused arguments
!     missing declarations and implicit none

      USE gas_common, ONLY :
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
     &     n,
     &     nm,
     &     nkc

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision dt

! Local scalars:
      integer k, kc, kp, j
! Local arrays:
      double precision c(n)

! Common blocks:
      common /blck01/ am3(n),cm3(n)
      double precision am3, cm3

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      double precision atke, atkh, atkm, tke, tkep, buoy

      common /cb45/ u(n),v(n),w(n)
      double precision u, v, w

      common /cb57/ xa(n),xb(n),xc(n),xd(n),xe(n),xf(n),oldf(n)
      double precision xa, xb, xc, xd, xe, xf, oldf

      common /liq_pl/ nkc_l
      integer nkc_l
!- End of header ---------------------------------------------------------------

! calculation of exchange coefficients
      xa(1)=atkh(1)*dt/(detw(1)*deta(1))
      xe(1)=0.
      do k=2,nm
         xa(k)=atkh(k)*dt/(detw(k)*deta(k))
         xc(k)=xa(k-1)*detw(k-1)/detw(k)*am3(k-1)/am3(k)
         xb(k)=1.+xa(k)+xc(k)
         xd(k)=xb(k)-xc(k)*xe(k-1)
         xe(k)=xa(k)/xd(k)
         c(k)=w(k)*dt/deta(k) ! large scale subsidence
      enddo

! gase phase species
! s1, s3, sl1, sion1 in mol/m^3
! for diffusion mixing ratio is needed --> factor 1./am3(k)
! (mol/m^3 --> mol/mol)
      do j=1,j1
         xf(1)=s1(j,2)/am3(2)
         do k=2,nm
            xf(k)=(s1(j,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
         enddo
         do k=nm,2,-1
            s1(j,k)=(xe(k)*s1(j,k+1)/am3(k+1)+xf(k))*am3(k)
         enddo
!        large scale subsidence
         do k=2,nm
            s1(j,k)=s1(j,k)-c(k)*(s1(j,k+1)-s1(j,k))
         enddo
      enddo
! dito for radicals
      do j=1,j5
         xf(1)=s3(j,2)/am3(2)
         do k=2,nm
            xf(k)=(s3(j,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
         enddo
         do k=nm,2,-1
            s3(j,k)=(xe(k)*s3(j,k+1)/am3(k+1)+xf(k))*am3(k)
         enddo
!        large scale subsidence
         do k=2,nm
            s3(j,k)=s3(j,k)-c(k)*(s3(j,k+1)-s3(j,k))
         enddo
      enddo

! aqueous phase species
!      if (lct.lt.2) return
!      if (ndt.lt.2) return

      do kc=1,nkc_l
         do j=1,j2
            xf(1)=sl1(j,kc,2)/am3(2)
!            do k=2,nf
            do k=2,nm
!            do k=lcl,lct
!            xf(ndb-1)=sl1(j,kc,ndb)
!            xf(ndt+1)=sl1(j,kc,ndt)
!            do k=ndb,ndt
               xf(k)=(sl1(j,kc,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
            enddo
!            do k=nf,2,-1
            do k=nm,2,-1
!            do k=lct,lcl,-1
!            do k=ndt,ndb,-1
               sl1(j,kc,k)=(xe(k)*sl1(j,kc,k+1)/am3(k+1)+
     &              xf(k))*am3(k)
            enddo
!           large scale subsidence
!            do k=2,nf        
            do k=2,nm        
               kp=k+1
               sl1(j,kc,k)=sl1(j,kc,k)-c(k)*(sl1(j,kc,kp)-sl1(j,kc,k))
            enddo
         enddo
      enddo
! aqueous phase ion species
      do kc=1,nkc_l
         do j=1,j6
            xf(1)=sion1(j,kc,2)/am3(2)
!            do k=2,nf
            do k=2,nm
!            do k=lcl,lct
!            xf(ndb-1)=sion1(j,ndb,kc)
!            xf(ndt+1)=sion1(j,ndt,kc)
!            do k=ndb,ndt
               xf(k)=(sion1(j,kc,k)/am3(k)+xc(k)*xf(k-1))/xd(k)
            enddo
!            do k=nf,2,-1
            do k=nm,2,-1
!            do k=lct,lcl,-1
!            do k=ndt,ndb,-1
               sion1(j,kc,k)=(xe(k)*sion1(j,kc,k+1)/
     &              am3(k+1)+xf(k))*am3(k)
            enddo
!           large scale subsidence
!            do k=2,nf         
            do k=2,nm
               kp=k+1
               sion1(j,kc,k)=sion1(j,kc,k)-c(k)*
     &              (sion1(j,kc,kp)-sion1(j,kc,k))
            enddo
         enddo
      enddo

      end subroutine difc

!
!--------------------------------------------------------------------
!

      subroutine atk0
! calculation of exchange coefficients, mixing length etc at model start

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka

      implicit double precision (a-h,o-z)

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho

! mixing length
      xl(1)=0.
      x1=(ug+vg)*2.7
      do k=2,n
         x2=0.4*etw(k)
         xl(k)=x2*x1/(x2+x1)
         xl(k)=dmin1(xl(k),deta(k))
      enddo
      atkm(1)=0.5*eta(2)*ustern/gclu
      atkh(1)=0.5*eta(2)*ustern/gclt
      do k=2,n-1
         vh=((u(k+1)-u(k))**2+(v(k+1)-v(k))**2)/deta(k)**2
         zz=etw(k)+z0
         x0=(0.4*zz/(1.+0.4*zz/xl(k)))**2
         st=g*(theta(k+1)-theta(k))/(deta(k)*theta(k))
         if (st.le.0.) then
! unstable case and neutral case
            atkm(k)=x0*sqrt(vh-11.*st)
            if ((vh-3.*st).eq.0.) then
               atkh(k)=atkm(k)
            else
               atkh(k)=1.35*atkm(k)*(vh-5.5*st)/(vh-3.*st)
            endif
         else
! stable case
            atkm(k)=x0*vh/sqrt(vh+6.*st)
            atkh(k)=1.35*atkm(k)*vh/(vh+6.*st)
         endif
         atkm(k)=dmax1(1.d-03,atkm(k))
         atkh(k)=dmax1(1.d-03,atkh(k))
      enddo
      atkm(n)=0.
      atkh(n)=0.

      end subroutine atk0

!
!-------------------------------------------------------------
!

      subroutine atk1

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nm,
     &     nka

! jjb work done
!     - removal of unused parameters, remaining one now in modules

! jjb questions
!     - why lct.le.lcl+2, not +1, or not lct.ne.lcl ?
!     - why hardcoded values in several loops: k=10,n (why 10, why n not nm), then lct-4,lct+4

      implicit double precision (a-h,o-z)
! turbulent exchange coefficients after 2.5 level model of Mellor
! and Yamada, JAS 1974.

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb42a/ tkeps(n),tkepb(n),tkepd(n)
      common /cb43/ gm(n),gh(n),sm(n),sh(n),xl(n)
      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb53a/ thet(n),theti(n)
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /kinv_i/ kinv
      dimension es(n),xlo(n),dmw(n),dthetl(n),xmw(n)
! statement functions
      qsatur(esat,ppp)=0.62198*esat/(ppp-0.37802*esat)
      vapsat(ttt)=610.7*dexp(17.15*(ttt-273.15)/(ttt-38.33))
      if (lct.le.lcl+2) then
!
! cloud free situation
!
! buoyancy term
         do k=2,nm
            x0=((1+0.61*xm1(k))*(theta(k+1)-theta(k))
     &         +0.61*theta(k)*(xm1(k+1)-xm1(k)))/deta(k)
            sm(k)=x0
            sh(k)=x0
            buoy(k)=0.8*buoy(k)+0.2*x0                    !filter
            thetl(k)=(1+0.61*xm1(k))*theta(k)             !for output only
         enddo
! calculate inversion level
         do k=10,n
            kinv=k
            if (buoy(k).gt.1.d-5) go to 2000
         enddo
 2000    continue
      else
!
! cloudy situation
!
! liquid water potential temperature
! l21/cp=2.4774d06/1005=2465.1
! l21/r0=2.4774d06/287.05=8630.6
! l21/r1=2.4774d06/461.51=5368.
! thet(k)=theta(k)/t(k); theti(k)=t(k)/theta(k)
         do k=1,nm
            thetl(k)=theta(k)-2465.1*thet(k)*xm2(k)/rho(k)
            xmw(k)=xm1(k)+xm2(k)/rho(k)
         enddo
         thetl(n)=thetl(nm)+1.
! buoyancy term
         do k=2,nm
            dthetl(k)=(thetl(k+1)-thetl(k))/deta(k)
            dmw(k)=(xmw(k+1)-xmw(k))/deta(k)
            x0=(1.+0.61*xmw(k))*dthetl(k)+0.61*thetl(k)*dmw(k)
            sh(k)=0.8*sh(k)+0.2*x0                        !filter
         enddo
         do k=2,lct-1
            ql=xm2(k)/rho(k)
            esat=vapsat(t(k))
            qsat=qsatur(esat,p(k))
            qslt=5368.*qsat/(t(k)*t(k))
            xa=1./(1.+2465.1*qslt)
            xb=xa*theti(k)*qslt
            betat=(1.+0.61*xm1(k)-ql)
            betaw=0.61*(thetl(k)+2465.1*thet(k)*ql)
            betal=(1.+0.61*xmw(k)-3.22*ql)*2465.1*thet(k)-1.61*thetl(k)
            x0=(betat-xb*betal)*dthetl(k)+(betaw+xa*betal)*dmw(k)
            sm(k)=0.8*sm(k)+0.2*x0                         !filter
            x1=60.*(dmin1(feu(k),1.d0)-1.)                 !calc. of alpha
            betal=betal*dexp(x1)
            x0=(betat-xb*betal)*dthetl(k)+(betaw+xa*betal)*dmw(k)
            buoy(k)=0.8*buoy(k)+0.2*x0                     !filter
         enddo
         do k=lct,nm
            buoy(k)=sh(k)
         enddo
! calculate inversion level
         kinv=lct+5
         do k=lct-4,lct+4
            if (buoy(k).gt.1.d-5) then
            kinv=min(kinv,k)
            endif
         enddo
         kinv=kinv-1
         print *,'cloud kinv ',kinv
      endif
         print *,'kinv : ',kinv
! mixing length
      do k=1,n
         xlo(k)=xl(k)
         es(k)=sqrt(2.*tke(k))
      enddo
      x0=0.
      x1=0.
      do k=2,kinv-1
         x2=es(k)*deta(k)
         x0=x0+x2*etw(k)
         x1=x1+x2
      enddo
      x2=x0/x1
      zinv=etw(kinv)
      x4=0.1-detw(kinv)/x2
      xl(1)=0.
      do k=2,kinv-1
         x0=0.4*etw(k)
         x1=dmax1(detw(k),x2*(0.1-x4*dexp((etw(k)-zinv)/15.)))
         xl(k)=x0*x1/(x0+x1)
      enddo
      do k=kinv,n
         x0=0.4*etw(k)
         x1=detw(k)
         xl(k)=x0*x1/(x0+x1)
      enddo
      do k=2,nm
         xl(k)=0.95*xlo(k)+0.05*xl(k)                    !filter
      enddo
! stability functions gh and gm
      gh(1)=0.
      gm(1)=0.
      do k=2,nm
         x1=xl(k)*xl(k)/(es(k)*es(k))
         ghn=-g*x1/theta(k)*buoy(k)
         ghn=dmin1(ghn,0.03d0)
         gh(k)=dmax1(ghn,-0.6d0)
         gmn=x1*((u(k+1)-u(k))**2+(v(k+1)-v(k))**2)/(deta(k)*deta(k))
         gm(k)=dmin1(gmn,25.d0*(0.03-gh(k)))
      enddo
! exchange coefficients
! turbulent kinetic energy production/dissipation
      atkh(1)=0.5*eta(2)*ustern/gclt
      atkm(1)=0.5*eta(2)*ustern/gclu
      atke(1)=atkm(1)
      do k=2,nm
         ghn=gh(k)
         gmn=gm(k)
         x0=1./(1.+ghn*(-36.7188+187.4408*ghn)
     &      +gmn*(5.0784-88.83949*ghn))
         smn=(.6992d0-9.339487d0*ghn)*x0
         shn=(.74d0-4.534128d0*ghn+.901924d0*gmn)*x0
         x1=es(k)**3/xl(k)
         tkeps(k)=x1*smn*gmn
         tkepb(k)=x1*shn*ghn
         tkepd(k)=-x1/16.6
         tkep(k)=tkeps(k)+tkepb(k)+tkepd(k)
         x2=es(k)*xl(k)
         atkh(k)=x2*shn
         atkm(k)=x2*smn
         atke(k)=dmin1(atkm(k),x2*0.2)
      enddo
      do k=1,nm
         atke(k)=0.5*(atke(k)+atke(k+1))
      enddo

      end subroutine atk1

!
!-------------------------------------------------------------
!

      subroutine soil (dt)

      USE constants, ONLY :
! Imported Parameters:
     &     rhow                  ! Water density [kg/m**3]

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nb,
     &     nbm,
     &     nka

      implicit double precision (a-h,o-z)
! fully implicit procedure for the solution of the diffusive heat
! and moisture transport within the soil.
! for further details see subroutine difm.

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb57/ xa(n),xb(n),xc(n),xd(n),xe(n),xf(n),oldu(n)
! soil temperature
      xe(1)=0.
      x0=dmax1(eb(1),ebc)
      akb=anu0*x0**bs0/((1.-ebs)*rhoc+eb(1)*rhocw)
      xa(1)=akb*dt/(dzbw(1)*dzb(1))
      do 1000 k=2,nbm
      x0=dmax1(eb(k),ebc)
      akb=anu0*x0**bs0/((1.-ebs)*rhoc+eb(k)*rhocw)
      xa(k)=akb*dt/(dzbw(k)*dzb(k))
      xc(k)=xa(k-1)*dzbw(k-1)/dzbw(k)
      xb(k)=1.+xa(k)+xc(k)
      xd(k)=xb(k)-xc(k)*xe(k-1)
 1000 xe(k)=xa(k)/xd(k)
      xf(1)=tb(1)
      do 1010 k=2,nbm
 1010 xf(k)=(tb(k)+xc(k)*xf(k-1))/xd(k)
      do 1020 k=nbm,2,-1
 1020 tb(k)=xe(k)*tb(k+1)+xf(k)
! volumetric moisture content
      x0=2.*bs+3.
      x1=bs+2.
      x2=-bs*aks*psis/ebs
      do 1030 k=2,nbm
      x3=(eb(k)+dzbw(k)*(eb(k+1)-eb(k))/(2.*dzb(k)))/ebs
      ak(k)=aks*x3**x0
 1030 d(k)=x2*x3**x1
      ak(1)=0.
      d(1)=0.
      if (abs(eb(2)-eb(1)).gt.1.d-5)
     & d(1)=ajm*dzb(1)/(rhow*(eb(2)-eb(1)))
      xa(1)=d(1)*dt/(dzbw(1)*dzb(1))
      do 1040 k=2,nbm
      xa(k)=d(k)*dt/(dzbw(k)*dzb(k))
      xc(k)=xa(k-1)*dzbw(k-1)/dzbw(k)
      xb(k)=1.+xa(k)+xc(k)
      xd(k)=xb(k)-xc(k)*xe(k-1)
 1040 xe(k)=xa(k)/xd(k)
      xf(1)=eb(1)
      do 1050 k=2,nbm
 1050 xf(k)=(eb(k)+dt/dzbw(k)*(ak(k-1)-ak(k))+xc(k)*xf(k-1))/xd(k)
      do 1060 k=nbm,2,-1
 1060 eb(k)=xe(k)*eb(k+1)+xf(k)

      end subroutine soil

!
!-------------------------------------------------------------
!

      subroutine surf0 (dt)

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka

      implicit double precision (a-h,o-z)
! lower boundary condition for water surface
! constant temperature and saturation specific humidity

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)

!      tw=tw-5.787d-6*dt
!      tw=tw-6.94444d-6*dt
      t(1)=tw
      pp21=610.7*dexp(17.15*(tw-273.15)/(tw-38.33))
!      xm1(1)=0.62198*pp21/(p(1)-0.37802*pp21)
!      xm1(1)=0.75*0.62198*pp21/(p(1)-0.37802*pp21) !# INDOEX
!      xm1(1)=0.75*0.62198*pp21/(p(1)-0.37802*pp21) !# Appledore
!      xm1(1)=0.7*0.62198*pp21/(p(1)-0.37802*pp21)  !# aerosol nosub run
!      xm1(1)=0.6*0.62198*pp21/(p(1)-0.37802*pp21)  !# aerosol nosub run
!      xm1(1)=0.55*0.62198*pp21/(p(1)-0.37802*pp21)  !# aerosol nosub run
!      xm1(1)=0.3*0.62198*pp21/(p(1)-0.37802*pp21)  !# aerosol nosub run
      xm1(1)=0.8*0.62198*pp21/(p(1)-0.37802*pp21) !# cloud nosub run
!      xm1(1)=0.982*0.62198*pp21/(p(1)-0.37802*pp21)  !# cloud sub run Raoult's law
!      xm1(1)=0.7*0.62198*pp21/(p(1)-0.37802*pp21)  !# aerosol sub run
      uu=u(2)
      vv=v(2)
      vqr=uu*uu+vv*vv
      vbt=sqrt(vqr)
      zp=0.5*eta(2)+z0
      zpdz0=dlog(zp/z0)
      xnvl=g*(theta(2)-tw)*2./(theta(2)+tw)
      zpdl=zp*xnvl/vqr
      call claf (zpdl,zpdz0,cu,ctq)
      ustern=dmax1(0.01d0,vbt/cu)
! charnock's relation: z0=0.015*ustern**2/g
      z0=0.015*ustern*ustern/g
      gclu=cu
      gclt=ctq

      end subroutine surf0

!
!-------------------------------------------------------------
!

      subroutine surf1 (dt)

      USE constants, ONLY :
! Imported Parameters:
     &     cp,                   ! Specific heat of dry air, in J/(kg.K)
     &     r1,                   ! Specific gas constant of water vapour, in J/(kg.K)
     &     rhow                  ! Water density [kg/m**3]


      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nb,
     &     nka

      implicit double precision (a-h,o-z)
! calculation of surface temperature and volumetric moisture content
! by means of balance of fluxes at the surface following Pielke,
! mesoscale meteorological modelling, chapter 11.
      logical l1

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb46/ ustern,gclu,gclt
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      dimension ebb(40),tss(40),ftss(40),fqss(40)
! Stefan-Boltzmann-constant
      parameter (sigma=5.6697d-8)
      parameter (al31=2.835d+6)
      parameter (t0=273.15)
      p21(tt)=610.7*dexp(17.15*(tt-273.15)/(tt-38.33))
      al21(tt)=3.1387818d+06-2335.5*tt
      cm(pp)=0.62198*pp/(ps-0.37802*pp)
      rrho=rho(1)
      uu=u(2)
      vv=v(2)
      vqr=dmax1(uu*uu+vv*vv,1.d-12)
      vbt=sqrt(vqr)
!     gamma=g/cp
      bs3=2.*bs+3.
!     bs2=bs+2.
      psi2=psis*(ebs/eb(2))**bs
      qq2=xm1(2)
      ts=t(1)
      ps=p(1)
      eb1=eb(1)
      l1=.false.
      zp=deta(1)+z0
      zpdz0=dlog(zp/z0)
      xnvl=g*(theta(2)-ts)*2./(theta(2)+ts)
      zpdl=zp*xnvl/vqr
      call claf (zpdl,zpdz0,cu,ctq)
      ust=dmax1(0.01d0,vbt/cu)
! specific humidity at the surface
      if (ts.ge.t0) then
      xm21s=cm(p21(ts))
      else
      xm21s=cm(p31(ts))
      endif
      psi1=psis*(ebs/eb1)**bs
      qs=xm21s*dexp(g*psi1/(r1*ts))
! tstern,qstern
      tst=(theta(2)-ts*(1.+0.608*qs))/ctq
      qst=(qq2-qs)/ctq
! heatflux from ground
      x1=dmax1(eb1,ebc)
      anu=anu0*x1**bs0
      ajb=anu*(tb(2)-ts)/dzb(1)
! microturbulent flux of water vapor
      ajq=rrho*ust*qst
! latent microturbulent enthalpy flux
! ajs is water flux due to droplet sedimentation;
      if (ts.lt.t0) then
      ajl=al31*ajq-(al31-al21(ts))*ajs
      else
      ajl=al21(ts)*ajq
      endif
! sensible microturbulent enthalpy flux
      ajt=rrho*cp*ust*tst
! ground moisture flux
!     rak1=rhow*aks*(eb1/ebs)**bs3
      xa=0.5
      rak1=rhow*aks*((xa*eb1+(1.-xa)*eb(2))/ebs)**bs3
      ajm=rak1*((psi2-psi1)/dzb(1)-1.)
!     debs=-bs*aks*psis/ebs
!     ajm=(debs*((eb1+eb(2))/(2.*ebs))**bs2*(eb(2)-eb1)/dzb(1)-
!    & aks*((eb1+eb(2))/(2.*ebs))**bs3)*rhow
! dew flux
      ajd=0.
      x0=ajq+ajm+ajs
      if (eb1.lt.ebs) go to 2000
      ajd=-x0
      ddew=tau/dt
      if (x0.lt.0.) ajd=dmin1(-x0,ddew)
      ddew=ddew-ajd
! flux balances
 2000 fts=sl+sk+ajb+ajl+ajt-sigma*ts**4
      fqs=x0+ajd
      if (l1) goto 2010
! 2-dimensional newton-raphson iteration for ts and eb(1)
      do 1000 iaa=1,20
      ebb(iaa)=eb1
      ftss(iaa)=fts
      fqss(iaa)=fqs
      tss(iaa)=ts
!     fts0=fts
!     fqs0=fqs
      ts0=ts
      eb10=eb1
! calculation of derivatives d(fluxes)/d(eta) and d(fluxes)/d(T)
! djbde=d(ajb)/d(eta); djbdt=d(ajb)/d(T) etc.
      djbde=0.
      if (eb1.gt.ebc) djbde=ajb*bs0/eb1
      djbdt=-anu/dzb(1)
      djqde=rrho*ust*qs*g*bs*psi1/(ctq*r1*ts*eb1)
      x0=p21(ts)
      djqdt=rrho*ust*qs/ctq*(g*psi1/(r1*ts*ts)+
     & x0*4027.163/((x0-.37802*ps)*(ts-38.33)**2))
      djtdt=-rrho*cp*ust/ctq
!     djmde=(ajm*bs3+rak1/dzb(1)*psi1*bs)/eb1
      djmde=rak1/dzb(1)*psi1*bs/eb1
!     djmde=(debs*(eb1/ebs)**bs2/dzb(1)*((eb(2)-eb1)*bs2/eb1-1.)-
!    & bs3*aks/eb1*(eb1/ebs)**bs3)*rhow
!     djmde=dmin1(djmde,0.)
! coefficients for solution of linear equation system
      x0=al21(ts)
      f1e=djbde+x0*djqde
      f1t=djbdt-2335.5*ajq+x0*djqdt+djtdt-4.*sigma*ts*ts*ts
      f2e=djqde+djmde
      f2t=djqdt
      det=f1e*f2t-f1t*f2e
!      if (abs(det).lt.1.d-10) x0=sqrt(1.d0-2.d0)
      if (abs(det).lt.1.d-10) stop 'SR surf1'
! new values of ts and eb1
      ts=ts+(fts*f2e-fqs*f1e)/det
      eb1=eb1+(fqs*f1t-fts*f2t)/det
      eb1=dmin1(eb1,ebs)
      eb1=dmax1(eb1,ebs/15.d0)
      if (ddew.gt.0.) eb1=ebs
      if (ts.gt.300..or.ts.lt.250.) ts=ts0-0.01
! new surface fluxes
      if (ts.ge.t0) then
      xm21s=cm(p21(ts))
      else
      xm21s=cm(p31(ts))
      endif
      psi1=psis*(ebs/eb1)**bs
      qs=xm21s*dexp(g*psi1/(r1*ts))
      qst=(qq2-qs)/ctq
      tst=(theta(2)-ts*(1.+0.608*qs))/ctq
      x1=dmax1(eb(1),ebc)
      anu=anu0*x1**bs0
      ajb=anu*(tb(2)-ts)/dzb(1)
      ajq=rrho*ust*qst
      if (ts.lt.t0) then
      ajl=al31*ajq-(al31-al21(ts))*ajs
      else
      ajl=al21(ts)*ajq
      endif
      ajt=rrho*cp*ust*tst
      ajm=rak1*((psi2-psi1)/dzb(1)-1.)
      ajd=0.
      x0=ajq+ajm+ajs
      if (eb1.lt.ebs) go to 2020
      ajd=-x0
      ddew=tau/dt
      if (x0.lt.0.) ajd=dmin1(-x0,ddew)
      ddew=ddew-ajd
! flux balances
 2020 fts=sl+sk+ajb+ajl+ajt-sigma*ts**4
      fqs=x0+ajd
! convergence criteria
      if(dabs(ts-ts0).le.1.d-2.and.dabs(eb1-eb10).le.1.d-3)go to 2030
      if (dabs(fts).le.1.d-1.and.dabs(fqs).le..1*dabs(ajq)) goto 2030
 1000 continue
      write (6,6000) eb1,ts,fts,fqs
 6000 format (10x,'no convergence of ts- and eb1-iteration:'/
     & 'eb1',f16.4,'ts',f16.4,'fts',f16.4,'fqs',f16.4)
      write (6,6010) (ebb(i),tss(i),ftss(i),fqss(i),i=1,20)
 6010 format (10x,3f16.4,e16.4)
 2030 if (tau.gt.0..and.ts.lt.t0) l1=.true.
      if (ts.gt.t0.and.reif.gt.0) l1=.true.
      if (.not.l1) go to 2010
      ts=t0
! change of dew or rime by evaporation/condensation
 2010 continue
      if (ts.ge.t0) tau=tau-ajd*dt
      if (ts.lt.t0) reif=reif-ajd*dt
      if (.not.l1) goto 2040
! coexistence of dew and rime
      uwr=dmin1(dt*fts/3.35d+5,reif)
      uwr=dmax1(uwr,-tau)
      tau=tau+uwr
      reif=reif-uwr
 2040 tau=dmax1(0.d0,tau)
      reif=dmax1(0.d0,reif)
! storage of surface values
      t(1)=ts
      tb(1)=ts
      xm1(1)=qs
      feu(1)=qs/xm21s
      eb(1)=eb1
      xnvl=g*(theta(2)-ts)*2./(theta(2)+ts)
      zpdl=zp*xnvl/vqr
      call claf (zpdl,zpdz0,cu,ctq)
      ust=dmax1(0.01d0,vbt/cu)
      ustern=ust
      gclu=cu
      gclt=ctq

      end subroutine surf1

!
!-------------------------------------------------------------
!

      function p31(t)
! water vapour pressure over ice
      implicit double precision (a-h,o-z)
      parameter (t1=273.16)
      xlog10=-9.09685*(t1/t-1.)-3.56654*dlog10(t1/t)
     & +0.87682*(1.-t/t1)+0.78614
      p31=100.*(10.**xlog10)

      end function p31

!
!-------------------------------------------------------------
!

      subroutine claf (zpdl,zpdz0,u,tq)
!
! Description:
!    interpolation of clarke functions by means of tabulated values
!    u: clarke function for momentum; tq: for temperature, humidity etc.
!

!
! History:
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.2      01/2017   Removed labels and goto's                     <Josue Bock>
!                    Use Fortran generic functions min and max
! 1.1      08/2016   Use Andreas Bott's updated version from str   <Josue Bock>
!                    Header
!                    Declaration of all variables and imput none
!
! 1.0       ?        Original code used in Mistra v741             <Andreas Bott>
!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!
! Declarations:

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision zpdl, zpdz0
! Scalar arguments with intent(out):
      double precision u, tq

! Local scalars:
      double precision dx, dy
      double precision zpdla, zpdz0a
      integer i,nl,nz

! Common blocks:
      common /cb61/ fu(18,7),ft(18,7),xzpdl(18),xzpdz0(7) ! Clarks table data
      double precision fu, ft, xzpdl, xzpdz0
!- End of header ---------------------------------------------------------------

! xzpdl tabled values range from -5.5 to 3.0. Here, zpdla is forced to be within this range
      zpdla=dmax1(zpdl,-5.5d0)
      zpdla=dmin1(zpdla,3.d0)
! xzpdz0 tabled values range from 3. to 17. zpdz0a is forced to be lower equal to the max value.
      zpdz0a=dmin1(zpdz0,17.d0)

! Find nl and nz: indexes of the tabled values greater or equal to the actual values
      do i=2,18                     ! note that zpdla > xzpdl(1), thus nl >= 2
         nl=i                       !  thus no out-of-bounds problem with index nl-1 in dx (see below)
         if (zpdla.lt.xzpdl(i)) exit
      enddo

      do i=1,7
         nz=i
         if (zpdz0a.lt.xzpdz0(i)) exit
      enddo

! dx: proportionality factor
      dx=(zpdla-xzpdl(nl-1))/(xzpdl(nl)-xzpdl(nl-1))

      if (nz.eq.1) then
         dy=zpdz0a/xzpdz0(1)
         u=(fu(nl,1)*dx+fu(nl-1,1)*(1.-dx))*dy
         if (zpdl.ge.0.) then
            tq=u/1.35
         else
            tq=(ft(nl,1)*dx+ft(nl-1,1)*(1.-dx))*dy/1.35
         endif
      else
         dy=(zpdz0a-xzpdz0(nz-1))/(xzpdz0(nz)-xzpdz0(nz-1))
         u=fu(nl-1,nz-1)+(fu(nl,nz-1)-fu(nl-1,nz-1))*dx+(fu(nl-1,nz)
     &     -fu(nl-1,nz-1))*dy+(fu(nl,nz)-fu(nl-1,nz)+fu(nl-1,nz-1)
     &     -fu(nl,nz-1))*dx*dy
         if (zpdl.ge.0.) then
            tq=u/1.35
         else
            tq=(ft(nl-1,nz-1)+(ft(nl,nz-1)-ft(nl-1,nz-1))*dx
     &         +(ft(nl-1,nz)-ft(nl-1,nz-1))*dy+(ft(nl,nz)-ft(nl-1,nz)
     &         +ft(nl-1,nz-1)-ft(nl,nz-1))*dx*dy)/1.35
         endif
      end if

      end subroutine claf

!
!-------------------------------------------------------------
!

      subroutine kon (dt,chem)

! diffusional droplet growth by condensation

!     jjb work done
!      removed one unused parameter; use module instead
!      removed k0=k, used only in equil argument. For the sake of clarity
      

      USE constants, ONLY :
! Imported Parameters:
     & pi

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

      logical chem!,chmic ! jjb defined below, but unused
      common /cb11/ totrad (mb,nrlay)
      double precision totrad

      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb48/ sk,sl,dtrad(n),dtcon(n)
      double precision sk, sl, dtrad, dtcon

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /cb60/ ffk(nkt,nka),totr(mb),dfdt,feualt,pp,to,tn,
     &              xm1o,xm1n,kr
      double precision ffk, totr, dfdt, feualt,pp,to,tn,xm1o,xm1n
      integer kr

      common /blck06/ kw(nka),ka
      common /blck07/ part_o_a(nka,n),part_o_d(nka,n),
     &                part_n_a(nka,n),part_n_d(nka,n),pntot(nkc,n)
      common /blck08/ vol1_a(nka,n),vol1_d(nka,n),vol2(nkc,n)
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
      !ommon /liq_pl/ nkc_l
      dimension potot(nkc,n)
      do 1000 k=2,nf+1
! initialize
         do ia=1,nka
            vol1_a(ia,k)=0.               
            vol1_d(ia,k)=0.               
            part_o_a(ia,k)=0.
            part_o_d(ia,k)=0.
            part_n_a(ia,k)=0.
            part_n_d(ia,k)=0.
         enddo
!        do kc=1,nkc_l ! jjb no reason to initialise only until nkc_l
         do kc=1,nkc
            vol2(kc,k)=0.
            potot(kc,k)=0.
            pntot(kc,k)=0.
         enddo
         dtcon(k)=0.
         tn=t(k)
         xm1n=xm1(k)
!         xm2n=xm2(k) ! jjb variable unreferenced
         pp=p(k)
         if (feu(k).ge.0.7) goto 2000
         p21=610.7*dexp(17.15*(tn-273.15)/(tn-38.33))
         feu(k)=xm1n*pp/((0.62198+0.37802*xm1n)*p21)
         call equil (1,k)
         go to 1000
 2000    continue
         dfdt=dfddt(k)
         to=talt(k)
         xm1o=xm1a(k)
         feualt=feu(k)
         do ib=1,mb
            totr(ib)=totrad(ib,k-1)
         enddo
         kr=nar(k)
         do ia=1,nka
            do jt=1,nkt
               ffk(jt,ia)=ff(jt,ia,k)
            enddo
         enddo
!         chmic=chem.and.(xm2(k).gt.1.d-05)
! aerosol chemistry:
!      go to 2010
! aerosol + cloud chemistry:
! vol2 and vol1 old liquid volume in class (1-nkc) and row/class (1-nkc,1-nka)
! (um^3/cm^3), used in SR konc to shift moles from aerosol to drop
! and vice versa. 
! part_o and part_n old and new part. conc. in row/class (1-nkc,1-nka) (cm^-3)
! take care if non-soluble parts are in aerosol
!      if (.not.chmic) go to 2010
         if (.not.chem) go to 2010
!         xpi=4./3.*3.1415927
         xpi=4./3.*pi
! small (dry) aerosol: aerosol (1) and droplet (3)
         do ia=1,ka
            do jt=1,kw(ia)
               vol1_a(ia,k)=vol1_a(ia,k)+ffk(jt,ia)*xpi*rq(jt,ia)**3
               part_o_a(ia,k)=part_o_a(ia,k)+ffk(jt,ia)
            enddo
            do jt=kw(ia)+1,nkt
               vol1_d(ia,k)=vol1_d(ia,k)+ffk(jt,ia)*xpi*rq(jt,ia)**3
               part_o_d(ia,k)=part_o_d(ia,k)+ffk(jt,ia)
            enddo
            vol2(1,k)=vol2(1,k)+vol1_a(ia,k)
            vol2(3,k)=vol2(3,k)+vol1_d(ia,k)
            potot(1,k)=potot(1,k)+part_o_a(ia,k)
            potot(3,k)=potot(3,k)+part_o_d(ia,k)
         enddo
! large (dry) aerosol: aerosol (2) and droplet (4)
         do ia=ka+1,nka
            do jt=1,kw(ia)
               vol1_a(ia,k)=vol1_a(ia,k)+ffk(jt,ia)*xpi*rq(jt,ia)**3
               part_o_a(ia,k)=part_o_a(ia,k)+ffk(jt,ia)
            enddo
            do jt=kw(ia)+1,nkt
               vol1_d(ia,k)=vol1_d(ia,k)+ffk(jt,ia)*xpi*rq(jt,ia)**3
               part_o_d(ia,k)=part_o_d(ia,k)+ffk(jt,ia)
            enddo
            vol2(2,k)=vol2(2,k)+vol1_a(ia,k)
            vol2(4,k)=vol2(4,k)+vol1_d(ia,k)
            potot(2,k)=potot(2,k)+part_o_a(ia,k)
            potot(4,k)=potot(4,k)+part_o_d(ia,k)
         enddo
! old version:
!      if (.not.chmic) go to 2010
!      ap1o(k)=0.
!      ap2o(k)=0.
!      do ia=1,nka
!         apo(k,ia)=0.
!         do jt=1,kg
!            apo(k,ia)=apo(k,ia)+ffk(jt,ia)
!         enddo
!         do jt=kgp,kd
!            ap1o(k)=ap1o(k)+ffk(jt,ia)
!         enddo
!         do jt=kdp,nkt
!            ap2o(k)=ap2o(k)+ffk(jt,ia)
!         enddo
!      enddo
 2010    continue
! condensation and evaporation of droplets
! to: temperature after condensation of last timestep;
! tn: temperature before condensation of current timestep;
! tn = to + diffusion, radiation etc.
! same with xm1o and xm1n
         call subkon (dt)
! new values
         t(k)=to
         talt(k)=to
         xm1(k)=xm1o
         xm1a(k)=xm1o
         p21=610.7*dexp(17.15*(to-273.15)/(to-38.33))
         feu(k)=xm1o*pp/((0.62198+0.37802*xm1o)*p21)
         dfddt(k)=(feu(k)-feualt)/dt
         xm2(k)=0.
         do ia=1,nka
            do jt=1,nkt
               xm2(k)=xm2(k)+ffk(jt,ia)*e(jt)
               ff(jt,ia,k)=ffk(jt,ia)
            enddo
         enddo
!         xl21=3138708.-2339.4*tn ! jjb variable unreferenced
         dtcon(k)=(to-tn)/dt
! aerosol chemistry:
!      go to 1000
! aerosol + cloud chemistry:
         if (.not.chem) go to 1000
! small (dry) aerosol: aerosol (1) and droplet (3)
         do ia=1,ka
            do jt=1,kw(ia)
               part_n_a(ia,k)=part_n_a(ia,k)+ffk(jt,ia)
            enddo
            do jt=kw(ia)+1,nkt
               part_n_d(ia,k)=part_n_d(ia,k)+ffk(jt,ia)
            enddo
            pntot(1,k)=pntot(1,k)+part_n_a(ia,k)
            pntot(3,k)=pntot(3,k)+part_n_d(ia,k)
         enddo
! large (dry) aerosol: aerosol (2) and droplet (4)
         do ia=ka+1,nka
            do jt=1,kw(ia)
               part_n_a(ia,k)=part_n_a(ia,k)+ffk(jt,ia)
            enddo
            do jt=kw(ia)+1,nkt
               part_n_d(ia,k)=part_n_d(ia,k)+ffk(jt,ia)
            enddo
            pntot(2,k)=pntot(2,k)+part_n_a(ia,k)
            pntot(4,k)=pntot(4,k)+part_n_d(ia,k)
         enddo


! old version:
!      if (.not.chmic) go to 1000
!      ap2n(k)=0.
!      do ia=1,nka
!         apn(k,ia)=0.
!         do jt=1,kg
!            apn(k,ia)=apn(k,ia)+ffk(jt,ia)
!         enddo
!         do jt=kdp,nkt
!            ap2n(k)=ap2n(k)+ffk(jt,ia)
!         enddo
!      enddo
 1000 continue
! cloudy region: lowest (lcl) and highest (lct) cloud layer
      do k=nf+1,1,-1
         lct=k
         if (xm2(k).gt.1.d-05) go to 2030
      enddo
      lcl=1
      lct=1
      if (chem) call konc
      return

 2030 continue
      do k=1,lct
         lcl=k
         if (xm2(k).gt.1.d-05) go to 2040
      enddo
 2040 continue
! update chemical species      
! aerosol chemistry:
!      if (chem.and.lct.gt.lcl) call konc
! aerosol + cloud chemistry:
      if (chem) call konc

      end subroutine kon

!
!-------------------------------------------------------------
!

      subroutine equil (ncase,kk)

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)
! equilibrium values for radius rg(i) and water mass eg(i)
! of humidified aerosol particles at given relative humidity
! new distribution of the particles on their equilibrium positions

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      dimension rg(nka),eg(nka)

! get equilibrium distribution for 
!     case 1: active layer (1D: k<=nf (called from SR kon), 0d: active layer)
!     case 2: all layers above nf
      kmin=kk
      kmax=kk
      if (ncase.eq.2) then
         kmin=nf+1
         kmax=n
      endif

      do k=kmin,kmax
         do ia=1,nka
            do jt=2,nkt
               ff(1,ia,k)=ff(1,ia,k)+ff(jt,ia,k)
               ff(jt,ia,k)=0.
            enddo
         enddo
! equilibrium radius
         a0=a0m/t(k)
         do ia=1,nka
            b0=b0m(ia)*2.
! b0=b0m*rho3/rhow; rho3=2000; rhow=1000
            rg(ia)=rgl(rn(ia),a0,b0,feu(k))
            eg(ia)=4.d-09*pi/3.*(rg(ia)**3-rn(ia)**3)
         enddo
! new distribution
         do 1020 ia=1,nka
            do jt=1,nkt
               if (eg(ia).le.ew(jt)) then
                  if (jt.eq.1) goto 1020
                  ff(jt,ia,k)=ff(1,ia,k)
                  ff(1,ia,k)=0.
                  goto 1020
               endif
            enddo
 1020    continue
! total liquid water
         xm2(k)=0.
         do ia=1,nka
            do jt=1,nkt
               xm2(k)=xm2(k)+ff(jt,ia,k)*e(jt)
            enddo
         enddo
      enddo

      end subroutine equil

!
!-------------------------------------------------------------
!

      subroutine subkon (dt)
! core of microphysics: growth in 2D particle grid

! jjb work done:
!     - reindexed all (nka,nkt) arrays for computing efficiency
!     - defined water diffusivity as an external function

      USE constants, ONLY :
! Imported Parameters:
     &     cp,                   ! Specific heat of dry air, in J/(kg.K)
     &     pi,
     &     r0,                   ! Specific gas constant of dry air, in J/(kg.K)
     &     r1,                   ! Specific gas constant of water vapour, in J/(kg.K)
     &     rhow                  ! Water density [kg/m**3]

      USE global_params, ONLY :
! Imported Parameters:
     &     nka,
     &     nkt,
     &     mb

      implicit double precision (a-h,o-z)

      ! jjb declarations
      integer :: kr0

      double precision :: xdv ! diffusivity of water vapour in air in [m2/s]
      double precision :: xka ! thermal conductivity of air in [J/(m*s*K)]

      double precision, external :: diff_wat_vap
      double precision, external :: therm_conduct_air

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb49/ qabs(18,nkt,nka,3),qext(18,nkt,nka,3), ! only qabs is used here
     &              asym(18,nkt,nka,3)
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb51/ dlgew,dlgenw,dlne

      double precision :: psi(nkt),u(nkt)

      common /cb60/ ffk(nkt,nka),totr(mb),dfdt,feualt,p,t,tn,xm1,xm1n,kr
      double precision ffk, totr, dfdt, feualt, p, t, tn, xm1, xm1n
      integer kr

      dimension cd(nkt,nka),cr(nkt,nka),sr(nkt,nka),falt(nkt,nka),c(nkt)


! xl21 latent heat of evaporation; p21t water vapor pressure at satur.
! xl mean free path;
! deltav water vapor jump; deltat temperature jump
! xdvs,xkas diffusivity, conductivity corrected for small droplets
! xdv0, xka0 factors used for calculating xdvs and xkas.
! all terms in MKSA units only droplet mass and radius in mg and microns
! a0[microns]=2*sigma/(r1*rhow*t)*1.e6; sigma=76.1e-3=surface tension
! b0 solution term of koehler equation: b0=xnue*xmol2/xmol3*xa*ma/mw
! xnue: number of ions; xmol2 (xol3) mol masses of pure water (aerosols)
! xa volume fraction of soluble to unsoluble part of aerosol
! mw (ma) masses of water (aerosol nucleus)
! p21(tr)=p21(t)*exp(a0/r-b0*ma/mw)
! all formulas and constants after Pruppacher and Klett Chapter 13.
! droplet growth equation after Davies, J. Atmos. Sci., 1987
      p21(t)=610.7*exp(17.15*(t-273.15)/(t-38.33))
      xl21=3138708.-2339.4*t
      xldcp=xl21/cp

!      xka=4.38e-03+7.1e-05*t ! (13-18a) converted ?
      xka=therm_conduct_air(t)

!      xdv=4.0122e-05/p*t**1.94
      xdv=diff_wat_vap(t,p)

      xl=24.483*t/p ! (10-140 ?)
      deltav=1.3*xl
      deltat=2.7*xl
      rho=p/(r0*t*(1+0.61*xm1))
      rho21=p21(t)/(r1*t)
      rho21s=(xl21/(r1*t)-1.)*rho21/t
      a0=a0m/t
      xdv0=xdv*sqrt(2.*pi/(r1*t))/3.6e-08
      xka0=xka*sqrt(2.*pi/(r0*t))/(7.d-07*rho*cp)
      kr0=kr

      if (totr(1).lt.1.) then
         ib0=7                  ! solar bands excluded, IR only
      else
         ib0=1
      end if

      do ia=1,nka
         do jt=1,nkt
            jtp=min(jt+1,nkt)
            de0=dew(jt)
            dep=dew(jtp)
            de0p=de0+dep

            rk=rw(jt,ia)
            sr(jt,ia)=max(0.1d0,exp(a0/rk-b0m(ia)*en(ia)/ew(jt)))
            xdvs=xdv/(rk/(rk+deltav)+xdv0/rk)
            xkas=xka/(rk/(rk+deltat)+xka0/rk)
            x1=rhow*(xl21+xkas/(xdvs*rho21s*sr(jt,ia)))
            cd(jt,ia)=3.d+12*rho21*xkas/(x1*rk*rk*rho21s*sr(jt,ia))
            if (kr0.eq.3.and.rn(ia).lt.0.5) kr0=2
            rad=0.
            do ib=ib0,mb
               rad=rad+totr(ib)*(qabs(ib,jt,ia,kr0)*de0+
     &             qabs(ib,jtp,ia,kr0)*dep)/de0p
            enddo
            cr(jt,ia)=rad*7.5e05/(rk*x1)-rhow*4190.*(tn-t)/(dt*x1)
         enddo
      enddo
      falt(:,:) = ffk(:,:)

      feuneu=feualt+dfdt*dt
      if (feualt.lt.0.95)
     & feuneu=xm1n*p/(p21(tn)*(.62198+.37802*xm1n))
      fquer=0.5*(feuneu+feualt)
      ! Initialisation
      res = 0.d0 ! residue
      aa0=1./dt
      do itk=1,10
         dwsum=0.
         do ia=1,nka
            do jt=1,nkt
               psi(jt)=falt(jt,ia)
               c(jt)=(cd(jt,ia)*(fquer-sr(jt,ia))-cr(jt,ia))/dlne
            enddo

            u(1)=max(0.d0,c(1))
            do jt=2,nkt-1
               u(jt)=0.5*(c(jt)+abs(c(jt))+c(jt-1)-abs(c(jt-1)))
            enddo
            u(nkt)=min(0.d0,c(nkt-1))

            !print*,'------- New call------',itk,ia
            call advec (dt,u,psi)

            do jt=1,nkt
               ffk(jt,ia)=psi(jt) ! jjb
               dwsum=dwsum+(psi(jt)-falt(jt,ia))*e(jt)
            enddo
         enddo
         dmsum=dwsum/rho
         dtsum=xldcp*dmsum
         xm1=xm1n-dmsum
         t=tn+dtsum
         p1=xm1*p/(0.62198+0.37802*xm1)
         feuneu=p1/p21(t)
         resold=res
         res=feuneu+feualt-2.*fquer
!         if (feuneu.lt.0.95.or.abs(res).lt.1.d-06) return
         if (abs(res).lt.1.d-06) return
! calculation of fquer by Newton-Interpolation
         dres=res-resold
         aa=aa0
         if (itk.gt.1.and.abs(dres).gt.1.d-8) aa=(fqa-fquer)/dres
         fqa=fquer
         fquer=fquer+aa*res
      enddo
      write (6,6000)
 6000 format (10x,'no convergence of condensation iteration')

      end subroutine subkon

!
!-------------------------------------------------------------
!

      function diff_wat_vap(temperature,pressure)

      ! Diffusivity of water vapour in air
      ! for temperatures between -40 and +40 Celsius

      ! Dv = 0.211d-4 *(T/T0)**1.94 *(P0/P)

      ! where T0=273.15 K and P0=101325 Pa

      ! see Pruppacher and Klett, Microphysics of clouds and precipitations
      ! equation (13-3) p. 503, here expressed in MKS units


! Author
! ------
!     Josue Bock

! History
! -------
!     12-01-2017

      implicit none

      double precision :: diff_wat_vap

      double precision, intent(in) :: temperature  ! in [K]
      double precision, intent(in) :: pressure     ! in [Pa]

      double precision, parameter :: cst=0.211d-4  ! in [m2/s]
      double precision, parameter :: exponent=1.94
      double precision, parameter :: T0=273.15     ! in [K]
      double precision, parameter :: P0=101325     ! in [Pa]

      double precision, parameter :: cst2=cst*P0/(T0**exponent)

      diff_wat_vap = cst2 * temperature**exponent / pressure

      if (temperature < T0-40.d0 .or. temperature > T0+40.d0) then
         print*,'Warning, in FN diff_wat_vap, the parameterisation'
         print*,'  is valid between -40 and +40 Celsius'
         print*,'  The temperature is: ',temperature,' K.'
         print*,'  The parameterisation has nonetheless been used'
         print*,'  The calculated diffusivity is: ',diff_wat_vap,' m2/s'
      end if

      end function diff_wat_vap
      
!
!-------------------------------------------------------------
!

      function therm_conduct_air(temperature)

      ! Thermal conductivity of air

      ! ka = 1.d-3 * (4.39 + 0.071*T)

      ! where T is in K, and ka is in J/(m.s.K)
      ! see Seinfeld and Pandis 2nd Ed., equation (17.71) p.786


! Author
! ------
!     Josue Bock

! History
! -------
!     13-01-2017

      implicit none

      double precision :: therm_conduct_air

      double precision, intent(in) :: temperature  ! in [K]
      
      double precision, parameter :: cst1 = 4.39d-3
      double precision, parameter :: cst2 = 7.1d-5

      therm_conduct_air = cst1 + cst2*temperature


      end function therm_conduct_air
      
!
!---------------------------------------------------------------------
!

      subroutine fbil (ij)
! [apparently checking whether particle number is ever negative]
! currently never called

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      do k=1,nf+5
         do ia=1,nka
            do jt=1,nkt
               if (ff(jt,ia,k).lt.-0.1) print *,ij,ff(jt,ia,k),k,ia,jt
            enddo
         enddo
      enddo

      end subroutine fbil

!
!-------------------------------------------------------------
!

      subroutine advec_old (dt)
! one of many implementations of Bott's advection scheme

      USE global_params, ONLY :
! Imported Parameters:
     &     nkt

       implicit double precision (a-h,o-z)

      parameter (ymin=1.d-32)  !1.d-08
      common /cb59/ y(nkt),u(nkt)
      dimension z(nkt)
      do 1000 i0=1,nkt
         z(i0)=y(i0)
         y(i0)=0.0
 1000    if (z(i0).ge.ymin) go to 2000
! return if not enough particles:
      return
! find lowest (=i0) and highest (=i1) bin where f(i)=y(i) > ymin:
 2000 do 1010 i1=nkt,i0+1,-1
         z(i1)=y(i1)
         y(i1)=0.
 1010    if (z(i1).ge.ymin) go to 2010
      go to 2020
 2010 if (i1-i0.le.1) go to 2020
      do 1020 i=i0+1,i1-1
         z(i)=y(i)
 1020    y(i)=0.
! do growth/advection calculation only in between smallest and highest bin with sufficient particles
 2020 do 1030 i=i0,i1
         if (z(i).lt.ymin) go to 1030
         k2=0
         dt1=dt
         k=i
 3000    dt0=dmin1(1.d0/(dabs(u(k))+1.d-15),dt1)
         k1=k
         x0=float(k)+u(k)*dt0
         dt1=dt1-dt0
         if (dt1.le.1.d-07) go to 2030
         k=k+1
! avoid index out of bounds (1):
         if (k.gt.nkt) then
            k=nkt
            print *,'SR advec: index out of bounds'
         endif
         if (u(k).lt.0.) k=k-2  
! avoid index out of bounds (2):
         if (k.le.0) then
            k=1
            print *,'SR advec: index out of bounds'
         endif
         if (k.ne.k2) go to 2040
         y(k)=y(k)+z(i)
         go to 1030
 2040    k2=k1
         go to 3000
 2030    k=idint (x0+0.999999d0)
         c0=x0-float(k-1)
! "arithmetic if" to choose the correct polynomial
         if (i-2) 2050,2060,2070
 2070    if (nkt-1-i) 2050,2060,2080
! at first and last grid point (i=1,nkt) first order polynomial
! 2050 x0=c0*(z(i)+.5*(1.-c0)*(z(i+1)-z(i)))
 2050    x0=c0*z(i)
         go to 2090
! at second and second last grid point (i=2,nkt-1) second order polynomial
 2060    al=1.-2.*c0
         al2=al*al
         a0=(26.*z(i)-z(i+1)-z(i-1))/24.
         a1=(z(i+1)-z(i-1))/16.
         a2=(z(i+1)+z(i-1)-2.*z(i))/48.
         x0=dmin1(z(i),a0*c0+a1*(1.d0-al2)+a2*(1.d0-al2*al))
         go to 2090
! i=3,nkt-2: fourth order polynomial
 2080    al=1.-2.*c0
         al2=al*al
         al3=al2*al
         a0=(9.*(z(i+2)+z(i-2))-116.*(z(i+1)+z(i-1))+2134.*z(i))/1920.
         a1=(-5.*(z(i+2)-z(i-2))+34.*(z(i+1)-z(i-1)))/384.
         a2=(-z(i+2)+12.*(z(i+1)+z(i-1))-22.*z(i)-z(i-2))/384.
         a3=(z(i+2)-2.*(z(i+1)-z(i-1))-z(i-2))/768.
         a4=(z(i+2)-4.*(z(i+1)+z(i-1))+6.*z(i)+z(i-2))/3840.
         x0=dmin1(z(i),a0*c0+a1*(1.d0-al2)+a2*(1.d0-al3)
     &        +a3*(1.d0-al2*al2)+a4*(1.d0-al2*al3))
 2090    x0=dmax1(0.d0,x0)
         y(k)=y(k)+x0
         k=max0(k-1,1)
         y(k)=y(k)+z(i)-x0
 1030 continue

      end subroutine advec_old
!
!-------------------------------------------------------------
!

      subroutine advec (dt,u,y)

! advection of particles for condensational growth

! one of many implementations of Bott's advection scheme
! specific for 2D particle spectrum

! jjb work done
!    - removed arithmetic if for polynomial order (1, 2 or 4) at the end
!    - rewritten using up to date features: do, do while, exit and cycle
!    - also rewritten the initial search for min/max indexes by reading/writting arrays in increasing indexes, for computing efficiency
!    - removed archaic forms of Fortran intrinsic functions
!    - passed y and u as arguments instead of common block

      USE global_params, ONLY :
! Imported Parameters:
     &     nkt

      implicit none

      double precision, intent(in) :: dt
      double precision, intent(in) :: u(nkt)
      double precision, intent(inout) :: y(nkt)

      double precision, parameter :: ymin=1.d-32  !1.d-08

      integer :: i0, i1    ! lowest and highest bin indexes where y(i) >= ymin
      integer :: i, jt     ! loop indexes
      integer :: k, k1, k2
      integer :: k_low, k_high

      double precision :: a0, a1, a2, a3, a4
      double precision :: al, al2, al3
      double precision :: c0
      double precision :: dt0, dt1
      double precision :: x0
      double precision :: z(nkt)


! find lowest (=i0) and highest (=i1) bin where y(i) > ymin:
      i0=0 ! initialise i0 to track the case where no bin has y(i) > ymin
      do jt=1,nkt
         z(jt)=y(jt)
         y(jt)=0.d0
         if (z(jt)>=ymin) then
            i0=jt
            exit
         end if
      end do

      ! return if not enough particles in all nkt classes:
      if (i0==0) return

      i1=i0
      do jt=i0+1,nkt        ! note that if i0=nkt, this loop won't run, on purpose
         z(jt)=y(jt)
         y(jt)=0.d0
         if (z(jt)>=ymin) then
            i1=jt
         end if
      end do


! do growth/advection calculation only in between smallest and highest bin with sufficient particles
      iloop: do i=i0,i1

         ! not enough particles in the current bin
         if (z(i).lt.ymin) cycle

         k2=0
         dt1=dt
         k=i
         do
            k1=k

            if(abs(u(k))>0.d0) then
               dt0=min(1.d0/(abs(u(k))),dt1)
            else ! u(k) = 0.d0
               y(k)=y(k)+z(i)
               cycle iloop
            end if

            x0=float(k)+u(k)*dt0

            dt1=dt1-dt0
            if (dt1.le.1.d-07) exit

! update k:
            ! the new k can be:
            !   k_new = k_old - 1   if u(k_old) <0,
            !   k_new = k_old       if u(k_old) is small,
            !   k_new = k_old + 1   if u(k_old) is large.
            ! but practically, the case k_new = k_old has already exited the loop (dt1 case above)
            k=int(floor(x0))

            if (k==k2) then ! avoid inifinite loop (would happen if two adjacent class have opposite sign for u)
               y(k)=y(k)+z(i)
               cycle iloop ! cycle external i loop
            end if

            if (k.lt.1) then
               k=1
               print *,'SR advec: index out of bounds (2)',i
            else if (k.gt.nkt) then
               k=nkt
               print *,'SR advec: index out of bounds (1)',i
            endif

            k2=k1
         end do

         k_low = int(floor(x0))
         k_high = k_low + 1

         c0 = x0 - float(k_low) ! = x0 - floor(x0)

         if (c0 < 0.) then
            print*,'c0 is negative',i,k_low,k_high,dt0,dt1,x0,c0
            stop 'SR advec: error with c0'
         else if (c0 == 0.d0) then
            ! This case will happen in two situations:
            !   - if u(k)*dt0 is exactly = +/- 1.000d0 (this really happen sometimes)
            !   - if u(k)*dt0 is so small that it has been numerically rounded to 0.
            !     (very unlikely, would imply that u(k) ~ tiny(0.d0) and dt1 ~ [1.1e-7 -- 1.e-1] )
            y(k)=y(k)+z(i)
            cycle iloop
         end if
         if (k_high .gt. nkt .or. k_low .lt. 1) then
            print*,k_high,k_low,x0,c0,i
            stop 'SR advec: error with k_high or k_low'
         end if

         if (i==1 .or. i==nkt) then
            ! at first and last grid point (i=1,nkt): first order polynomial
            x0=c0*z(i)
         else if (i==2 .or. i==nkt-1) then
            ! at second and second last grid point (i=2,nkt-1): second order polynomial
            al=1.-2.*c0
            al2=al*al
            a0=(26.*z(i)-z(i+1)-z(i-1))/24.
            a1=(z(i+1)-z(i-1))/16.
            a2=(z(i+1)+z(i-1)-2.*z(i))/48.
            x0=min(z(i),a0*c0+a1*(1.d0-al2)+a2*(1.d0-al2*al))
         else
            ! i=3,nkt-2: fourth order polynomial
            al=1.-2.*c0
            al2=al*al
            al3=al2*al
            a0=(9.*(z(i+2)+z(i-2))-116.*(z(i+1)+z(i-1))+2134.*z(i))
     &         /1920.
            a1=(-5.*(z(i+2)-z(i-2))+34.*(z(i+1)-z(i-1)))/384.
            a2=(-z(i+2)+12.*(z(i+1)+z(i-1))-22.*z(i)-z(i-2))/384.
            a3=(z(i+2)-2.*(z(i+1)-z(i-1))-z(i-2))/768.
            a4=(z(i+2)-4.*(z(i+1)+z(i-1))+6.*z(i)+z(i-2))/3840.
            x0=min(z(i),a0*c0+a1*(1.d0-al2)+a2*(1.d0-al3)
     &           +a3*(1.d0-al2*al2)+a4*(1.d0-al2*al3))
         end if

         x0=max(0.d0,x0)

         y(k_low)  = y(k_low) + z(i) - x0
         y(k_high) = y(k_high)+ x0

      end do iloop

      end subroutine advec
!
!-------------------------------------------------------------
!

      subroutine advsed0
!     one of many implementations of Bott's advection scheme

! advection for sedimentation, upstream procedure
!
!     jjb removed declaration of 6 unused variables
! changed internal parameter n=nf, misleading since in most places n=nf+50
! (obviously, led to a bug when including erroneously "n" from global_params

! jjb cleaning of unused varaibles 15/12/16
! jjb declaration of all variables, implicit none 15/12/16
! jjb fortran generic functions min and max 15/12/16
! jjb checked identical to AB str code

      USE global_params, ONLY :
! Imported Parameters:
     &     nf

      implicit none

      integer :: i
      double precision :: fm(nf),fp(nf)
      
      common /cb58/ c(nf),y(nf)
      double precision c, y

! upstream procedure
      do i=1,nf-1
         fm(i)=-min(0.d0,c(i))*y(i+1)
         fp(i)=max(0.d0,c(i))*y(i)
      end do
      do i=2,nf-1
         y(i)=y(i)-fm(i-1)+fp(i-1)+fm(i)-fp(i)
      end do

      end subroutine advsed0

!
!-------------------------------------------------------------
!

      subroutine advsed1
! area preserving flux form; Bott (1989): Monthly Weather Review.
! fourth order monotone version.
! y(i) is transport quantity, input and output.
! boundary conditions are y(1)=const, y(n)=const.
! c(i) is Courant number satisfying the CFL criterion, input.
! fm(i), fp(i) are fluxes for u(i)<0 and u(i)>0, respectively.
! a0, a1, a2, a3, a4 are coefficients of polynomials in gridbox i.
! At i=1 and i=n first order polynomial,
! at i=2 and i=n-1 second order polynomial,
! at 3<=i<=n-2 fourth order polynomial.
! w(i) are weighting factors.
! the numerical grid is equidistant.
! the procedure is one dimensional.
! for multidimensional applications time splitting has to be used.
! the quantities c(i), fm(i), fp(i)  are given at the right
! boundary of grid cell i.
! Thus, fm(i) is flux from gridbox i+1 into gridbox i for c(i)<0,
! fp(i) is flux from gridbox i into gridbox i+1 for c(i)>0.

      USE global_params, ONLY :
! Imported Parameters:
     &     nf

      implicit double precision (a-h,o-z)

      common /cb58/ c(nf),y(nf)
      dimension a0(nf),a1(nf),a2(nf),a3(nf),a4(nf),fm(nf)
      a0(2)=(26.*y(2)-y(3)-y(1))/24.
      a1(2)=(y(3)-y(1))/16.
      a2(2)=(y(3)+y(1)-2.*y(2))/48.
      a3(2)=0.
      a4(2)=0.
      do i=3,nf-2
         a0(i)=(9.*(y(i+2)+y(i-2))-116.*(y(i+1)+y(i-1))
     &         +2134.*y(i))/1920.
         a1(i)=(-5.*(y(i+2)-y(i-2))+34.*(y(i+1)-y(i-1)))/384.
         a2(i)=(-y(i+2)+12.*(y(i+1)+y(i-1))-22.*y(i)-y(i-2))/384.
         a3(i)=(y(i+2)-2.*(y(i+1)-y(i-1))-y(i-2))/768.
         a4(i)=(y(i+2)-4.*(y(i+1)+y(i-1))+6.*y(i)+y(i-2))/3840.
      enddo
      a0(nf-1)=(26.*y(nf-1)-y(nf)-y(nf-2))/24.
      a1(nf-1)=(y(nf)-y(nf-2))/16.
      a2(nf-1)=(y(nf)+y(nf-2)-2.*y(nf-1))/48.
      a3(nf-1)=0.
      a4(nf-1)=0.
      cl=-c(nf-1)
      fm(nf-1)=dmin1(y(nf),cl*(y(nf)-(1.-cl)*(y(nf)-y(nf-1))*0.5))
      clm=cl
      do i=nf-1,2,-1
         cl=clm
         clm=-c(i-1)
         x1=1.-2.*cl
         x2=x1*x1
         x3=x1*x2
         ymin=dmin1(y(i),y(i+1))
         ymax=dmax1(y(i),y(i+1))
         fmim=dmax1(0.d0,a0(i)*cl-a1(i)*(1.-x2)+a2(i)*(1.-x3)
     &        -a3(i)*(1.-x1*x3)+a4(i)*(1.-x2*x3))
         fmim=dmin1(fmim,y(i)-ymin+fm(i))
         fmim=dmax1(fmim,y(i)-ymax+fm(i))
         fmim=dmax1(0.d0,fmim-(cl-clm)*y(i))
         w=y(i)/dmax1(fmim+1.d-15,y(i))
         fm(i-1)=fmim*w
      enddo
      y(1)=y(1)+fm(1)
      do i=2,nf-1
         y(i)=y(i)-fm(i-1)+fm(i)
      enddo
      y(nf)=y(nf)-fm(nf-1)

      end subroutine advsed1

!
!-------------------------------------------------------------
!

      subroutine advseda
! one of many implementations of Bott's advection scheme

! jjb arithmetic do replaced by do / enddo 15/12/2016

      USE global_params, ONLY :
! Imported Parameters:
     &     nf

      implicit double precision (a-h,o-z)

      common /cb58/ c(nf),y(nf)
      dimension flux(nf)
      do i=1,nf
         flux(i)=0.
      end do
      a0=(26.*y(2)-y(3)-y(1))/24.
      a1=(y(3)-y(1))/16.
      a2=(y(3)+y(1)-2.*y(2))/48.
      cl=-c(1)
      x1=1.-2.*cl
      x2=x1*x1
      fm=dmax1(0.d0,a0*cl-a1*(1.d0-x2)+a2*(1.d0-x1*x2))
      w=y(2)/dmax1(fm+1.d-15,a0+2.*a2)
      flux(1)=fm*w
      do i=3,nf-2
         a0=(9.*(y(i+2)+y(i-2))-116.*(y(i+1)+y(i-1))+2134.*y(i))/1920.
         a1=(-5.*(y(i+2)-y(i-2))+34.*(y(i+1)-y(i-1)))/384.
         a2=(-y(i+2)+12.*(y(i+1)+y(i-1))-22.*y(i)-y(i-2))/384.
         a3=(y(i+2)-2.*(y(i+1)-y(i-1))-y(i-2))/768.
         a4=(y(i+2)-4.*(y(i+1)+y(i-1))+6.*y(i)+y(i-2))/3840.
         cl=-c(i-1)
         x1=1.-2.*cl
         x2=x1*x1
         x3=x1*x2
         fm=dmax1(0.d0,a0*cl-a1*(1.-x2)+a2*(1.-x3)-a3*(1.-x1*x3)
     &        +a4*(1.-x2*x3))
         w=y(i)/dmax1(fm+1.d-15,a0+2.d0*(a2+a4))
         flux(i-1)=fm*w
      end do
      a0=(26.*y(nf-1)-y(nf)-y(nf-2))/24.
      a1=(y(nf)-y(nf-2))/16.
      a2=(y(nf)+y(nf-2)-2.*y(nf-1))/48.
      cl=-c(nf-2)
      x1=1.-2.*cl
      x2=x1*x1
      fm=dmax1(0.d0,a0*cl-a1*(1.d0-x2)+a2*(1.d0-x1*x2))
      w=y(nf-1)/dmax1(fm+1.d-15,a0+2.d0*a2)
      flux(nf-2)=fm*w
      cl=-c(nf-1)
      flux(nf-1)=dmin1(y(nf),cl*(y(nf)-(1.d0-cl)*(y(nf)-y(nf-1))*0.5d0))
      y(1)=y(1)+flux(1)
      do i=2,nf-1
         y(i)=y(i)-flux(i-1)+flux(i)
      end do
      y(nf)=y(nf)-flux(nf-1)

      end subroutine advseda


!
!------------------------------------------------------------------------
!

!     subroutine stem_kpp (dd,xra,z_box,n_bl,box)     ! jjb
      subroutine stem_kpp (dd,xra,z_box,n_bl,box,nuc) ! jjb nuc is needed in 2 IF tests

      USE constants, ONLY :
! Imported Parameters:
     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc

      implicit double precision (a-h,o-z)
! chemical reactions
! aerosol mass change due to chemical reactions
!
! NOTE: the aerosol processing as implemented here is mass conserving
!       but shifts too many particles into the smallest bins, therefore
!       the particle spectra are a bit odd
!       this has to be improved in future versions

      logical box,nuc ! jjb nuc has to be declared after adding it as an argument
      integer tix,tixp

      parameter (lsp=9)
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /blck06/ kw(nka),ka
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /liq_pl/ nkc_l
      ! Thus commented here, and the lines "feeding" it commented as well to save cpu time
      ! If used again, the indexes should be re-ordered to save cpu-time (n=last, lsp first for dss)

      dimension lj2(lsp)
      dimension sion1o(lsp,nkc,n) ! sion1o old ion conc. [mole m**-3]
      dimension dsion1(lsp,nkc,n) ! dsion1 change in ion conc. [mole/part.]
      ! jjb fs is not used over all its dimensions (only in an old line commented)
      ! maybe worth to reduce its dimensions after double check, to save cpu time
      dimension fs(nka,nkc,n)
      dimension sap(nkc,n) ! sap total number of aerosols [cm**-3]
      dimension smp(nkc,n) ! smp total aerosol mass [mg cm**-3]
      dimension vc(nkc,nkc,n)
      data lj2/1,2,8,9,13,14,19,20,30/
      common /nucfeed/ ifeed ! jjb added so that ifeed is known in this SR (used in 2 IF tests)

!      fpi=4./3.*3.1415927
      fpi=4./3.*pi
      
      nmin = 2
      nmax = nf
      if (box) then
         nmin = n_bl
         nmax = n_bl
      endif
!      print*,'in stem_kpp call aer_source'
      call aer_source (box,dd,z_box,n_bl)
!      print*,'in stem_kpp call liq_parm'
      call liq_parm (xra,box,n_bl)

!*************************** no aerosol processing *****************
!      call kpp_driver (box,dd,n_bl)
!      return
!*************************** no aerosol processing *****************

! sion1 ion conc. [mole m**-3], sion1o old ion conc. [mole m**-3], 
! dsion1 change in ion conc. [mole/part.]
! sap total number of aerosols [cm**-3], smp total aerosol mass [mg cm**-3]
! fs total aerosol mass integrated over all liquid classes [mg cm**-3]
! lj2 species that define the aerosol mass

!      if (.not.box) then

! loop over aqueous chemistry layers to get values of sion1, fs, smp, sap
! before chemistry integration (at t=t_0)
      do k=nmin,nmax
         do kc=1,nkc_l          !initialize variables
            do ic=1,nkc_l
               vc(ic,kc,k)=0.
            enddo
            sap(kc,k)=0.
            smp(kc,k)=0.
            if (cm(kc,k).eq.0.) goto 900
!           define upper and lower limits of ia loop
            if (kc.eq.1.or.kc.eq.3) then
!               if ((nuc).and.(ifeed.eq.1)) then ! jjb corrected (?) below CHECH WITH Susanne PECHTL
!                 ial=1
!               else
!                 ial=1+1
!               endif
               if ((nuc).and.(ifeed.eq.2)) then
                 ial=2
               else
                 ial=1
               endif
               iau=ka
            else
               ial=ka+1
               iau=nka-1
            endif
!           loop over all ia that are in the current kc bin
            do ia=ial,iau
               fs(ia,kc,k)=0.
!              define upper and lower limits of jt loop
               if (kc.eq.1.or.kc.eq.2) then
                  jtl=1
                  jtu=kw(ia)
               else
                  jtl=kw(ia)+1
                  jtu=nkt
               endif
!              loop over all jt that are in the current kc bin
               do jt=jtl,jtu
                  fs(ia,kc,k)=fs(ia,kc,k)+ff(jt,ia,k)*en(ia)
                  sap(kc,k)=sap(kc,k)+ff(jt,ia,k)
               enddo
               smp(kc,k)=smp(kc,k)+fs(ia,kc,k)
            enddo
            do l=1,lsp
               ll=lj2(l)
               sion1o(l,kc,k)=sion1(ll,kc,k)
            enddo
 900        continue
         enddo   ! kc
      enddo      ! k

!      endif ! .not.box
! chemistry SR
!      print*,'in stem_kpp call kpp_driver'
      call kpp_driver (box,dd,n_bl)
!      if (box) return

! redistribution of particles along aerosol grid due to modified aerosol mass
      do k=nmin,nmax
         do kc=1,nkc_l
            if (cm(kc,k).eq.0.) goto 1000
            if (sap(kc,k).gt.1.e-6) then
!              change in mass determining chemical species
               do l=1,lsp
                  ll=lj2(l)
                  dsion1(l,kc,k)=(sion1(ll,kc,k)-sion1o(l,kc,k))*1.e-06
     &                 /sap(kc,k) 
!                  dss(k,l,kc)=dss(k,l,kc)+dsion1(l,kc,k) ! jjb output no longer used
               enddo
! den: new aerosol mass in mg/particle due to chemical reactions
! mole masses for l=1,2,8,9,13,14,19,20,30: H+=1g/mole, NH4=18g/mole, SO4(2-)=96g/mole,
! HCO3-=61g/mole, NO3=62g/mole, Cl-=35.5g/mole, HSO4-=97g/mole, Na+=23g/mole, CH3SO3-=95g/mole
! HCO3-=61g/mole --> 44 g/mole as water remains in particle when CO2 degasses due to acidification
               den=(dsion1(1,kc,k)*1.+dsion1(2,kc,k)*18.+dsion1(3,kc,k)
!     &              *96.+dsion1(4,kc,k)*61.+dsion1(5,kc,k)*62.+
     &              *96.+dsion1(4,kc,k)*44.+dsion1(5,kc,k)*62.+
     &              dsion1(6,kc,k)*35.5+dsion1(7,kc,k)*97.+
     &              dsion1(8,kc,k)*23.+dsion1(9,kc,k)*95.)*1000.
               
!              define upper and lower limits of ia loop
               if (kc.eq.1.or.kc.eq.3) then
                  if ((nuc).and.(ifeed.eq.1)) then
                    ial=1
                  else
                    ial=1+1
                  endif
                  iau=ka
               else
                  ial=ka+1
                  iau=nka-1
               endif
!              loop over all ia that are in the current kc bin
!              c0: courant number for redistribution
! if growing reverse loop order to avoid increasing the mass of some
! particles twice
               istart=ial
               iend=iau
               iinkr=1
               if (den.ge.0.) then
                  istart=iau
                  iend=ial
                  iinkr=-1
               endif
!               do ia=ial,iau
               do ia=istart,iend,iinkr
                  if (den.gt.0.) then
                     x0=en(ia)+den*en(ia)/smp(kc,k)*sap(kc,k)
                  else
!                     x0=en(ia)+den*fs(ia,kc,k)/smp(kc,k)*sap(kc,k)
                     x0=en(ia)+den*en(ia)/smp(kc,k)*sap(kc,k)
!                     x0=dmax1(en(ia)+den,0.d0)
                     if (x0.le.0.d0) print *,k,kc,ia,'aerosol growth'
                  endif
                  do iia=1,nka-1
                     if (en(iia).le.x0.and.en(iia+1).gt.x0) then
                        ix=iia
                        c0=(en(iia+1)-x0)/(en(iia+1)-en(iia))
!                        c0=(dlog10(en(iia+1))-dlog10(x0))/dlgenw
                        go to 2000
                     endif
                  enddo
                  if (en(1).gt.x0) then
                     ix=1
                     c0=1.
                  else
                     ix=nka-1
                     c0=0.
                  endif
 2000             continue 
!              define upper and lower limits of jt loop               
               if (kc.eq.1.or.kc.eq.2) then
                  jtl=1
                  jtu=kw(ia)
               else
                  jtl=kw(ia)+1
                  jtu=nkt
               endif
!              loop over all jt that are in the current kc bin
                  do jt=jtl,jtu
                     if (ff(jt,ia,k).gt.0.) then
                        x1=ff(jt,ia,k)
                        ff(jt,ia,k)=0.
                        ff(jt,ix,k)=ff(jt,ix,k)+x1*c0  
                        ff(jt,ix+1,k)=ff(jt,ix+1,k)+x1*(1.-c0) 
! find "targetbin" for ix and ix+1:
! ix
                        if (ix.gt.ka) then
                           if (jt.gt.kw(ix)) then
                              tix=4
                           else
                              tix=2
                           endif
                        else
                           if (jt.gt.kw(ix)) then
                              tix=3
                           else
                              tix=1
                           endif
                        endif
! ix+1
                        if (ix+1.gt.ka) then
                           if (jt.gt.kw(ix+1)) then
                              tixp=4
                           else
                              tixp=2
                           endif
                        else
                           if (jt.gt.kw(ix+1)) then
                              tixp=3
                           else
                              tixp=1
                           endif
                        endif
! store change of volume for ix --> ix and ia --> ix+1
                        if (tix.ne.kc)  vc(tix,kc,k) =vc(tix,kc,k) +x1*
     &                       c0*fpi*rq(jt,ia)**3
                        if (tixp.ne.kc) vc(tixp,kc,k)=vc(tixp,kc,k)+x1*
     &                       (1.-c0)*fpi*rq(jt,ia)**3
! output for control
!                       if (tix.ne.kc)  fss(k,kc,tix)=fss(k,kc,tix) +x1* ! jjb output no longer used
!     &                       c0
!                       if (tixp.ne.kc) fss(k,kc,tixp)=fss(k,kc,tixp)+x1* ! jjb output no longer used
!     &                       (1.-c0)
                     endif
                  enddo   ! jt loop
               enddo      ! ia loop
            endif
            
 1000       continue   
         enddo  ! kc loop
      enddo     ! k  loop
! move chemical species if transport above chemistry bins took place
      do k=nmin,nmax
         do kc=1,nkc_l
            do kkc=1,nkc_l
               if (kkc.eq.kc) goto 1001
! "kc" (from) bin and "kkc" (to) bin, i.e. bin that loses moles and bin that gains moles
               if (vc(kkc,kc,k).eq.0.) goto 1001
               if (cw(kc,k).gt.0.) then
                  nfrom=kc
                  nto=kkc
! cw in m^3(aq)/m^3(air), vc in um^3/cm^3: 10^-12
                  vol_ch=vc(kkc,kc,k)*1.d-12
                  xfact=0.
!            if (vol_ch.gt.1.e-4*cw(k,nfrom)) then
                  xfact=1.-(cw(nfrom,k)-vol_ch)/cw(nfrom,k)
                  do l=1,j2
                     xch=sl1(l,nfrom,k)*xfact
                     sl1(l,nfrom,k)=sl1(l,nfrom,k)-xch
                     sl1(l,nto,k)  =sl1(l,nto,k) +xch
                  enddo
                  do l=1,j6
                     xch=sion1(l,nfrom,k)*xfact
                     sion1(l,nfrom,k)=sion1(l,nfrom,k)-xch
                     sion1(l,nto,k)  =sion1(l,nto,k) +xch
                  enddo
               end if
!            endif  ! if vol_ch > 10^-4* cw
 1001          continue
!         enddo  ! ic loop
            enddo               ! kkc loop
! 1002       continue ! jjb statement label unreferenced
         enddo                  !kc loop
      enddo     ! k loop

!      do k=2,n
!         do kc=1,nkc_l
!            do kcc=1,nkc_l
!               svc(k,kc,kcc)=svc(k,kc,kcc)+vc(kcc,kc,k) ! jjb output no longer used
!            enddo
!         enddo
!      enddo

      end subroutine stem_kpp

!
!----------------------------------------------------------------------
!

      subroutine adjust_f
! adjustment of initial aerosol size distribution for specific scenarios

!      USE constants, ONLY :
!! Imported Parameters:
!     & pi

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      dimension f_inter(nka)

      x0=1.

      do k=2,n
         do ia=1,nka
            f_inter(ia)=ff(1,ia,k)
         enddo

         fsum(k)=0.
         do ia=nka,1,-1
            if (rn(ia).gt.0.5) x0=0.1
! init spectrum:aerosols in equilibrium with rH=76.2 %
! "dry" them !(this is not really exact..)
! equilibrium radius at 76 %
            a0=a0m/t(k)
            b0=b0m(ia)*2.
! b0=b0m*rho3/rhow; rho3=2000; rhow=1000
            rg=rgl(rn(ia),a0,b0,.762d0)
!            eg=4.d-09*pi/3.*(rg**3-rn(ia)**3)
            do kl=1,nka
               if (rg.le.rn(kl)) then
!               if (eg.le.enw(kl)) then   
                  ff(1,ia,k)=x0*f_inter(kl)
                  goto 1000
               endif
            enddo
 1000       continue
            fsum(k)=fsum(k)+ff(1,ia,k)
         enddo

      enddo

      end subroutine adjust_f

!----------------------------------------------------------------------

      subroutine partdep (ra)
! calculate particle dry deposition velocity after Seinfeld and Pandis,
! 1998, p.958ff

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

      implicit double precision (a-h,o-z)

      common /blck06/ kw(nka),ka
      common /blck12/ cw(nkc,n),cm(nkc,n)
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb46/ ustern,gclu,gclt
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kinv_i/ kinv
      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
!      double precision ra,rb,vs,vd,z,xD,sc,st,rx,Cc
      dimension xx1(nkc)

      integer ia,jt
      call monin (phi)

! particle dry deposition velocity:v_d=1/(ra + rb + ra rb v_s)+ v_s
!    ra=1/(kappa ustar) (ln (z/z0) +Phi) ;where z=height of surface (constant flux) layer
!                                         Phi takes stratification into account
!    rb=1/(ustar(Sc^(2/3) + 10^(-3/St))) ;Sc=nu/D  St=v_s*ustar^2/(g nu)

      xk=1.38066d-23 !Boltzmann number
      z=0.1*eta(kinv) !surface layer height: 10% of BL (Stull), insensitive parameter
      ra=1./(0.4*ustern)*(dlog(z/z0)+phi)  !ra:aerodynamic resistance; kappa=0.4

      k=2 ! only in lowest model layer
      xeta=1.8325e-5*(416.16/(t(k)+120.))*((t(k)/296.16)**1.5) !dynamic viscosity of air, Jacobsen p. 92
      xnu=xeta/rho(k)  !kinematic viscosity of air
!      write (110,*) 'nu eta',xnu,xeta

! set xx1(kc)=0.
      do kc=1,nkc
         xx1(kc) = 0.d0
         vdm(kc) = 0.d0 ! jjb added this initialisation. If not, old values are still used when cw goes to 0 from t to t+1
      enddo
! free path length
      xlam = 2.28e-5 * t(k) / p(k)

      do ia=1,nka
         do jt=1,nkt
            rx=rq(jt,ia)*1.e-6
            vs=vterm(rx,t(k),p(k)) ! Stokes fall velocity incl Cc
            Cc=1.+xlam/rx*(1.257+.4*dexp(-1.1*rx/xlam)) ! Cunningham slip flow corr.
            xD=xk*t(k)*Cc/(6*pi*xeta*rx) ! aerosol diffusivity
            Sc=xnu/xD !Schmidt number
            St=vs*ustern**2/(g*xnu) !Stokes number
            rb=1./(ustern*(sc**(-2./3.)+10**(-3./st))) !quasi laminar resistance
            vd(jt,ia)=1./(ra+rb+ra*rb*vs)+vs !deposition velocity
!            write (110,20) ia,jt,rq(jt,ia),vs*100.,vd*100.,100./
!     &           (ra+rb+ra*rb*vs)
! calculate mass weighted mean dry deposition velocities
! loop over the nkc different chemical bins
            do kc=1,nkc
               if (cw(kc,k).eq.0.) goto 1001
! define kc limits  ---
               if (kc.eq.1.and.(ia.gt.ka.or.jt.gt.kw(ia))) goto 1001
               if (kc.eq.2.and.(ia.lt.(ka+1).or.jt.gt.kw(ia))) goto 1001
               if (kc.eq.3.and.(ia.gt.ka.or.jt.lt.(kw(ia)+1))) goto 1001
               if (kc.eq.4.and.(ia.lt.(ka+1).or.jt.lt.(kw(ia)+1))) goto 
     &              1001
! LWC weighted deposition velocity
               xx1(kc)=xx1(kc)+rx*rx*rx*vd(jt,ia)*ff(jt,ia,k)*1.e6
! deposition velocity:
!!              vdm(kc)=4.*3.1415927/(3.*cw(kc,k))*xx1(kc)
!              vdm(kc)=4.*pi/(3.*cw(kc,k))*xx1(kc)
 1001          continue
            enddo               ! kc
         enddo
      enddo

      do kc=1,nkc
         if (cw(kc,k).gt.0.d0)
!     &              vdm(kc)=4.*3.1415927/(3.*cw(kc,k))*xx1(kc) ! don't make this calculation nka * nkt * nkc times!
     &              vdm(kc)=4.*pi/(3.*cw(kc,k))*xx1(kc) ! don't make this calculation nka * nkt * nkc times!
      end do

!      do kc=1,nkc
!         rx=rc(2,kc)
!         if (rx.gt.0.) vs=vterm(rx,t(2),p(2))
!         write (110,21) kc,vdm(kc),vs,rc(2,kc)
!      enddo
! 20   format (1p,2i3,5d16.8)
! 21   format (1p,i3,3d16.8)

      end subroutine partdep

!
!-------------------------------------------------------
!

      subroutine monin (phi)
! calculate the Monin-Obukhov length after Seinfeld and Pandis, 1998, p.862

! jjb work done
!     - integer exponents (**2 instead of **2.)
!     - removed archaic forms of intrinsic functions (dlog, datan)

      USE constants, ONLY :
! Imported Parameters:
     &     cp              ! Specific heat of dry air, in J/(kg.K)

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka

      implicit double precision (a-h,o-z)

      double precision, intent(out) :: phi ! see S & P 1st Ed, p. 963, equation (19.14)
      double precision :: xeta

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb42/ atke(n),atkh(n),atkm(n),tke(n),tkep(n),buoy(n)
      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb46/ ustern,gclu,gclt
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /kinv_i/ kinv

! check inversion height - it is diagnosed in SR atk1 but might be zero after restart
      if (kinv.eq.0)  then
         kinv = 70
         print *,'SR monin: kinv = 0., set to kinv=70'
      endif

! reference height = 10 % of BL height
      z=0.1*eta(kinv)
      do k=1,kinv
         if (eta(k).ge.z) goto 100
      enddo
 100  continue

! L=-rho c_p T_0 Ustar^3/(kappa g \bar(q_3))
! with \bar(q_3)=rho c_p \bar(w'theta')=rho c_p (-1.) atkh d theta/d z
      
! dtheta'/dz in height k
      if (k.eq.1) then
         k=2
         print *,'SR monin: index out of bounds (1)'
      endif
      if (k.eq.n) then
         k=n-1
         print *,'SR monin: index out of bounds (2)'
      endif
      dtdz=((theta(k+1)-theta(k))/deta(k)+(theta(k)-theta(k-1))/
     &     deta(k-1))/2.
      q3=rho(k)*cp*(-1.)*atkh(k)*dtdz

      xmo=-1.*rho(k)*cp*t(1)*ustern**3/(0.4*g*q3) ! Seinfeld 2, p. 747, (16.70)

! effect on ra
      zeta=z/xmo
      zeta0=z0/xmo
! |L|>10^5  (neutral) => phi=0.
      if (abs(xmo) > 1.d5) then
         phi=0.
      else
! stable
         if (xmo.gt.0.) then
            phi=4.7*(zeta-zeta0)
! unstable      
         else if (xmo.lt.0) then
            xeta0=(1.-15.*zeta0)**0.25
            xeta=(1.-15.*zeta)**0.25
            phi=log( (xeta0**2+1.)*(xeta0+1.)**2
     &               /((xeta**2+1.)*(xeta+1.)**2) ) 
     &       +2.*(atan(xeta)-atan(xeta0))
         else            ! jjb: note that the case xmo=0 is not explained in S & P
            print*,'Warning, in SR monin, xmo=0',xmo
         endif
      endif



      print *,'L, Ri',xmo,z/xmo
      print *,'Phi= ',phi

      end subroutine monin

!
!----------------------------------------------------------------------
!


      subroutine ion_mass (srname)
! calculation of ion mass for ion balance checks

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc

      implicit double precision (a-h,o-z)

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /liq_pl/ nkc_l
      character *10   srname
      dimension xsum(nf)

      xHp=0.
      xNHp=0.
      xSOm=0.
      xHCOm=0.
      xNOm=0.
      xClm=0.
      xHSOm=0.
      xNap=0.
      xCHSO=0.
      xxsum=0.
!      xsumi=0. !  redefined later before referenced; no initialisation needed
      do k=2,nf
         do j=1,nkc_l
            xHp=xHp+     sion1(1,j,k) *detw(k)*1.  *1.d6
            xNHp=xNHp+   sion1(2,j,k) *detw(k)*19. *1.d6
            xSOm=xSOm+   sion1(8,j,k) *detw(k)*96. *1.d6
            xHCOm=xHCOm+ sion1(9,j,k) *detw(k)*61. *1.d6
            xNOm=xNOm+   sion1(13,j,k)*detw(k)*62. *1.d6
            xClm=xClm+   sion1(14,j,k)*detw(k)*35.5*1.d6
            xHSOm=xHSOm+ sion1(19,j,k)*detw(k)*97. *1.d6
            xNap=xNap+   sion1(20,j,k)*detw(k)*23. *1.d6
            xCHSO=xCHSO+ sion1(30,j,k)*detw(k)*95. *1.d6
         enddo
         xsum(k)=0.
         do ia=1,nka
            do jt=1,nkt
               xsum(k)=xsum(k)+ff(jt,ia,k)*en(ia)
            enddo
         enddo
         xsum(k)=xsum(k)*1.e+09
         xxsum=xxsum+xsum(k)*detw(k)
      enddo

      xsumi=xHp+xNHp+xSOm+xHCOm+xNOm+xClm+xHSOm+xNap+xCHSO

      write (*,21) srname
      write (*,22)xxsum,xsumi,xNHp,xSOm,xHCOm,xNOm,xClm,xHSOm,xNap,xCHSO

 21   format (a10)
 22   format (10d14.6)

      end subroutine ion_mass


!
!-------------------------------------------------------
!

      subroutine box_init (nlevbox,nz_box,n_bl,BL_box)
!     initialisation for box models runs
      
      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n

      implicit double precision (a-h,o-z)
      logical BL_box

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
!      common /boxdat/ t0, xm10 ! this CB was fed here, but used nowhere else
      common /kinv_i/ kinv

! initialize kinv (needed in SR kpp_driver)
      kinv=nf
! important for restart only: call init_konc only to start with "fresh" aerosol
!      call init_konc

      if (.not.BL_box) then
! init vals from certain level:
         t0   = t(nlevbox)
         xm10 = xm1(nlevbox)
         print *,' level used for box run: ',nlevbox
      else
!  arithmetic average over box
         tsum  = 0.
         xmsum = 0.
         do k=2,nz_box
            tsum  = tsum + t(k)
            xmsum = xmsum + xm1(k)
         enddo
         t0   = tsum/(nz_box-1)
         xm10 = xmsum/(nz_box-1)
      endif
      t(n_bl)  = t0
      xm1(n_bl)= xm10
! p21 after magnus formula
      p21=610.7*dexp(17.15*(t(n_bl)-273.15)/(t(n_bl)-38.33))
      feu(n_bl)=xm1(n_bl)*p(n_bl)/((0.62198+0.37802*xm1(n_bl))*p21)
      print *,"box temp and hum: ",t(n_bl),xm1(n_bl),feu(n_bl)

! for strange reasons it didn't work to init the whole Mistra column and 
! not producing a crash or strange results, therefore this 1D column
! init is repeated every hour in SR box_update

! if smogchamber run: adjust roughness length z0
!     z0 = to be determined; maybe scale with the decline of the particle population
!     in the chamber; make sure that this value is not overwritten later during the 
!     run

! also make sure for smogchamber that deposition occurs also to the sidewalls of the chamber

      end subroutine box_init

!
!-------------------------------------------------------
!

!      subroutine box_update (box_switch,ij,nlevbox,nz_box,n_bl,chem, ! jjb unused arguments CHEM, HALO, IOD
      subroutine box_update (box_switch,ij,nlevbox,nz_box,n_bl,
!      &     halo,iod,BL_box)
     &     BL_box)

! update of astronomical, aerosol and meteorological properties for box model runs

      USE constants, ONLY :
! Imported Parameters:
     & pi
      
      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nrlay,
     &     nkc,
     &     mbs

      implicit double precision (a-h,o-z)

!      logical chem,halo,iod,fa_lse,BL_box ! jjb unused arguments removed: chem, halo, iod
      logical fa_lse,BL_box

      common /cb16/ u0,albedo(mbs),thk(nrlay)
      double precision u0, albedo, thk

      common /cb18/ alat,declin                ! for the SZA calculation
      double precision alat,declin

!      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb40/ xtime,lday,lst,lmin,it,lcl,lct
      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb54/ xm1(n),xm2(n),feu(n),dfddt(n),xm1a(n),xm2a(n)
      common /blck12/ cw(nkc,n),cm(nkc,n)

!     dimension rc(nf,nkc),freep(nf) ! jjb rc now in blck11
      dimension freep(nf)

!      nmin = n_bl ! jjb variable unreferenced
!      nmax = n_bl ! jjb variable unreferenced

      fa_lse = .false. ! jjb was missing, thus undefined when calling SRs fast_k_mt_* below

! calculate u0 
      rlat=alat*1.745329e-02
      rdec=declin*1.745329e-02
      zeit=lst*3600.+dfloat(lmin-1)*60.
! greater intervals and variable dtg is known only in the old chemical module
!      horang=7.272205e-05*zeit-3.1415927
      horang=7.272205e-05*zeit-pi
      u00=dcos(rdec)*dcos(rlat)*dcos(horang)+dsin(rdec)*dsin(rlat)
      ru0=6371.*u00
      u0=8./(sqrt(ru0**2+102000.)-ru0)

! new photolysis rates are calculated from main program! 

! get data for whole column for cw, rc, xkmt for box runs with averaging 
! over BL:

! if T, rH, .. vary in time update every hour or at "reasonable times" if not,
! just call this at the beginning of the run

!      if (lmin.eq.1.and.ij.eq.1) then
      if (box_switch.eq.1..and.lmin.eq.1.and.ij.eq.1) then
!        free path length (lambda=freep):
         do k=2,nf
            freep(k)=2.28e-5 * t(k) / p(k)
         enddo

         call cw_rc (nf) 
!!         do k=2,nf
!!            do kc=1,nkc ! jjb reordered
!         do kc=1,nkc
!            do k=2,nf
!               cwm(k,kc)=cw(k,kc)
!               rcm(k,kc)=rc(k,kc)
!            enddo
!         enddo
         stop 'jjb: box version has to be updated'
         call v_mean (t(:nf))
         !call v_mean_a  (t,nf)  
!        call henry_a (t,p,nf) ! jjb second argument (p) not used
         call henry_a (t,nf)   ! jjb removed
!        call fast_k_mt_a(freep,cw,fa_lse,nf) ! jjb cw now passed as a CB
         call fast_k_mt_a(freep,fa_lse,nf)
!        call equil_co_a (cw,t,nf)
         call equil_co_a (t,nf) ! jjb cw now passed as a CB
         call activ (fa_lse,nf)
         call dry_cw_rc (nf)
         call dry_rates_g (t,p,nf)
         call dry_rates_a (freep,nf)
         xph3=0.
         xph4=0.
         if (cm(3,n_bl).gt.0.) xph3 = 1.
         if (cm(4,n_bl).gt.0.) xph4 = 1.
         if (xph3.eq.1..or.xph4.eq.1.) then
            !call v_mean_t  (t,nf) 
!           call henry_t (t,p,nf) ! jjb second argument (p) not used
            call henry_t (t,nf)   ! jjb removed
!           call fast_k_mt_t(freep,cw,fa_lse,nf) ! jjb cw now passed as a CB
            call fast_k_mt_t(freep,fa_lse,nf)
!           call equil_co_t (cw,t,nf) ! jjb cw now passed as a CB
            call equil_co_t (t,nf)
         endif

         if (BL_box) then 
!           average parameters over depth of BL if BL_box=.true.
            call ave_parms (n_bl,nz_box) 
            call ave_aer (n_bl,nz_box) 
            if (xph3.eq.1..or.xph4.eq.1.) 
     &           call ave_tot (n_bl,nz_box) 
         else
            call set_box_gas (nlevbox,n_bl) 
            call set_box_lev_a (nlevbox,n_bl) 
! p21 after magnus formula
            p21=610.7*dexp(17.15*(t(n_bl)-273.15)/(t(n_bl)-38.33))
            feu(n_bl)=xm1(n_bl)*p(n_bl)/((0.62198+0.37802*xm1(n_bl))*
     &           p21)
            call equil (1,n_bl)
            if (xph3.eq.1..or.xph4.eq.1.) 
     &           call set_box_lev_t (nlevbox,n_bl) 
         endif
!         call print_vals (nlevbox,n_bl)
         box_switch = 0.
      endif


! the following parameters are used at the beginning of SR kpp_driver, so they
! have an effect on the chemistry

! prescribe a variation in T, rh as function of u0
!     put SINUS here, t0 xm10
!      t(n_bl)    = xsinus * t0        ! temp [K]
!      xm1(n_bl)  = xsinus * xm10      ! spec humidity [kg/kg]
!c p21 after magnus formula
!      p21=610.7*dexp(17.15*(t(n_bl)-273.15)/(t(n_bl)-38.33))
!      feu(n_bl)=xm1(n_bl)*p(n_bl)/((0.62198+0.37802*xm1(n_bl))*p21)

!----------tests from run bug_fix_n: str.f SR surf0 (see also SST.f in bug_fix_n)
!      common /tw_0/ tw0
!c      tw=tw-5.787d-6*dt
!      tw0=tw0-6.94444d-6*dt
!c assume diurnal variation in SST:
!      th      = mod(time/3600.,24.) ! time in [h] (0..24)
!!     pi05    = 0.5 * 3.1425927     ! 1/2 pi
!      pi05    = 0.5 * pi            ! 1/2 pi
!      tmax    = 21.                 ! shift to get maximum of SST at 15:00
!      sst_amp = 1.                  ! amplitude of SST diu var
!      sst_var = sst_amp * sin(pi05*(tmax - th)/6.)
!      tw      = tw0 + sst_var
!      t(1)=tw 
!----------
! change only for different LWC
!      conv2(n_bl,1-4)    ! conversion factor [1/(1000 m3/m3)]
!      cm(k,1-4)          ! LWC in particle class [m3/m3]
!      cloud(n_bl,1-4)    ! cloud present in layer [true/false] - should be false everywhere
! don't have to be changed
!      cm3(n_bl,1)        ! air dens. [mlc/cm3]
!      am3(n_bl,1)        ! air dens. [mol/m3]
!      am3(n_bl,2)        ! [CO] [mol/m3]
!      rho(n_bl)          ! density [kg/m3]
!      p(n_bl)            ! pressure [Pa]

      end subroutine box_update


!
!---------------------------------------------------------------
!

      subroutine sedc_box (dt,z_box,n_bl)
! dry deposition and emission of gaseous species for box runs

! jjb work done = implicit none, missing declarations, little cleaning, modules including constants

      USE constants, ONLY :
! Imported Parameters:
     & Avogadro

      USE gas_common, ONLY :
! Imported Parameters:
     &     j1,
! Imported Array Variables with intent (in):
     &     es1,
     &     ind_gas_rev,
! Imported Array Variables with intent (inout):
     &     s1,
     &     vg
      
      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision dt, z_box
      integer n_bl

! Local scalars:
      integer j
      double precision s12old

! Common blocks
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      double precision time
      integer lday, lst, lmin, it, lcl, lct

!- End of header ---------------------------------------------------------------

!      vg(4)=0.27e-2 ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or 
!      vg(34)=vg(30) ! N2O5=HCl
!      vg(37)=0.     ! DMS           emission is net flux 
!      vg(38)=vg(30) ! HOCl = HCl
!      vg(43)=vg(30) ! HOBr = HCl
!      vg(50)=vg(49) ! I2O2=HOI
!      vg(51)=vg(49) ! INO2=HOI
!      vg(56)=0.     ! CH3I          emission is net flux
!      vg(57)=0.     ! CH2I2         emission is net flux
!      vg(58)=0.     ! CH2ClI        emission is net flux
!      vg(59)=0.     ! C3H7I         emission is net flux
!      vg(63)=vg(30) ! CH3SO3H = HCl
!      vg(71)=0.     ! CH2BrI         emission is net flux
!      vg(72)=0.     ! CHBr2I         emission is net flux
!      vg(73)=0.     ! C2H5I          emission is net flux

      if(ind_gas_rev(4) /= 0)
     & vg(ind_gas_rev(4))=0.27e-2              ! NH3 old value, that fitted "nicely" in model    !=0. ! emission is net flux or 
      if(ind_gas_rev(34) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(34))=vg(ind_gas_rev(30)) ! N2O5=HCl
      if(ind_gas_rev(37) /= 0)
     &  vg(ind_gas_rev(37))=0.                  ! DMS           emission is net flux 
      if(ind_gas_rev(38) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(38))=vg(ind_gas_rev(30)) ! HOCl = HCl
      if(ind_gas_rev(43) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(43))=vg(ind_gas_rev(30)) ! HOBr = HCl
      if(ind_gas_rev(50) /= 0 .and. ind_gas_rev(49) /= 0)
     &  vg(ind_gas_rev(50))=vg(ind_gas_rev(49)) ! I2O2=HOI
      if(ind_gas_rev(51) /= 0 .and. ind_gas_rev(49) /= 0)
     &  vg(ind_gas_rev(51))=vg(ind_gas_rev(49)) ! INO2=HOI
      if(ind_gas_rev(56) /= 0)
     &  vg(ind_gas_rev(56))=0.                  ! CH3I          emission is net flux
      if(ind_gas_rev(57) /= 0)
     &  vg(ind_gas_rev(57))=0.                  ! CH2I2         emission is net flux
      if(ind_gas_rev(58) /= 0)
     &  vg(ind_gas_rev(58))=0.                  ! CH2ClI        emission is net flux
      if(ind_gas_rev(59) /= 0)
     &  vg(ind_gas_rev(59))=0.                  ! C3H7I         emission is net flux
      if(ind_gas_rev(63) /= 0 .and. ind_gas_rev(30) /= 0)
     &  vg(ind_gas_rev(63))=vg(ind_gas_rev(30)) ! CH3SO3H = HCl
      if(ind_gas_rev(71) /= 0)
     &  vg(ind_gas_rev(71))=0.                  ! CH2BrI         emission is net flux
      if(ind_gas_rev(72) /= 0)
     &  vg(ind_gas_rev(72))=0.                  ! CHBr2I         emission is net flux
      if(ind_gas_rev(73) /= 0)
     &  vg(ind_gas_rev(73))=0.                  ! C2H5I          emission is net flux


      if (lst/4*4.eq.lst.and.lmin.eq.1) then
         print *,lday,lst,lmin
         print *,' dry deposition velocities'
         do j=1,j1
            print *,j,vg(j)
         enddo
      endif

      do j=1,j1
! deposition, vg in m/s
         if (vg(j).ge.1.e-05) then
            s12old=s1(j,n_bl)
            s1(j,n_bl)=s1(j,n_bl)*exp(-dt/z_box*vg(j))
            s1(j,1)=s1(j,1)+(s12old-s1(j,n_bl))*z_box
         endif

! emission
! es1: emission rates in molec./cm**2/s, s1 in mol/m**3            
         s1(j,n_bl)=s1(j,n_bl)+es1(j)*dt*1.e+4/(z_box*Avogadro)
      enddo

      end subroutine sedc_box


!
!----------------------------------------------------------------
!

      subroutine box_partdep (dt, z_box, n_bl)
      
      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nka,
     &     nkt,
     &     nkc

      implicit double precision (a-h,o-z)

! dry deposition of particles and aqueous constituents in box
      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /kpp_vt/ vt(nkc,nf),vd(nkt,nka),vdm(nkc)

! calculation of deposition velocity is done in SR partdep; for smog chamber runs
! the roughness length z0 has to be adjusted in SRs box_init 

! apply v_dep to particles
      do ia = 1, nka
         do jt = 1, nkt
!            x = x0 * exp(-dt/z_box*vdep(l))
            ff_old = ff(jt,ia,n_bl)
            ff(jt,ia,n_bl) =  ff_old * exp(-dt/z_box*vd(jt,ia))
! add deposited numbers
            ff(jt,ia,1) = ff(jt,ia,1) + (ff_old - ff(jt,ia,n_bl))*z_box
         enddo
      enddo
! apply v_dep to non-ionic aqueous constituents
      do kc = 1, nkc
         x_depterm = exp(-dt/z_box*vdm(kc))
         do l = 1, j2
            s_old = sl1(l,kc,n_bl)
            sl1(l,kc,n_bl) = s_old * x_depterm
! add deposited numbers
            sl1(l,kc,1) = sl1(l,kc,1) + (s_old - sl1(l,kc,n_bl))*z_box
         enddo
      enddo
! apply v_dep to ionic aqueous constituents
      do kc = 1, nkc
         x_depterm = exp(-dt/z_box*vdm(kc))
         do l = 1, j6
            s_old = sion1(l,kc,n_bl)
            sion1(l,kc,n_bl) = s_old * x_depterm
! add deposited numbers
            sion1(l,kc,1) = sion1(l,kc,1) + (s_old - sion1(l,kc,n_bl))
     &           *z_box
         enddo
      enddo

      end subroutine box_partdep


!
!---------------------------------------------------------------------
!

      subroutine mic_init (iaertyp,fogtype)
! initialization of microphysics for restart (cut and paste from SR initm)
      
      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nb,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)
! initial profiles of meteorological variables
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb44/ g,a0m,b0m(nka),ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw
      double precision g,a0m,b0m,ug,vg,z0,ebs,psis,aks,
     &              bs,rhoc,rhocw,ebc,anu0,bs0,wmin,wmax,tw

      common /cb45/ u(n),v(n),w(n)
      common /cb47/ zb(nb),dzb(nb),dzbw(nb),tb(nb),eb(nb),ak(nb),d(nb),
     &              ajb,ajq,ajl,ajt,ajd,ajs,ds1,ds2,ajm,reif,tau,trdep
      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)
      double precision theta, thetl, t, talt, p, rho
      common /cb63/ fcs(nka),xmol3(nka)
      dimension sr(nka,nkt)
!      dimension xb0m(nka) ! jjb used for a test at the end, probably useless
      character *1 fogtype
      character *10 fname
!      data xmol2 /18./ ! jjb variable unreferenced (commented below)

! these write statements are in SR initm, has to be done here as well, to be consistent
! with plot progs
! initial output for plotting
      fname='pi .out'
      fname(3:3)=fogtype
      open (97, file=fname,status='unknown')
      write (97,6000) (eta(k),etw(k),rho(k),p(k),w(k),k=1,n)
 6000 format (5e16.8)
      close (97)
      write (19) zb
      close (19)
      do jt=1,nkt
!         jtp=min0(jt+1,nkt)  ! jjb variable unreferenced
!         de0=dew(jt)  ! jjb variable unreferenced
!         dep=dew(jtp) ! jjb variable unreferenced
!         de0p=de0+dep ! jjb variable unreferenced
         do ia=1,nka
            rk=rw(jt,ia)
            sr(ia,jt)=dmax1(.1d0,dexp(a0m/(rk*t(2))
     &                -b0m(ia)*en(ia)/ew(jt)))
         enddo
      enddo
      fname='fi .out'
      fname(3:3)=fogtype
      open (44, file=fname,status='unknown')
      write (44,6010) rn,en,rq,e,sr
      close (44)
 6010 format (5e16.8)

! parameters a0m, b0m of koehler curve of subroutine subkon: 
! sr(ia,jt)=exp(a0m/(r(ia,jt)*t)-b0m(ia)*en(ia)/ew(jt)))
! a0m see p. 141 pruppacher and klett a0m=2 sigma/(r1*t*rhow*a)
! 152200= 2 sigma*10**6 with sigma is surface tension = 76.1*10**-3
! see pruppacher and klett p. 104
!      a0m=152200./(r1*rhow) - is read in in SR startm
! aerosol types: 1=urban 2=rural 3=ocean 4=tropospheric
      k0=iaertyp
      do ia=1,nka
!         write (77,*)  ia
         if (k0-2) 2000,2010,2020
! b0m=fcs*xnue*xmol2/xmol3: fcs(ia) fraction of soluble species
! xnue number of ions; xmol2 (xmol3) mol masses of water (aerosol)
! NH4NO3 mole mass 80; (NH4)2SO4 mole mass 132
! soluble part of urban aerosol: 2 mole NH4NO3 and 1 mole (NH4)2SO4
 2000    continue
         fcs(ia)=.4-rn(ia)*(.4-.1)
         if (rn(ia).gt.1.) fcs(ia)=.1
         xmol3(ia)=(132.+80.*2.)/3.
         xnue=(3.+2.*2.)/3.
         go to 1030
! soluble part of rural aerosol: pure (NH4)2SO4
! 2010 fcs(ia)=.5-rn(ia)*(.5-.1)
 2010    continue
         fcs(ia)=.9-rn(ia)*(.9-.5)
         if (rn(ia).gt.1.) fcs(ia)=.5
!     fcs(ia)=0.1
         xnue=3.
         xmol3(ia)=132.
         go to 1030
!c soluble part of ocean aerosol: small pure (NH4)2SO4; large pure NaCl
! soluble part of ocean aerosol: pure (NH4)2SO4; 
 2020    continue
         fcs(ia)=1.
!      xnue=3.                   
!      xmol3(ia)=132.
! 32% (NH4)2SO4, 64% NH4HSO4, 4% NH4NO3
         xnue=0.32*3.+0.64*2+0.04*2                   
         xmol3(ia)=0.32*132.+0.64*115+0.04*80
! large are NaCl
         if (rn(ia).lt.0.5) go to 1030
!         xnue=2.                !no change in microphysics due to halogen chemistry ! jjb variable unreferenced
         xmol3(ia)=58.4
 1030    continue
!         xb0m(ia)=fcs(ia)*xnue*xmol2/xmol3(ia)
! 1030    b0m(ia)=fcs(ia)*xnue*xmol2/xmol3(ia)
!         write (77,*) xb0m(ia),b0m(ia) ! test shows: xb0m and b0m are the same
      enddo

      end subroutine mic_init

!
!----------------------------------------------------------------------
!

      subroutine out_mass
! subroutine to print aerosol and ion mass
      
      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     n,
     &     nka,
     &     nkt,
     &     nkc

       implicit double precision (a-h,o-z)

      parameter (lsp=9)
      common /cb40/ time,lday,lst,lmin,it,lcl,lct
      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      common /liq_pl/ nkc_l
      dimension xsum(n),lj2(lsp),xionmass(n,nkc),xion(n),xmm(lsp),
     &     zion(n)
      data lj2/1,2,8,9,13,14,19,20,30/
! molar mass of mass-determining ions in g/mole
      data xmm /1.,18.,96.,61.,62.,35.5,97.,23.,95./

      write (13,*) ' '
      write (13,*) 'output:', lday,lst,lmin

! Initialisation
      xxsum = 0.d0

! output of aerosol mass-------------
      do k=2,n
         xsum(k)=0.
         do ia=1,nka
            do jt=1,nkt
               xsum(k)=xsum(k)+ff(jt,ia,k)*en(ia)
            enddo
         enddo
         xsum(k)=xsum(k)*1.e+09
         xxsum=xxsum+xsum(k)*detw(k)
      enddo

      write (13,6240) 
      write (13,6250) xsum
      write (13,6260) xxsum

! output of ion mass-------------
      do k=1,n
         do kc=1,nkc_l
! calculate the mass
            xionmass(k,kc)=0.
            do l=1,lsp
               ll=lj2(l)
               xionmass(k,kc)=xionmass(k,kc)+ sion1(ll,kc,k)*xmm(l)
            enddo
          enddo

         xion(k)=xionmass(k,1)+xionmass(k,2)+xionmass(k,3)+
     &        xionmass(k,4)
      enddo

      do k=1,n
         zion(k) = xion(k)*1.d6
      enddo

      write (13,6280) 
      write (13,6250) zion

 6240 format (/,6x,'aerosol mass in ug m**-3 in layers 2 - nf')
 6250 format (1x,15f8.3)
 6260 format(6x,'total aerosol mass in ug m**-2 of layers 2 - nf',f12.3)
! 6270 format(4d16.8)
 6280 format(/,6x,'ion mass in ug m**-3 in layers 2 - nf')

      end subroutine out_mass

!
!-----------------------------------------------------------------------------
!

      subroutine get_n_box (z_box,nz_box)
      
      USE global_params, ONLY :
! Imported Parameters:
     &     n

      implicit double precision (a-h,o-z)
! get grid level for box height

      common /cb41/ detw(n),deta(n),eta(n),etw(n)
      double precision detw, deta, eta, etw


      nz_box=0
      do k=1,n
         if (etw(k).ge.z_box) then
            nz_box=k
            goto 6543
         endif
      enddo
 6543 continue

      if (nz_box.eq.0) then
         nz_box=70
         print *,"WARNING!! box height too high, set to ",etw(nz_box)
      endif

      print *,"box height set to ",etw(nz_box)
      print *,"box height corresponds to level ",nz_box

      end subroutine get_n_box

!
!-------------------------------------------------------------------------
!

      subroutine ffssuumm
! calculation of particle number

!      jjb cleaning
!          do end do without label

      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /blck06/ kw(nka),ka
!      common /kpp_kg/ vol2(nkc,n),vol1(n,nkc,nka),part_o
!     &     (n,nkc,nka),part_n(n,nkc,nka),pntot(nkc,n),kw(nka),ka
 

      do k=2,n
         fsum1=0.
         fsum2=0.
! small aerosol
         do ia=1,ka
            do jt=1,kw(ia)
               fsum1=fsum1+ff(jt,ia,k)
            enddo
         enddo
! large aerosol
         do ia=ka+1,nka
            do jt=1,kw(ia)
               fsum2=fsum2+ff(jt,ia,k)
            enddo
         enddo
         print *,k,fsum1,fsum2,fsum(k)
      end do

      end subroutine ffssuumm

!
!----------------------------------------------------------------
!

      subroutine oneD_dist_old
!  calculate 1D size distribution of 2D particles dist.

      USE constants, ONLY :
! Imported Parameters:
     & pi
      
      USE global_params, ONLY :
! Imported Parameters:
     &     n,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)
      double precision Np

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
      common /oneDs_0/partN(n,nka+nkt),partA(n,nka+nkt),partV(n,nka+nkt)
     &     ,partr(n,nka+nkt),drp(nka+nkt),nrp
      dimension rp(nka+nkt)  !particle radius [um]
      dimension Np(nka+nkt)  !particle number [part cm-3]
      dimension Ap(nka+nkt)  !particle surface [um2 cm-3]
      dimension Vp(nka+nkt)  !particle volume [um3 cm-3]

!     set up radius range problem: if rq(nkt,1) is used to map all other
!     radii on, every now and then 2 rq(jt,i) bins will fall into the
!     same rq(1:nkt,1) bin and produce artifial spikes in the size
!     distribution. The reason is that the radius resolution for small
!     ia is very fine (see plot of 2D particle grid). To avoid this use
!     rq(1,ia) until ia=nka and then rq(1:nkt,nka); the dimension of rq
!     is nka+nkt; to avoid having drp=0. the 


      rnmax = rn(nka)
      do ia=1,nka
         rp(ia)=rq(1,ia)
      enddo
      do jt=1,nkt
         jtt=jt
         if (rq(jt,nka).gt.rp(nka)*1.05) goto 1021
      enddo
 1021 continue
      do jt=1,nkt-jtt
         rp(nka+jt)=rq(jt+jtt,nka)
      enddo
      do ij=1,nka+nkt-1
         drp(ij)=rp(ij+1)-rp(ij)
!         print *,ij,rp(ij),drp(ij)
      enddo
      nrp=nka+nkt-jtt
      drp(nrp)=drp(nrp-1)

      do k=2,n
      do ij = 1,nrp
       Np(ij) = 0.
       Ap(ij) = 0.
       Vp(ij) = 0.
       fpi=4.*pi
       do ia = 1, nka
          if (rn(ia).gt.rp(ij)) goto 2001
          do jt = 1,nkt
             if (rq(jt,ia).le.rp(ij)) then
                if (ij.gt.1) then
                   if (rq(jt,ia).gt.rp(ij-1)) then
                      Np(ij) = Np(ij) + ff(jt,ia,k)
                      Ap(ij) = Ap(ij) + ff(jt,ia,k)*fpi*rp(ij)*rp(ij)
                      Vp(ij) = Vp(ij) + ff(jt,ia,k)*fpi*rp(ij)*
     &                     rp(ij)*rp(ij)/3.
                   endif
                else if ((ij.eq.1).and.(rq(jt+1,ia).gt.rp(ij))) then           ! jjb BUG here, jt+1 leads to out of bounds index when jt = nkt
!                write (*,100) ij,ia,jt,rq(jt,1),rp(ij),rq(jt,ia),ff(jt,ia,k)
                   Np(ij) = Np(ij) + ff(jt,ia,k)
                   Ap(ij) = Ap(ij) + ff(jt,ia,k)*fpi*rp(ij)*rp(ij)
                   Vp(ij) = Vp(ij) + ff(jt,ia,k)*fpi*rp(ij)*
     &                  rp(ij)*rp(ij)/3.
                endif
             else
                goto 2002          
             endif
          enddo                 !jt
 2002     continue
        enddo                    !ia
 2001  continue

! to plot dN/dr, dA/dr, dV/dr:
       partN(k,ij) = Np(ij)/drp(ij)
       partA(k,ij) = Ap(ij)/drp(ij)
       partV(k,ij) = Vp(ij)/drp(ij)
       partr(k,ij) = rp(ij)

      enddo                     !ij
      enddo                     ! k

! 100  format(3i4,4d16.8)

      end subroutine oneD_dist_old


!
!----------------------------------------------------------------
!

      subroutine oneD_dist
!  calculate 1D size distribution of 2D particles dist.
      
! jjb rewritten, but needs improvements.
!     idea: map all rq onto a specific which would range between rq(1,1) and rq(nkt,nka) with XXX values
!     (XXX could be equal to any of the already used parameters, nka, nkt, nka+nkt for finer resolution, or whatever)
!     currently, many particles are likely to be in the last rp class

      USE global_params, ONLY :
! Imported Parameters:
     &     nf,
     &     n,
     &     nka,
     &     nkt

      implicit double precision (a-h,o-z)
      double precision Np

      common /cb50/ enw(nka),ew(nkt),rn(nka),rw(nkt,nka),en(nka),
     &              e(nkt),dew(nkt),rq(nkt,nka)
      double precision enw,ew,rn,rw,en,e,dew,rq

      common /cb52/ ff(nkt,nka,n),fsum(n),nar(n)
!     common /oneDs/  partN(n,nkt,2),partr(n,nkt),drp(nkt),xlogdrp(nkt) ! jjb last variable shouldn't be in this CB
      common /oneDs/  partN(n,nkt,2),partr(n,nkt),drp(nkt)              ! jjb removed and declared below
      dimension xlogdrp(nkt) ! jjb

      dimension rp(nkt)  !particle radius [um]
      dimension Np(nkt)  !particle number [part cm-3]

! Ap, Vp can easily be calculated in ferret, so reduce output file size

! ignore radius range problem as outlined in SR oneD_dist_old and use rq(1:nkt,1)
! to map all particles onto; this might lead to adding a few bins of 2D spectrum 
! into one on 1D but I think this is just a "cosmetic" problem which in the end 
! might be overcome by smoothing the plot. Benefit: everything is easier and there 
! are no "empty" bins as when using the old radius 2D --> 1D mapping nka+nkt 

! define 1D grid
      do jt=1,nkt
         rp(jt)=rq(jt,1)
      enddo
! calculate width of each bin in 1D grid
!      drp(1)=rp(1)
!      xlogdrp(1)=log10(rp(1))
      do ij=1,nkt-1
         drp(ij)=rp(ij+1)-rp(ij)
! unit: "implicit" division of rp by 1um to get a unit-less property to be able to use log
         xlogdrp(ij)=log10(rp(ij+1))-log10(rp(ij))
!         print *,ij,rp(ij),drp(ij)
!         print *,drp(ij),xlogdrp(ij),log10(drp(ij))
      enddo
      drp(nkt)=rq(nkt,nka)-rp(nkt)
      xlogdrp(nkt)=log10(rq(nkt,nka))-log10(rp(nkt))
! map 2D spectrum on 1D spectrum
      do k=2,nf
         Np(:)=0.
         do ij = 1,nkt-1
            do ia = 1, nka
               if (rn(ia).gt.rp(ij+1)) goto 2001 ! save time
               do jt = 1,nkt
                  if (rq(jt,ia).lt.rp(ij+1)) then    
                     if (rq(jt,ia).ge.rp(ij)) then
                        Np(ij) = Np(ij) + ff(jt,ia,k)
                     endif
                  else
                     if (ij == nkt-1) then ! special case
                        Np(nkt) = Np(nkt) + ff(jt,ia,k)
                     else
                        goto 2002 ! save time
                     end if
                  endif
               enddo           !jt
 2002          continue
            enddo              !ia
 2001       continue
         end do                !ij

       do ij = 1,nkt
! to plot dN/dr, dA/dr, dV/dr:
          partN(k,ij,1) = Np(ij)/drp(ij)      ! #/um/cm3
! to plot dN/dlogr, dA/dlogr, dV/dlogr:
          partN(k,ij,2) = Np(ij)/xlogdrp(ij)  ! #/cm3
          partr(k,ij) = rp(ij)
!         print *,partN(k,ij),partA(k,ij),partV(k,ij),partr(k,ij)
       enddo                    !ij

! check if bins were "missed":
       xnsum=0.
       do jt=1,nkt
          xnsum=xnsum+Np(jt)
       enddo
       if (xnsum.lt..99*fsum(k)) print *,'N too small',k,fsum(k),xnsum
       if (xnsum.gt.1.01*fsum(k)) print *,'N too big',k,fsum(k),xnsum


      enddo                     ! k

! 100  format(3i4,4d16.8)

      end subroutine oneD_dist


!
!--------------------------------------------------------------------------
!


