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

      module kpp_KPP_ROOT_Parameters
      include 'KPP_ROOT_Parameters.h' ! KPP parameters
      end module kpp_KPP_ROOT_Parameters

      module kpp_KPP_ROOT_Global
      include 'KPP_ROOT_Parameters.h' ! KPP parameters
      include 'KPP_ROOT_Global.h'     ! KPP common blocs and additional user common blocks and other definitions
      end module kpp_KPP_ROOT_Global


      subroutine KPP_ROOT_drive
     &     (tkpp,dt_ch,k,xcvv1,xcvv2,yhal,yiod,
     &      yliq1,yliq2,yhet1,yhet2,air,h2o,xph_rat)

      USE constants, ONLY :
! Imported Parameters:
     &     xconv1=>conv1 ! multiply by conv1 to get cm^3(air)/mlc --> m^3(air)/mol

      USE gas_common, ONLY :
     &     j1, j5,
     &     s1,s3,
     &     gas_k2m_a, gas_m2k_a,
     &     rad_k2m_a, rad_m2k_a

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j3,
     &     j6,
     &     n,
     &     nf,
     &     nkc,
     &     nlev,
     &     nrxn

      implicit none

      include 'KPP_ROOT_Parameters.h' ! KPP parameters
      include 'KPP_ROOT_Global.h'     ! KPP common blocs and additional user common blocks and other definitions

! Subroutine arguments
      double precision tkpp,dt_ch,xcvv1,xcvv2,yhal,yiod,
     &      yliq1,yliq2,yhet1,yhet2,air,h2o

      double precision xph_rat(nphrxn)
      integer k

! Local scalar:
      integer kl ! do loop index to search if k belongs to il(nlev) list
      integer j

! Common blocks:
      common /blck12/ cw(nkc,n),cm(nkc,n)
      double precision cw, cm

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)
      double precision sl1, sion1

      common /budg/ bg(2,nrxn,nlev),il(nlev)
      double precision bg ! reaction rates (bg(1,:,:): instantaneous, bg(2,:,:): cumulative)
      integer il          ! indexes of the selected levels for reaction rates output

      common /kpp_laer/ henry(NSPEC,nf),xkmt(nf,nkc,NSPEC),
     &     xkef(nf,nkc,NSPEC),xkeb(nf,nkc,NSPEC)
      double precision henry, xkmt,xkef,xkeb

      common /kpp_drya/ xkmtd(nf,2,NSPEC),xeq(nf,NSPEC)
      double precision xkmtd,xeq

!      common /k_surf/ xkmt_OHClm(nf,nkc) ! jjb gamma_surf now commented, but keep this!

c parameters for /kpp_rate_a/
      cvv1=xcvv1
      cvv2=xcvv2
      xliq1=yliq1
      xliq2=yliq2
      xhet1=yhet1
      xhet2=yhet2
      xhal=yhal
      xiod=yiod
      conv1=xconv1
      ph_rat=xph_rat ! jjb

c the following data is needed only for heterogeneous reactions on dry aerosol
! jjb matrix below
c$$$      ycwd(1)=cwd(k,1)
c$$$      ycwd(2)=cwd(k,2)
c$$$      do l=1,nspec
c$$$         yxkmtd(1,l)=xkmtd(k,1,l)
c$$$         yxkmtd(2,l)=xkmtd(k,2,l)
c$$$         yxeq(l)    =xeq(k,l)
c$$$      enddo
      ycwd(:)=cw(:2,k) ! jjb matrix
      yxkmtd(:,:)=xkmtd(k,:,:) ! jjb matrix
      yxeq(:)    =xeq(k,:) ! jjb matrix

c liquid phase rates
      if (xliq1.eq.1..or.xliq2.eq.1.) then ! jjb this test is probably useless
! jjb matrix below
c$$$         do l=1,nspec
c$$$            yhenry(l)=henry(k,l)
c$$$         enddo
         yhenry(:)=henry(:,k) ! jjb matrix
      else
! jjb matrix below
c$$$         do l=1,nspec
c$$$            yhenry(l)=0.
c$$$         enddo
         yhenry(:)=0.d0 ! jjb matrix
      endif
      if (xliq1.eq.1.) then
         ycw(1)=cw(1,k)
! jjb matrix below
c$$$         do l=1,nspec
c$$$            yxkmt(1,l)=xkmt(k,1,l)
c$$$            ykef(1,l)=xkef(k,1,l)
c$$$            ykeb(1,l)=xkeb(k,1,l)
c$$$            ykmt_OHClm(1) = xkmt_OHClm(k,1) ! jjb shouldn't be in the do loop !
c$$$         enddo
         yxkmt(1,:)=xkmt(k,1,:) ! jjb matrix
         ykef(1,:)=xkef(k,1,:) ! jjb matrix
         ykeb(1,:)=xkeb(k,1,:) ! jjb matrix
!         ykmt_OHClm(1) = xkmt_OHClm(k,1) ! jjb gamma_surf now commented, but keep this!
      else
         ycw(1)=0.
! jjb matrix below
c$$$         do l=1,nspec
c$$$            yxkmt(1,l)=0.
c$$$            ykef(1,l)=0.
c$$$            ykeb(1,l)=0.
c$$$            ykmt_OHClm(1) = 0. ! jjb shouldn't be in the do loop !
c$$$         enddo
         yxkmt(1,:)=0. ! jjb matrix
         ykef(1,:)=0. ! jjb matrix
         ykeb(1,:)=0. ! jjb matrix
!         ykmt_OHClm(1) = 0. ! jjb gamma_surf now commented, but keep this!
      endif
      if (xliq2.eq.1.) then
         ycw(2)=cw(2,k)
! jjb matrix below
c$$$         do l=1,nspec
c$$$            yxkmt(2,l)=xkmt(k,2,l)
c$$$            ykef(2,l)=xkef(k,2,l)
c$$$            ykeb(2,l)=xkeb(k,2,l)
c$$$            ykmt_OHClm(2) = xkmt_OHClm(k,2) ! jjb shouldn't be in the do loop !
c$$$         enddo
         yxkmt(2,:)=xkmt(k,2,:) ! jjb matrix
         ykef(2,:)=xkef(k,2,:) ! jjb matrix
         ykeb(2,:)=xkeb(k,2,:) ! jjb matrix
!         ykmt_OHClm(2) = xkmt_OHClm(k,2) ! jjb gamma_surf now commented, but keep this!
      else
         ycw(2)=0.
! jjb matrix below
c$$$         do l=1,nspec
c$$$            yxkmt(2,l)=0.
c$$$            ykef(2,l)=0.
c$$$            ykeb(2,l)=0.
c$$$            ykmt_OHClm(2) = 0. ! jjb shouldn't be in the do loop !
c$$$         enddo
            yxkmt(2,:)=0. ! jjb matrix
            ykef(2,:)=0. ! jjb matrix
            ykeb(2,:)=0. ! jjb matrix
!            ykmt_OHClm(2) = 0. ! jjb gamma_surf now commented, but keep this!
      endif


c concentrations are handed over HERE (and not in separate SRs) because the 
c parameter (I_XXX) are different for each KPP block



c PRN2,PRPN,OZID are products, that don't react further, so no transport is needed
c maybe they are interesting as output ?? #

c include C(ind_)=s1/3(k,)
      do j=1,j1
         C(gas_m2k_a(1,j)) = s1(gas_m2k_a(2,j),k)
      end do

      do j=1,j5
         C(rad_m2k_a(1,j)) = s3(rad_m2k_a(2,j),k)
      end do

C#DEFFIX
      FIX(indf_O2)  = 0.21*air
      FIX(indf_H2O) = h2o
      FIX(indf_N2) = 0.79*air

c liquid phase
! jjb matrix below
c$$$         do kc=1,2 !nkc
c$$$            do l=1,j2
c$$$c               if (sl1(l,kc,k).lt.0) print *,k,'sl1(',l,kc,') < 0 !'
c$$$               sl1(l,kc,k)=max(0.d0,sl1(l,kc,k)) ! eliminate negative values
c$$$            enddo
c$$$            do l=1,j6
c$$$               sion1(l,kc,k)=max(0.d0,sion1(l,kc,k)) ! eliminate negative values
c$$$            enddo
c$$$         enddo
      sl1(:,:,k)=max(0.d0,sl1(:,:,k)) ! eliminate negative values
      sion1(:,:,k)=max(0.d0,sion1(:,:,k)) ! eliminate negative values

      if (cvv1.gt.0) then 
         FIX(indf_H2Ol1)=55.55/cvv1
      else
         FIX(indf_H2Ol1)=0.
      endif     
      if (cvv2.gt.0) then 
         FIX(indf_H2Ol2)=55.55/cvv2
      else
         FIX(indf_H2Ol2)=0.
      endif

c aerosol
c include C(ind_)=sl1/sion1(k,,1/2)
      include 'aer_mk.dat'


c integrate

      dt=dt_ch
      call Update_RCONST ()
      call INTEGRATE (tkpp,tkpp+dt_ch)
!        call bud_a (h2o,co,air,dt_ch,k) ! jjb h2o, co, air unused
!         call bud_a (dt_ch,k)

! Call budget subroutine only for selected levels
      do kl=1,nlev
         if(k.eq.il(kl)) then
            call bud_KPP_ROOT (dt_ch,kl)
            exit
         end if
      end do

! Call specific budget subroutine for all levels
      call bud_s_KPP_ROOT (dt_ch,k)

c hand-over concentrations: KPP --> MISTRA
c#DEFVAR
      do j=1,j1
         s1(j,k) = C(gas_k2m_a(j))
      end do
C#DEFRAD
      do j=1,j5
         s3(j,k) = C(rad_k2m_a(j))
      end do

c aerosol
c include sl1/sion1(k,,1/2)=C(ind_)
      include 'aer_km.dat'

      end subroutine KPP_ROOT_drive

