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

! activity.f:  calculation of activity coefficients after Pitzer

!  called by SR activ (in kpp.f)
!              |
!              |___SR pitzer
!                    |______SR calpar
!                    |______FN gammann
!                    |        |_______FN g
!                    |        |_______FN gs
!                    |        |_______SR efunct
!                    |
!                    |______FN gammasn
!                             |_______SR efunct
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Authors:
! --------
!    < Beiping Luo >

!      + contributions from:
!         < Roland von Glasow >  adapted for Mistra
!         < Josue Bock >         partly re-written, corrected, cleaned and improved*
!
!         * : technical improvements only, scientific content neither changed nor checked.

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine pitzer (k,kc,wact)
!
! Description:

!      This only a simple version, which only take the unsymmetrical
!      factor e-theta and E-theta ' into account.

!      Anion list:
!         1 HSO4-
!         2 SO4=
!         3 NO3-
!         4 Cl-

!      Cation list:
!         1 H+
!         2 NH4+
!         3 Na+
!

! Current Code Owner: released under GNU General Public License
!
! History:
! Version   Date     Comment
! -------   ----     -------
!           1990     Original code.                <Beiping Luo>
!
!           ?        Adapted for Mistra            <Roland von Glasow>
!
!         2016/2017  Multiple code improvement     <Josue Bock>
!                     to prepare release version:
!                     - general file header
!                     - all subroutines headers
!                     - removal of unused variables
!                     - all missing declarations and implicit none everywhere
!                     - lots of cleaning, properly indenting
!                     - integer vs real conflicts corrected in the whole module
!                     - use module for parameters
!                     - updated Fortran features such as intent, external
!                     - several bugfix

! Recent changes:
! ---------------
!    03/02/2017  <jjb>  xs array size is only 11 instead of 100, corrected in calpar as well

!    03/02/2017  <jjb>  bugfix: c1(3,:) was not defined along with b0, b1, and c0 (see after call calpar).
!                           Now set to 0. This was causing troubles when running the model with gfortran

!    03/02/2017  <jjb>  defined the dimension of all arrays using NC, NA instead of nmax (removed)

!    03/02/2017  <jjb>  use Mistra indexes for species from hardcoded "ionind" list

!    03/03/2017  <jjb>  more cleaning / consistency of variable naming with other places
!                       computing efficiency : if / elseif / endif  insead of repeated if ... endif
!                       za and zc are integers: changed declaration, use "real" to convert kind when necessary
!                       calculation of ZI, I and I2 only once, passed to SRs as arguments

!
! Code Description:
!   Language:          Fortran 77 (with Fortran 90 features)
!

! Declarations:
! Modules used:

      USE global_params, ONLY :
! Imported Parameters:
     &     j2,
     &     j6,
     &     nf,
     &     n,
     &     nkc

      implicit none

! Subroutine arguments
      integer, intent(in) :: k, kc ! Layer index, liquid bin index
      double precision, intent(out) :: wact

! Common blocks
      common /blck12/ cw(nkc,n),cm(nkc,n)                       ! LWC input
      double precision cw, cm

      common /blck17/ sl1(j2,nkc,n),sion1(j6,nkc,n)             ! ion mixing ratio
      double precision sl1, sion1

      common /cb53/ theta(n),thetl(n),t(n),talt(n),p(n),rho(n)  ! get temp from MISTRA
      double precision theta, thetl, t, talt, p, rho

      common /kpp_mol/ xgamma(nf,j6,nkc)                        ! gamma output
      double precision xgamma

! External functions:
      double precision, external :: gammann
      double precision, external :: gammasn   ! water activity

! Local parameters:
      integer,parameter :: NA = 4 ! number of anions
      integer,parameter :: NC = 3 ! number of cations
      double precision, parameter :: T0 = 298.15

! Local scalars:
      integer :: ia, ic                    ! loop indexes
      integer :: NI                        ! loop index, index of the ion whose activity is calculated
      double precision :: g
      double precision :: I, I2            ! I = ionic strength
      double precision :: rhmix1, r2mix1
      double precision :: TK               ! temperature in layer k (in [K])
      double precision :: xam, xap         ! sum of negative / positive charges (signed value)
      double precision :: xmix3, xu2, xu3
      double precision :: ZI               ! ZI = sum (mc*zc) + sum(ma*za)

! Local arrays:
      integer Iflag(nc,na)
      double precision :: b0(nc,na),B1(nc,na)  ! Pitzer coefficients
      double precision :: c0(nc,na),C1(nc,na)
      double precision :: MC(nc),MA(na)        ! concentration of cations and anions
      integer          :: ZC(nc),ZA(na)        ! ionic charge of cations and anions
      double precision :: omega(nc,na)
      double precision :: xs(11)

!      character ionn*6                                  ! jjb currently unused, but left commented for later use
!      dimension ionn(7)                                 !     could be tested along with Mistra ion names (index)
!      data ionn /'H+', 'NH4+', 'Na+', 'HSO4-','SO4--','NO3-', 'Cl-'/

      integer, parameter :: ionind(NC+NA) = (/1,2,20,19,8,13,14/) !indices of ions in sion1

      data ZC /1, 1, 1/
      data ZA /1, 2, 1, 1/

!- End of header ---------------------------------------------------------------


! -- Iflag definition, used in calpar
      Iflag(1,1)=3
      Iflag(1,2)=4
      Iflag(1,3)=1
      Iflag(1,4)=2

      Iflag(2,1)=5
      Iflag(2,2)=6
      Iflag(2,3)=7
      Iflag(2,4)=8

      Iflag(3,1)=0 ! <jjb> added this for security, even if never used.
      Iflag(3,2)=0 !        If a mistake was done in calpar, this is safer.
      Iflag(3,3)=0 !        See notes below for complete story
      Iflag(3,4)=0

      TK=t(k)            ! temperature in layer k (in [K])

! concentrations --> molality:
!  mol/m^3_air * m^3_air/m^3_solvent * 10^-3 * 1 dm^3/kg (density of water) = mol/kg_solvent
!   sion1      * cm^-1               * 1.d-3

! input molality of cations
      ! NB : cm /= 0. (already checked in the calling SR: activ, in kpp.f)
      do ic=1,nc
         mc(ic) = sion1(ionind(ic),kc,k)*1.d-3/cm(kc,k)
      enddo

! input of molality of anions
      do ia=1,na
         ma(ia) = sion1(ionind(NC+ia),kc,k)*1.d-3/cm(kc,k)
      enddo

! -- Calculate xap, xam, ZI and I
      I = 0.
      xap = 0.
      DO ic=1,NC
         I = I + mc(ic)*REAL(ZC(ic)**2)
         xap = xap + mc(ic)*REAL(zc(ic))
      ENDDO
      xam=0.
      DO IA=1,NA
         I = I + ma(ia)*REAL(ZA(ia)**2)
         xam = xam - ma(ia)*REAL(za(ia))
      ENDDO
      I=0.5*I
      I2=sqrt(I)

      ZI = xap - xam

!x! charge balance
!x      nadd=0
!x!     if( abs(xam+xap).ge.1.D-10) then
!x         if( abs(xam+xap).ge.1.D-5) then
!x 2010       continue
!x            nadd=nadd+1
!xcd            print *,xam,xap,abs(xam+xap)
!x            if (abs(xam).gt.abs(xap)) then
!x               xc=abs(xap/xam)
!x               do kn=1,na
!x                  ma(kn)=ma(kn)*xc
!x               enddo
!x            else
!x               xc=abs(xam/xap)
!x               do kn=1,nc
!x                  mc(kn)=mc(kn)*xc
!x               enddo
!x            endif
!x            do I=1,nc
!x               xap=xap+zc(I)*mc(I)
!xcd               print *,ii, mc(I)
!x             enddo
!x             do I=1,na
!x                xam=xam-za(I)*ma(I)
!xcd                print *,ii, mA(I)
!x             enddo
!x!             print *,xam,xap,abs(xam+xap),xc
!x!             print*, ' THE SOLUTION IS NOT NEUTRAL !'
!x!             if( abs(xam+xap).ge.1.D-10) goto 2010
!x             if( (1-xc).ge.1.D-4) goto 2010
!x!             goto 1000
!x         endif


      call calpar(Iflag,NC,NA,TK, b0,b1,c0,c1,omega,xs) ! <jjb> note NC is used here, to get array size in calpar
                                                        !      however, NC=3 values are set below, not in calpar
                                                        !      thus, NC-1=2 is used in calpar do loop.

! get the parameters for Na-anions (NA=3)
      ! -- Na-HSO4 ! Die Temperatureabhaengigkeit fuer Na-HSO4 fehlt noch.
      b0(3,1)=.0454
      b1(3,1)=.398
      c0(3,1)=0.d0
      c1(3,1)=0.d0
      omega(3,1)=2.d0
      ! -- Na-SO4
      b0(3,2)=0.0261 + (TK-T0)* 2.36d-3
      b1(3,2)=1.484 + (TK-T0)* 5.63d-3
      c0(3,2)=.00938 - (TK-T0)* .172d-3
      c1(3,2)=0.d0
      omega(3,2)=2.d0
      ! -- Na-NO3
      b0(3,3)=.0068 + (TK-T0)* 12.66d-4
      b1(3,3)=.1783 + (TK-T0)* 20.6d-4
      c0(3,3)=-.00072/2d0 - (TK-T0)* 23.16d-5/2d0
      c1(3,3)=0.d0
      omega(3,3)=2.d0
      ! -- NaCl
      b0(3,4)=0.0765d0 + (TK-T0)* 7.159d-4
      b1(3,4)=.2664d0 + (TK-T0)* 7d-4
      c0(3,4)=.00127/2d0 - (TK-T0)* 10.5d-5/2d0
      c1(3,4)=0.d0
      omega(3,4)=2.d0

! -- Calculate ions activities
      DO NI=1, NC+NA
         ! gammann calculate on the binary terms
         g=gammann(TK,NI,NC,NA,mC,mA,zC,zA,ZI,I,I2,b0,b1,C0,C1,omega)

         ! calculate the ternary terms of H+
         if( NI.eq.1) then
            xmix3= ma(2)*mc(2)*xs(7)+mC(2)*mA(1)*xs(9)
            rhmix1=ma(1)*ma(3)*xs(1)+xs(3)*ma(2)*ma(3)
            g=g*exp(rhmix1)

         ! calculate the ternary terms of NH4+
         else if(NI.eq.2) then
            xmix3 = ma(2)*mc(1)*xs(8)+mA(1)*ma(2)*xs(7)
     &              +ma(1)*mc(1)*xs(9)
            g=g*exp((xmix3)+2*mc(1)*xs(10))


         ! calculate the ternary terms of HSO4-
         else if(NI-NC.eq. 1) then
            xu2=mc(1)*ma(3)*xs(1)
            xu2=xu2+mc(1)*ma(4)*xs(4) + ma(4)*xs(5)*2
            ! mixture of NH4,SO4,HSO4
            xmix3= ma(2)*mc(2)*xs(7)+mc(1)*MC(2)*xs(9)
            g=g*exp(xu2+xmix3)

         ! calculate the ternary terms of SO4--
         else if(NI-NC.eq. 2) then
            xu3=ma(3)*mc(1)*xs(3) + ma(3)*xs(2)*2
            xu3=xu3+ma(4)*mc(1)*xs(6)
            xmix3= ma(1)*mc(2)*xs(7) +mc(1)*mc(2)*xs(8)
            g=g*exp(xu3+xmix3)

         ! calculate the ternary terms of NO3-
         else if(NI-NC.eq. 3) then
            r2mix1=ma(1)*mc(1)*xs(1) !H,HSO4,NO3
     &             + 2.*ma(2)*xs(2)      !SO4,NO3
     &             + ma(2)*mc(1)*xs(3) ! H,SO4,NO3
     &             + mc(2)*ma(2)*xs(11) ! NH4,SO4,NO3
            g=g*exp(r2mix1)
         endif

         xgamma(k,ionind(NI),kc)=g

      END DO

! -- Calculate water activity      
      wact=GAMMASN(TK,NC,NA,MC,MA,ZC,ZA,ZI,I,I2,B0,B1,C0,C1,omega)


      end subroutine pitzer

!------------------------------------------------------------------

      subroutine calpar(Iflag,NC,NA,T,b0,b1,C0,C1,omega,xs)

! 03/02/2017 <jjb> corrected array sizes:
!     - bb(27) instead of bb(50)
!     - b2(21) instead of b2(50)
!     - b3(42) instead of b3(50)
!     - b5(31) instead of b5(50) ! notice unused indexes
!     - b6(13) instead of b6(50)
!     - b7(13) instead of b7(50)
!     - b8(9)  instead of b8(50)
!     - b9(17) instead of b9(50) ! unused, thus commented
!
! 03/02/2017 <jjb> aligned data blocks for clarity and ease futur check
!                  removed duplicated definitions of dt inside the if cases
!                  declared all variables, implicit none
!
! 03/02/2017 <jjb> defined the dimension of all arrays using NC, NA instead of a larger nmax (removed)
!                  changed the main do-loop bound, from NC to NC-1 (previously, calpar was called with
!                  a "fake" NC = NC-1). Now, called with the "real" NC value.
!                  Data for NC = 3 are in SR Pitzer.
!
! 03/03/2017 <jjb> xs(1:3) and xs(11) are defined twice the same (H-NO3 and NH4-Cl cases)
!                  commented the second, no need to overwrite with same values
!                  Also removed the initialisation (data xs /11*0./) in pitzer header: useless
!
! 04/03/2017 <jjb> introduced dt2, dt3, dt4 to avoid multiple ** calculations

      implicit none

! Subroutine arguments
! Scalar arguments with intent(in):
      integer, intent(in) :: NA, NC
      double precision, intent(in)  :: T
! Array arguments with intent(in):
      integer, intent(in) :: Iflag(NC,NA)

! Array arguments with intent(out):
      double precision, dimension(nc,na), intent(out) :: b0, b1
      double precision, dimension(nc,na), intent(out) :: c0, c1

      double precision, intent(out) :: omega(nc,na)
      double precision, intent(out) :: xs(11)

! Local parameters:
      double precision, parameter :: T0 = 298.15

! Local scalars:
      double precision :: dt, dt2, dt3, dt4
      integer :: i, j

! Local arrays:
      double precision :: bb(27)
      double precision :: b2(21)
      double precision :: b3(42)
      double precision :: b5(31)
      double precision :: b6(13)
      double precision :: b7(13)
      double precision :: b8(9)
!      double precision :: b9(17) ! <jjb> unused

!- End of header ---------------------------------------------------------------


! ----------- data for H-NO3 -------------------------------
      data (bb(I),I=1,27)/
     *  3.895835D-03, -1.55571D-02,   1.703729D-02, -5.6173712D-03,
     *  5.732047D-03,  0.91622,       0.613523,     -0.68489,
     *  0.3038,       -0.32888,       7.6086113D-07, 7.2714678D-05,
     * -1.0037D-04,    3.475D-05,    -3.62927D-05,   5.380465D-02,
     * -2.2163D-02,   -1.0166D-02,    6.5423D-03,   -8.80248D-03,
     *  0.907342,     -6.78428D-4,    9.576D-4,      0.D0,
     *  0.D0,          7.769D-3,     -5.819D-4 /

! ---------- data for H-Cl ---------------------------------
      data (b2(I),i=1,21) /
     *  0.23378,      -7.21238D-02,  -1.7335667D-02, 5.760665D-03,
     * -8.29279D-03,   0.2897,        7.575434D-02, -1.1474D-03,
     *  0.38038,      -0.309442,     -2.794885D-03,  2.309349D-04,
     *  9.322982D-04, -2.398D-04,     2.85959D-04,  -0.21154,
     *  0.101481,      5.945618D-02, -0.107864,      8.81749D-02,
     *  1.9916 /

! ----------- data for H-HSO4, H-SO4 -----------------------
      data (b3(I),i=1,42) /           0.148843,     -7.769D-2,
     *  2.8062D-2,     4.7903D-4,     7.25D-4,       0.17843,
     *  0.678,         8.7381D-2,    -0.57881,       7.58D-2,
     * -9.878D-4,      5.447651D-4,  -2.58798D-4,    1.8466527D-5,
     *  1.23457D-5,    0.37138,      -9.24874D-2,   -9.21372D-3,
     * -1.065158D-2,   5.4987733D-2,  0.2726312,    -1.34824D-3,
     * -0.24711,       1.25978D-2,    0.11919,       0.7397,
     * -3.01755,      -4.5305,       -3.1072,       -0.8555842,
     *  9.2223D-4,    -4.1694532D-3,  7.141266D-3,   2.32984D-3,
     * -6.98191D-4,   -2.242,         0.71925,       2.52,
     * -0.7391,       -1.548503,      1.5452,        2. /

! ----------- data for NH4-HSO4 ----------------------------
! -- without chan's data
      data (b5(I),I=1,31)/
     * -8.746D-4,  -2.3125,    -9.56785D-6, 2.58238,    2.38,
     * -3.1314D-4,  1.6896D-2, -0.7351,     0.6883,     1.813D-3,  ! indexes 8 & 9 never used
     * -0.1012515, -2.66D-2,   -2.86617D-3, 0.22925,    0.438188,
     *  2.522D-4,  -2.90117D-5, 0.9014,     0.41774,   -1035.9,    ! indexes 20 - 23 never used
     *  0.0,       -299.69,     0.0,       -4.9687D-4,  0.0,       !
     *  1.21485D-2, 0.0,       -1.0334D-3,  0.0,        8.48374D-2,
     *  0.0 /

! -- with chan's data
!      data (b5(I),I=1,31)/
!     * -7.8224D-3,  -1.722,      7.5882D-5,  0.962,      1.8914,
!     * -4.285D-4,    3.99D-4,   -0.76151,    0.53133,    1.16726D-3, ! indexes 8 & 9 never used
!     * -9.45855D-2, -2.007D-2,   6.0482D-3, -0.1392,     2.9544,
!     *  2.45133D-4, -5.492D-5,   0.8358,    -0.3791,    -980.81,     ! indexes 20 - 23 never used
!     *  0.,          265.7,      0.,        -1.7568D-3, -7.72D-4,    !
!     *  3.437D-4,   -1.6456D-2, -1.2034D-3, -2.2943D-3,  7.9281D-2,
!     *  3.772D-2 /

! ----------- data for NH4-SO4 -----------------------------
      data (b6(I),I=1,13) /     -1.2058223D-2, 1.1043,     4.79018D-5,
     *  2.14346D-2, 0.58,       -2.9146D-2,    1.9631D-4,  1.1378,
     *  0.9283,     1.28548D-4,  1.684D-5,     2.6267D-2, -2.6D-4 /    ! wt=1


!      data (b6(I),I=1,13)/     -1.2058223D-2,  1.1043,    4.79018D-5,
!     *  2.14346D-2,  0.58,     -0.1188,        8.5D-2,    2.10514,
!     *  0.5942,      7.888D-4, -5.503D-4,      5.815D-2, -3.766D-2 /  ! wt=5d-5


! ----------- data for NH4-NO3 -----------------------------
      data (b7(I),I=1,13) /-2.3275D-2, 0.15, 1.1634D-4, 1.62D-3,
     *  0.43, ! 0.107264,  0.221, -3.8442D-4, -1.3872D-2, 4*0D0/
     *  8.78D-2 ,  0.2753645, -3.349D-4 , -1.093D-2 ,
     * -4.769D-2 , 0.1776,     1.25D-4 ,   6.9751D-3 /

!     *  2.15438D-2,0.67073,  2*0.0, -1.3662D-2,3.4747D-002,2*0.0/


! ----------- data for NH4-Cl ------------------------------
      data (b8(I),I=1,9) / -6.333D-4,  -3.99546D-4, 0.3155,
     *  0.1414, -3.837D-5,  1.08331D-4, 5.2436D-2,  1.6827D-2, 1.19 /

!! ----------- data for Na-SO4 ------------------------------              ! <jjb> unused
!      data (b9(I),I=1,17) /       1.63D-03,    2.092D-03,  3.484156D-02,
!     * -1.057D-02,   1.0775,      0.9,         0.8206,    -9.6425D-02,
!     *  2.7492D-03, -9.2838D-04, -6.8268D-04,  4.8126D-05, 7.182D-02,
!     *  0.7586,     -0.30291,     0.2311,1.7/

      xs(4)=0.d0
      xs(5)=0.d0
      xs(6)=0.d0

! -- Common definition
      dt=(T-T0)/100.D0
      dt2 = dt*dt
      dt3 = dt*dt2
      dt4 = dt2*dt2

      DO I=1,NC-1   ! Note the use of NC-1 here. Case NC=3: data are already in Pitzer subroutine
         DO J=1,NA
            ! -- Iflag == 1 : H-NO3
            if ( Iflag(I,J) .eq. 1) then
               B0(i,j)=bb(1)+dt*bb(2)+dt2*bb(3)+dt3*bb(4)+dt4*bb(5)
               B1(i,j)=bb(6)+dt*bb(7)+dt2*bb(8)+dt3*bb(9)+dt4*bb(10)
               c0(i,j)=bb(11)+dt*bb(12)+dt2*bb(13)+dt3*bb(14)+dt4*bb(15)
               c1(i,j)=bb(16)+dt*bb(17)+dt2*bb(18)+dt3*bb(19)+dt4*bb(20)

               omega(i,j)=bb(21)

               xs(1)=bb(22)+bb(23)*dt
               xs(2)=bb(24)+bb(25)*dt
               xs(3)=bb(26)+bb(27)*dt
               xs(11)=4.75458D-4 - 4.0577D-3*dt
                     !4.633D-4-4.093D-3*dt

            ! -- Iflag == 2 : H-Cl
            else if (Iflag(I,J) .eq. 2) then
               B0(i,j)=b2(1)+dt*b2(2)+dt2*b2(3)+dt3*b2(4)+dt4*b2(5)
               B1(i,j)=b2(6)+dt*b2(7)+dt2*b2(8)+dt3*b2(9)+dt4*b2(10)
               c0(i,j)=b2(11)+dt*b2(12)+dt2*b2(13)+dt3*b2(14)+dt4*b2(15)
               c1(i,j)=b2(16)+dt*b2(17)+dt2*b2(18)+dt3*b2(19)+dt4*b2(20)

               omega(i,j)=b2(21)


            ! -- Iflag == 3 : H-HSO4  ( and Iflag == 4 : H-SO4 )
            else if (Iflag(I,j).eq.3) then ! NB: Iflag(i,j) == 3 when i=1, j=1.
                                           !    this case also fills arrays(1,2) which is Iflag == 4
               B0(i,j)=b3(1)+dt*b3(2)+dt2*b3(3)+dt3*b3(4)+dt4*b3(5)
               B1(i,j)=b3(6)+dt*b3(7)+dt2*b3(8)+dt3*b3(9)+dt4*b3(10)
               c0(i,j)=b3(11)+dt*b3(12)+dt2*b3(13)+dt3*b3(14)+dt4*b3(15)
               c1(i,j)=b3(16)+dt*b3(17)+dt2*b3(18)+dt3*b3(19)+dt4*b3(20)

               B0(i,j+1)=b3(21)+dt*b3(22)+dt2*b3(23)+dt3*b3(24)
     &                   +dt4*b3(25)
               B1(i,j+1)=b3(26)+dt*b3(27)+dt2*b3(28)+dt3*b3(29)
     &                   +dt4*b3(30)
               c0(i,j+1)=b3(31)+dt*b3(32)+dt2*b3(33)+dt3*b3(34)
     &                   +dt4*b3(35)
               c1(i,j+1)=b3(36)+dt*b3(37)+dt2*b3(38)+dt3*b3(39)
     &                   +dt4*b3(40)

               omega(I,J)=b3(41)
               omega(I,J+1)=b3(42)

            ! -- Iflag == 5 : NH4-HSO4
            else if ( Iflag(I,J) .eq. 5) then
               B0(I,J)=b5(1)+b5(11+1)*dt+b5(11+2)*dt2
               B1(I,J)=b5(2)+b5(11+3)*dt+b5(11+4)*dt2
               c0(I,J)=b5(3)+b5(11+5)*dt+b5(11+6)*dt2
               c1(I,J)=b5(4)+b5(11+7)*dt+b5(11+8)*dt2

               omega(I,J)=b5(5)

               xs(7)=b5(6)+b5(11+13)*dt+dt2*b5(11+14)
               xs(8)=b5(7)+b5(11+15)*dt+dt2*b5(11+16)
               xs(9)=b5(10)+b5(11+17)*dt+dt2*b5(11+18)
               xs(10)=b5(11)+b5(11+19)*dt+dt2*b5(11+20)

            ! -- Iflag == 6 : NH4-SO4
            else if ( Iflag(I,J) .eq. 6) then
               b0(I,J)=b6(1)+b6(6)*dt+b6(7)*dt2
               b1(I,J)=b6(2)+b6(8)*dt+b6(9)*dt2
               C0(I,J)=b6(3)+b6(10)*dt+b6(11)*dt2
               C1(I,J)=b6(4)+b6(12)*dt+b6(13)*dt2

               omega(I,J)=b6(5)

            ! -- Iflag == 7 : NH4-NO3
            else  IF (Iflag(I,j).eq.7) then
               b0(I,J)=b7(1)+b7(6)*dt+b7(10)*dt2
               b1(I,J)=b7(2)+b7(7)*dt+b7(11)*dt2
               c0(I,J)=b7(3)+b7(8)*dt+b7(12)*dt2
               c1(I,J)=b7(4)+b7(9)*dt+b7(13)*dt2

               omega(I,j)=b7(5)

               !xs(1)=bb(22)+bb(23)*dt            ! <jjb> already defined (identical) in H-NO3 case
               !xs(2)=bb(24)+bb(25)*dt
               !xs(3)=bb(26)+bb(27)*dt
               !xs(11)=4.75458D-4-4.0577D-003*dt
               !      !4.633D-4-4.093D-3*dt

            ! -- Iflag == 8 : NH4-Cl
            else  IF (Iflag(I,j).eq.8) then
               b0(I,J)=b8(1)+b8(2)*dt
               b1(I,J)=b8(3)+b8(4)*dt
               c0(I,J)=b8(5)+b8(6)*dt
               c1(I,J)=b8(7)+b8(8)*dt

               omega(I,j)=b8(9)

!            else  IF (Iflag(I,J).eq.9) then                       ! <jjb> unused. See SR pitzer
!               B0(i,j)=b9( 1)+dt*b9( 2)+dt2*b9( 3)+dt3*b9( 4)
!               B1(i,j)=b9( 5)+dt*b9( 6)+dt2*b9( 7)+dt3*b9( 8)
!               c0(i,j)=b9( 9)+dt*b9(10)+dt2*b9(11)+dt3*b9(12)
!               c1(i,j)=b9(13)+dt*b9(14)+dt2*b9(15)+dt3*b9(16)
!
!               omega(I,j)=b9(17)

            endif
         enddo
      enddo

      end subroutine calpar

!-------------------------------------------------------------

      function gammann(T,N,NC,NA,mC,mA,ZC,ZA,ZI,I,I2,B0,B1,C0,C1,omega)
!
!      This only a simple version, which only take the unsymmetrical
!      factor e-theta and E-theta ' into account.
!

! 25/01/2017 <jjb> cleaned, properly indented. Later: rewritten without labels, goto
!
! 03/02/2017 <jjb> changed nested do loops order for efficiency (4 places)
!                  za and zc are integers: changed declaration, use "real" to convert kind when necessary
!                  imported ZI from calling subroutine, instead of recalculating it
!
! 03/03/2017 <jjb> added a test to avoid division by zero: ABS(XO4) > TINY(0.)
!                  this case should already be handled in calling subroutine activ, but left here for security


      implicit none

! Function declaration
      double precision :: gammann

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) :: T      ! temperature, in K
      integer, intent(in) :: N               ! the index of the ion in the C+A list
      integer, intent(in) :: NA, NC          ! number of anion and cation species
! Array arguments with intent(in):
      double precision, intent(in) :: mC(nc), mA(na)      ! molality of cations and anions
      integer,          intent(in) :: ZC(nc), ZA(na)      ! ionic charge of cations and anions
      double precision, intent(in) :: ZI                  ! ZI = sum (zc*mc) + sum(za*ma)
      double precision, intent(in) :: I, I2               ! I = sum (mc*zc**2) + sum(ma*za**2) and I2=sqrt(I)
      double precision, intent(in) :: B0(nc,na),B1(nc,na) ! Pitzer coefficients
      double precision, intent(in) :: C0(nc,na),C1(nc,na)
      double precision, intent(in) :: omega(nc,na)

! External functions:
      double precision, external :: g, gs

! Local scalars:
      double precision :: alpha, Aphi
      double precision :: a1, a2, a3
      double precision :: E, ED
      double precision :: F1, F2, F3, F4, F
      double precision :: gam
      double precision :: gg, ggs
      double precision :: omega1
      double precision :: x, xo, xo4, xhx, xhxs
      integer :: icharg, jcharg                 ! arguments passed to SR efunc
      integer :: IA,IA1,IA2, IC,IC1,IC2         ! loop indexes
      integer :: N1
      integer :: Z1, Z2

! Local arrays:
      double precision :: B(nc,na),Bs(nc,na)
      double precision :: C(nc,na),Cs(nc,na)

!- End of header ---------------------------------------------------------------


! --  Calculate B, Bs, C, Cs
      alpha=2.
      x=I2*alpha
      gg=g(x)
      ggs=gs(x)

      DO IA = 1, Na
         DO IC = 1, NC
            B(IC,IA) = b0(IC,Ia) + gg* b1(IC,Ia)
            BS(IC,IA) = ggs*b1(IC,Ia)/I

            omega1=OMEGA(IC,IA)
            xo=omega1*I2
            xo4 = xo*xo*xo*xo

            IF(ABS(XO4) > TINY(0.d0)) THEN
               !xhx=1./xo**4*(6.-exp(-xo)*(6.+6.*xo+3.*xo**2+xo**3))
               xhx=1./xo4*(6.-exp(-xo)*(6.+6.*xo+3.*xo**2+xo**3))
            ELSE
               PRINT*,'Warning: in SR gammann, xo4 too small'
               xhx=0.d0
            ENDIF
            xhxs=exp(-xo)/2.-2.*xhx

            C(IC,Ia)=(C0(IC,Ia)+4.*C1(IC,Ia)*xhx)
            Cs(IC,Ia)= C1(IC,Ia)/I*xhxs
         ENDDO
      ENDDO

      Aphi=.377+4.684d-4*(T-273.15)+3.74d-6*(T-273.15)**2

      F1= -Aphi*(I2/(1.+1.2*I2) +2./1.2*log(1.+1.2*I2))

      F2=0.
      DO IA = 1, NA
         DO IC = 1, NC
            F2=F2+mc(IC)*Ma(Ia)*(BS(IC,Ia)+2.*ZI*CS(IC,IA) )
         ENDDO
      ENDDO

      ICHARG=1
      JCHARG=2
      CALL EFUNC(ICHARG,JCHARG,Aphi,I,E,ED)

      F3= 0.
      DO IC1 = 1, NC
         DO IC2 = IC1+1, NC
            z1=ZC(IC1)
            z2=ZC(IC2)
            IF(Z1.NE.Z2) F3=F3+ED*MC(IC1)*MC(IC2)
         ENDDO
      ENDDO

      F4=0.
      DO IA1 = 1, NA
         DO IA2 = IA1+1, NA
            Z1=ZA(IA1)
            Z2=ZA(IA2)
            IF(Z1.NE.Z2) F4=F4+ED*MA(IA1)*MA(IA2)
         ENDDO
      ENDDO

      F=F1+F2+F3+F4

      ! -- Cations
      IF (N.LE.NC) THEN
         a1=REAL(ZC(N)**2) *F

         a2=0.
         DO Ia=1,NA
            a2=a2+MA(IA)*( 2*B(N,IA)+ZI*C(N,IA) )
         ENDDO

         a3=0.
         DO IA=1,NA
            DO IC=1,NC
               a3=a3+MA(IA)*MC(IC)*C(IC,IA)
            ENDDO
         ENDDO
         a3=REAL(ZC(N))*a3

         F3=0.
         DO IC = 1, NC
            Z1=ZC(N)
            Z2=ZC(IC)
            IF(Z1.NE.Z2) F3=F3+E*MC(IC)
         ENDDO

         gam=a1+a2+a3 +F3

      ! -- Anions
      ELSE
         N1=N-NC

         a1=REAL(ZA(N1)**2) *F

         a2=0.
         DO IC=1,NC
            a2=a2+MC(IC)*( 2.*B(IC,N1)+ZI*C(IC,N1) )
         ENDDO

         a3=0.
         DO IA=1,NA
            DO IC=1,NC
               a3=a3+MA(IA)*MC(IC)*C(IC,IA)
            ENDDO
         ENDDO
         a3=REAL(ZA(N1))*a3

         F4=0.
         DO IA = 1, NA
            Z1=ZA(N1)
            Z2=ZA(IA)
            IF(Z1.NE.Z2) F4=F4+E*MA(IA)
         ENDDO

         gam=a1+a2+a3+F4

      ENDIF

      gammann= exp(gam)

      end function gammann

!-------------------------------------------------------------

      function g(x)

! 01/2017 <jjb> implicit none

      implicit none
      double precision :: g,x

      g=2*(1.-(1.+x)*exp(-x))/x**2
      end function g

!-------------------------------------------------------------

      function gs(x)

! 01/2017 <jjb> implicit none

      implicit none
      double precision :: gs, x

      gs=2.*(-1.+(1.+x+x**2/2.)*exp(-x) )/x**2
      end function gs

!-------------------------------------------------------------

      SUBROUTINE EFUNC(ICHARG,JCHARG,APHI,XI,E,ED)

! 03/02/2017 <jjb> declared all variables + organised header
!                  removed internal variables XE, XED which were identical to E, ED
!                  renamed A -> APHI to be consistent with calling name
!                  in calling SRs, renamed J1 -> ICHARG and J2 -> JCHARG

      IMPLICIT NONE

! Subroutine arguments
! Scalar arguments with intent(in):
      INTEGER, INTENT(IN) :: ICHARG, JCHARG
      DOUBLE PRECISION, INTENT(IN) :: APHI, XI
! Scalar arguments with intent(out):
      DOUBLE PRECISION, INTENT(OUT) :: E, ED

! Local scalars:
      INTEGER :: I
      DOUBLE PRECISION :: DUM
! Local arrays:
      DOUBLE PRECISION :: X(3)
      DOUBLE PRECISION :: J0(3), J1(3)
!- End of header ---------------------------------------------------------------

      IF((ICHARG.EQ.JCHARG) .OR. (XI .LE. 1.D-30) )THEN
         E=0.
         ED=0.
      ELSE
         X(1)=6.*ICHARG*JCHARG*SQRT(XI)*APHI
         X(2)=6.*ICHARG*ICHARG*SQRT(XI)*APHI
         X(3)=6.*JCHARG*JCHARG*SQRT(XI)*APHI

         DO I=1,3
            DUM=-1.2D-2*X(I)**.528
            J0(I)=X(I)/(4.+4.581*X(I)**(-.7238)*EXP(DUM))
            J1(I)=(4.+4.581*X(I)**(-.7238)*EXP(DUM)*(1.7238-DUM*.528))
     &            /(4.+4.581*X(I)**(-.7238)*EXP(DUM))**2
         ENDDO

         E=ICHARG*JCHARG/4./XI*(J0(1)-.5*J0(2)-.5*J0(3))
         ED=ICHARG*JCHARG/8./XI**2*(X(1)*J1(1)-.5*X(2)*J1(2)-.5*
     &   X(3)*J1(3))-E/XI
      ENDIF

      END SUBROUTINE EFUNC

!------------------------------------------------------------------

      FUNCTION GAMMASN(T,NC,NA,MC,MA,ZC,ZA,ZI,I,I2,B0,B1,C0,C1,OMEGA)

!      GAMMASN === WATER ACTIVITY = PH2O/PH2O(PURE WATER)
!      ALL ARGUMENTS ARE INPUTS
!      C1 == 0.
!      OMEGA =1.


! 03/02/2017 <jjb> changed nested do loops order for efficiency (2 places)
!                  removed unreferenced variable xmh0
!                  za and zc are integers: changed declaration, use "real" to convert kind when necessary
!                  imported ZI, I and I2 from calling subroutine, instead of recalculating it


      IMPLICIT NONE

! Function declaration:
      double precision :: gammasn

! Subroutine arguments
! Scalar arguments with intent(in):
      double precision, intent(in) :: T      ! temperature, in K
      integer, intent(in) :: NA, NC          ! number of anion and cation species
! Array arguments with intent(in):
      double precision, intent(in) :: MC(nc), MA(na)      ! molality of cations and anions
      integer,          intent(in) :: ZC(nc), ZA(na)      ! ionic charge of cations and anions
      double precision, intent(in) :: ZI                  ! ZI = sum (mc*zc) + sum(ma*za)
      double precision, intent(in) :: I, I2               ! I = sum (mc*zc**2) + sum(ma*za**2) and I2=sqrt(I)
      double precision, intent(in) :: B0(nc,na),B1(nc,na) ! Pitzer coefficients
      double precision, intent(in) :: C0(nc,na),C1(nc,na)
      double precision, intent(in) :: OMEGA(nc,na)


! Local scalars:
      double precision :: ALPHA, APHI, AS
      double precision :: E, ED
      double precision :: F3, F4, FPHI1
      double precision :: OMEGA1
      double precision :: PP, PHI, PHIX
      double precision :: X, XO, XMI, XS
      integer :: icharg, jcharg            ! arguments passed to SR efunc
      integer :: IA,IA1,IA2, IC,IC1,IC2    ! loop indexes
      integer :: Z1, Z2

! Local arrays:
      double precision :: BPHI(NC,NA),C(NC,NA)

!- End of header ---------------------------------------------------------------


! --  CALCULATE BPHI, C
      ALPHA=2.
      X=I2*ALPHA
      DO IA = 1, NA
         DO IC = 1, NC
            BPHI(IC,IA) = B0(IC,IA) + EXP(-X)*B1(IC,IA)
            OMEGA1 = OMEGA(IC,IA)                       ! see comment in fct header
            XO = OMEGA1*I2
            C(IC,IA) = C0(IC,IA) + C1(IC,IA)*EXP(-XO)   ! see comment in fct header
        ENDDO
      ENDDO

      APHI=.377+4.684D-4*(T-273.15)+3.74D-6*(T-273.15)**2

! --  CALCULATE SUM MI
      XMI=0.
      DO IC=1,NC
         XMI=XMI+MC(IC)
      ENDDO
      DO IA=1,NA
         XMI=XMI+MA(IA)
      ENDDO

      FPHI1 = -APHI*I**(3./2.)/(1.+1.2*I2)
      XS = 0.
      DO IA=1,NA
         DO IC=1,NC
            XS = XS + MA(IA)*MC(IC)*(ZI*C(IC,IA)+BPHI(IC,IA))
         ENDDO
      ENDDO

      ICHARG=1
      JCHARG=2
      CALL EFUNC(ICHARG,JCHARG,APHI,I,E,ED)

      PP=E+I*ED

      F3 = 0.
      DO IC1 = 1, NC
         DO IC2 = IC1+1, NC
            Z1=ZC(IC1)
            Z2=ZC(IC2)
            IF(Z1.NE.Z2) F3=F3 + PP*MC(IC1)*MC(IC2)
         ENDDO
      ENDDO

      F4 = 0.
      DO IA1 = 1, NA
         DO IA2 = IA1+1, NA
            Z1 = ZA(IA1)
            Z2 = ZA(IA2)
            IF(Z1.NE.Z2) F4 = F4 + PP*MA(IA1)*MA(IA2)
        ENDDO
      ENDDO

      PHIX =  FPHI1 + XS  + F3 + F4
      PHI = 1.+ PHIX*2./XMI
      AS = (-PHI)*18./1000.*XMI
      GAMMASN = EXP(AS)

      END FUNCTION GAMMASN
