!
! Copyright 2015-2017 Josue Bock
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


      module global_params

      implicit none

      save
! param_met
      ! number of layers, constant height
      integer nf
      parameter ( nf = 100 )

      ! number of layers, constant height + exponentially increasing height
      integer n
      parameter ( n = nf+50 )

      integer nm
      parameter ( nm = n-1 )

! param_mic
      integer nka
      parameter ( nka = 70 )

      integer nkt
      parameter ( nkt = 70 )

! param_rad
      ! number of layers including the extra ones for radiation calculation
      !   see SR initr (in radinit.f) for height definition
      !   note that the infinitesimally thin layer is ignored, thus only
      !   (n-1) layers are taken from the main model to the radiative code
      integer, parameter :: nrlay = (n -1) + 11

      integer, parameter :: nrlev = nrlay + 1


      integer mb
      parameter ( mb = 18 )   ! Total number of wavelength bands in the radiative code

      integer mbs
      parameter ( mbs = 6 )   ! Wavelength bands for the solar region

      integer mbir
      parameter ( mbir = 12 ) ! Wavelength bands for the infrared region

      integer ncw
      parameter ( ncw = 8 )

! parameters soil grid
      integer nb
      parameter ( nb = 20 )

      integer nbm
      parameter ( nbm = nb-1 )

! param_che

!      integer, parameter :: j1 = 96
      integer, parameter :: j1_fake = 96

      integer, parameter :: j3 = 25

      integer, parameter :: j2 = j1_fake + j3

      integer, parameter :: j6 = 55

      integer, parameter :: nkc = 4

      integer, parameter :: nphrxn = 47

      integer, parameter :: nlev = 15

      integer, parameter :: nrxn = 1627

      integer, parameter :: nmax_chem_aer = nf

      end module global_params
