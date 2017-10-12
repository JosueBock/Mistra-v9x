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


      module constants

      implicit none

      ! Avogadro constant
      double precision Avogadro
      parameter ( Avogadro = 6.02214d+23 ) ! [mol-1]

      ! Conversion factor 1:
      !   conv1 = Avogadro / 10^6
      !   multiply by conv1 to convert cm^3(air)/mlc --> m^3(air)/mol
      double precision conv1
      parameter ( conv1 = Avogadro * 1.d-6 )

      ! Air molar mass
      double precision m_air
      parameter ( m_air = 28.97d-3 ) ! [kg/mol]

      ! Pi
      double precision pi
      parameter ( pi = 3.1415926535897932d0 )

      ! Molar gas constant
      !   ref: Gavioso et al. (2015), Metrologia 52 (S274-S304) doi: 10.1088/0026-1394/52/5/S274
      double precision gas_const
      parameter ( gas_const = 8.3144743 ) ! [J/K/mol]

! Water density                             ! This could be improved by a parameterisation as a function of T and P
      double precision water_density
      parameter ( water_density = 1000.d0 ) ! [kg/m**3]

      ! Gravity
      double precision g ! [m/s**2]
      parameter ( g = 9.806 )

      ! Specific heat of dry air
      !    ref: Seinfeld & Pandis, 2nd Ed., Table A.7 p. 1178
      double precision, parameter :: cp = 1005.d0  ! [J/(kg*K)]


      end module constants
