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

! Declare and initialise all physical and chemical constants used in Mistra


! Author:
! -------
!     Josue Bock

USE precision, only: &
  &   dp                   ! kind double precision real

implicit none

save

! Avogadro constant
real(kind=dp), parameter :: Avogadro = 6.022140857e+23_dp        ! [1/mol]

! Conversion factor 1:
!     conv1 = Avogadro / 10^6
!     multiply by conv1 to convert cm^3(air)/mlc --> m^3(air)/mol
real(kind=dp), parameter :: conv1 = Avogadro * 1.e-6_dp          ! [m3/cm3/mol]

! Air molar mass
real(kind=dp), parameter :: M_air = 28.96546e-3_dp               ! [kg/mol]

! Water molar mass
real(kind=dp), parameter :: M_wat = 18.01528e-3_dp               ! [kg/mol]

! Pi
real(kind=dp), parameter :: pi = 3.1415926535897932_dp

! Molar gas constant
!     ref: Gavioso et al. (2015), Metrologia 52 (S274-S304) doi: 10.1088/0026-1394/52/5/S274
real(kind=dp), parameter :: gas_const = 8.3144743_dp             ! [J/K/mol]

! Specific gas constant of dry air
!     r0 = R / M_air
real(kind=dp), parameter :: r0 = gas_const / M_air               ! [J/(kg.K)]

! Specific gas constant of water vapour
!     r1 = R / M_water
real(kind=dp), parameter :: r1 = gas_const / M_wat               ! [J/(kg.K)]

! Water density
!     This could be improved by a parameterisation as a function of T and P
real(kind=dp), parameter :: rhow = 1000.0_dp                     ! [kg/m3]

! Gravity
real(kind=dp), parameter :: g = 9.80665_dp                       ! [m/s2]

! Specific heat capacity of dry air at constant pressure
!     ref: Seinfeld & Pandis, 2nd Ed., Table A.7 p. 1178
real(kind=dp), parameter :: cp = 1005.0_dp                       ! [J/kg/K]



end module constants
