!
! Copyright 2016-2017 Josue Bock
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


! In this file are all modules which replace old common block

      module gas_common

      USE global_params, ONLY :
     &     n

      implicit none

      save

! GASES (non radical)
      ! Numbers of gas species
      integer :: j1                             ! Number of non radical gases actually used
      integer :: j1_br, j1_cl, j1_iod, j1_halo  ! Number of brominated / chlorinated / ...  gases

      ! Indexes
      integer, allocatable :: ind_gas (:)        ! User defined index for non radical gas species
      integer, allocatable :: ind_gas_rev (:)    ! revert old hard coded indexes. Should disappear soon...
      integer, allocatable :: ind_gas_br  (:,:)  ! Same for halogenated species:
      integer, allocatable :: ind_gas_cl  (:,:)  !   - first column holds stoechiometric nb
!     integer, allocatable :: ind_gas_iod (:,:)  !   - second column holds the (compressed) index

      ! Logical
      logical, allocatable :: gas_is_halo (:)    ! True is halogenated (Br, Cl, I) gas. Used to initialise differently these gases.

      ! Names
      character(len=12),  allocatable :: gas_name (:)
      character(len=100), allocatable :: gas_name_long (:)

      ! Data:
      !   - molar mass of non radical gas [kg/mol]
      double precision, allocatable :: gas_mass (:)
      !   - concentration of non radical gas [mol/m3]
      double precision, allocatable :: s1 (:,:)
      !   - user defined initial concentrations at ground and top level [ppb]
      double precision, allocatable :: s1_init_grd(:)
      double precision, allocatable :: s1_init_top(:)
      !   - ground emission of the gas []
      double precision, allocatable :: es1 (:)
      double precision, allocatable :: vg (:)

! Conversion tables, gas species
      ! KPP to Mistra: in Mistra order
      integer, allocatable :: gas_k2m_g (:)
      integer, allocatable :: gas_k2m_a (:)
      integer, allocatable :: gas_k2m_t (:)
      ! Mistra to KPP: in KPP order
      integer, allocatable :: gas_m2k_g (:,:)
      integer, allocatable :: gas_m2k_a (:,:)
      integer, allocatable :: gas_m2k_t (:,:)



! RADICALS (gas phase)
      ! Numbers of radical species
      integer :: j5                            ! Number of radicals actually used
      integer :: j5_br, j5_cl, j5_iod, j5_halo ! Number of brominated / chlorinated / ...  radicals

      ! Indexes
      integer, allocatable :: ind_rad (:)        ! User defined index for radical species
      integer, allocatable :: ind_rad_br  (:,:)  ! Same for halogenated species:
      integer, allocatable :: ind_rad_cl  (:,:)  !   - first column holds stoechiometric nb
!     integer, allocatable :: ind_rad_iod (:,:)  !   - second column holds the (compressed) index

      ! Logical
!     logical, allocatable :: rad_is_halo (:)

      ! Names
      character(len=12),  allocatable :: rad_name (:)
      character(len=100), allocatable :: rad_name_long (:)

      ! Data:
      !   - molar mass of radical [kg/mol]
      double precision, allocatable :: rad_mass (:)
      !   - concentration of radical [mol/m3]
      double precision, allocatable :: s3 (:,:)

! Conversion tables, raedical species
      ! KPP to Mistra: in Mistra order
      integer, allocatable :: rad_k2m_g (:)
      integer, allocatable :: rad_k2m_a (:)
      integer, allocatable :: rad_k2m_t (:)
      ! Mistra to KPP: in KPP order
      integer, allocatable :: rad_m2k_g (:,:)
      integer, allocatable :: rad_m2k_a (:,:)
      integer, allocatable :: rad_m2k_t (:,:)


! FIXED species which are exchanged
      integer, parameter :: j4 = 1
      double precision :: fix_mass(j4)
      data fix_mass /32.d-3/

      double precision, allocatable :: vmean_init(:)
      double precision, allocatable :: vmean(:,:)

      end module gas_common
