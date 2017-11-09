module file_unit

! Description :
! -----------
!   Declare and define the unit of each I/O file used in Mistra


! Author :
! ------
!   Josue Bock


  ! JPFUNERR: unit number for error messages
  ! JPFUNOUT: unit number for standard output

  ! JPFUNNAM: unit number for namelist
  ! JPFUNCFGOUT: unit number for config output

  ! JPFUNAERRAD: unit number for aerosol (particles) radiation parameters
  implicit none

  public
  save

  integer, parameter :: jpfunerr = 0
  integer, parameter :: jpfunout = 6

  integer, parameter :: jpfunnam = 1
  integer, parameter :: jpfuncfgout = 2

  integer, parameter :: jpfunaerrad = 51

end module file_unit
