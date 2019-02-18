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


module config

! Description :
! -----------
!   General configuration of the model
!   Declaration of all user-defined parameters
!   Contains utilitary subroutines to read/write these parameters

! Interface :
! ---------

! Input :
! -----

! Output :
! ------

! Externals :
! ---------

! Method :
! ------

! Author :
! ------
!   Josue Bock

! Modifications :
! -------------
  !   12-Mar-2017  Josue Bock  First version based on the initial parameter reading (file istart) which was done in MAIN program
  !
  !   03-Nov-2017  Josue Bock  Defined namelist /mistra_cfg/ that contains the former "istart" file parameters
!-----------------------------------------------------------------------------------------------------------------------------------



! Declarations:
! ------------
! Modules used:
use precision, only : &
  dp                    ! double precision kind


implicit none

public
save

logical :: &
  rst,     & ! rst      : restart or initialization of program run
  netCDF,  & ! netCDF   : output in netCDF format
  binout,  & ! binout   : output in binary format
  mic,     & ! mic      : microphysics included
  chem,    & ! chem     : chemistry included
  halo,    & ! halo     : halogen chemistry on/off
  iod,     & ! iod      : iodine chemistry on/off
  box,     & ! box      : box model run
  BL_box,  & ! BL_box   : box only, average init cond over BL and/or mean of J-values over BL
  nuc,     & ! nuc      : nucleation on/off
  Napari,  & ! Napari   : nuc only, Napari = ternary H2SO4-H2O-NH3 nucleation
  Lovejoy    ! Lovejoy  : nuc only, Lovejoy = homogeneous OIO nucleation

integer :: &
  iaertyp, & ! iaertyp  : type of aerosol; 1=urban, 2=rural, 3=ocean, 4=background
  ifeed,   & ! ifeed    : retroaction over microphysics, and/or chemistry. See manual.
  lstmax,  & ! lstmax   : integration time in hours
  neula,   & ! neula    : eulerian (0) or lagrangian (1) view
  nlevbox, & ! nlevbox  : box only, level to be used for init cond of box if  BL_box=false
  nkc_l      ! nkc_l    : number of output classes for aq. chem.

real (KIND=dp) :: &
  scaleo3_m,      & ! scaleo3_m: total O3 in DU (for photolysis only)
  z_box             ! z_box    : height of MBL (if box run)

character (len=100) :: cnmlfile

character (len=100) :: cinpdir      ! input directory: general data files for Mistra
character (len=109) :: cinpdir_phot ! input directory for photolysis data files
character (len=100) :: coutdir      ! output directory
character (len=100) :: cmechdir     ! mechanism directory

namelist /mistra_cfg/ &
     rst,             &
     lstmax,          &
     netcdf,          &
     binout,          &
     mic,             &
     iaertyp,         &
     chem,            &
     halo,            &
     iod,             &
     nkc_l,           &
     neula,           &
     box,             &
     bl_box,          &
     nlevbox,         &
     z_box,           &
     nuc,             &
     scaleo3_m

contains


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine read_config

! Description :
! -----------
!   Read configuration file


! Interface :
! ---------
!   SR read_config is called by main program during initialisation

! Input :
! -----
!   configuration file 'config.txt'

! Output :
! ------

! Externals :
! ---------
!   none

! Method :
! ------

! Author :
! ------
!   Josue Bock

! Modifications :
! -------------
  !   12-Mar-2017  Josue Bock  First version based on the initial parameter reading (file istart)
  !                            which was done in MAIN program
  !
  !      Nov-2017  Josue Bock  Rewritten (almost) from scratch, with a different structure (env variables, ...)
!--------------------------------------------------------------------------------------------------------------

  use file_unit, only : &
! Imported Parameters:
       jpfunnam,        &
       jpfunout,        &
       jpfuncfgout

  use precision, only : &
! Imported Parameters:
       dp

  implicit none

  integer :: istat

! =======================================================
! -- 1. -- get I/O directories from environment variables
! =======================================================
  call getenv ('INPDIR',cinpdir)
  if (trim(cinpdir) == '') then
     call abortM ('Error during initialisation: no input directory specified as environment variable INPDIR')
  end if
  cinpdir_phot = trim(cinpdir)//'photolys/'

  call getenv ('OUTDIR',coutdir)
  if (trim(coutdir) == '') then
     call abortM ('Error during initialisation: no output directory specified as environment variable OUTDIR')
  end if

  call getenv ('MECHDIR',cmechdir)
  if (trim(cmechdir) == '') then
     call abortM ('Error during initialisation: no mechanism directory specified as environment variable MECHDIR')
  end if

! ===============================================================
! -- 2. -- Mistra_cfg namelist: default values, and read namelist
! ===============================================================
! Default values
rst = .false.
lstmax = 1
netCDF = .false.
binout = .false.
mic = .false.
iaertyp = 3
chem = .false.
halo = .false.
iod = .false.
nkc_l = 4
neula = 1
box = .false.
bl_box = .false.
nlevbox = 2
z_box = 700._dp
nuc = .false.
scaleo3_m = 300._dp

call getenv ('NAMELIST',cnmlfile)

if (trim(cnmlfile) /= '') then
   open (UNIT=jpfunnam, FILE=cnmlfile, STATUS='old', FORM='formatted', IOSTAT=istat)
   if (istat /= 0) call abortM ('Error in SR read_config: cannot open namelist file: '//cnmlfile)

   read (UNIT=jpfunnam, NML=mistra_cfg, IOSTAT=istat)
   if (istat /= 0) call abortM ('Error in SR read_config: cannot read namelist mistra_cfg in file: '//cnmlfile)

   close (UNIT=jpfunnam)
else
   write(jpfunout,'(a)') 'Warning: no namelist specified, only hardcoded default settings will be used'
end if

! =====================================================
! -- 3. -- Export current configuration in file cfg.out
! =====================================================
  open (unit=jpfuncfgout, FILE=trim(coutdir)//'cfg.out', STATUS='new', FORM='formatted', IOSTAT=istat)
  if (istat /= 0) call abortM ('Error in SR read_config: cannot open cfg.out file in dir: '//coutdir)

  write (jpfuncfgout,'(a)') 'Mistra configuration:'
  write (jpfuncfgout,'(a)') '  input directory:',cinpdir
  write (jpfuncfgout,'(a)') '  output directory:',coutdir
  write (jpfuncfgout,'(a)') '  mechanism directory:',cmechdir
  write (jpfuncfgout,'(a)') '  mistra_cfg namelist:'

  write (unit=jpfuncfgout, NML=mistra_cfg, IOSTAT=istat)
  if (istat /= 0) call abortM ('Error in SR read_config: cannot write current mistra_cfg namelist values in cfg.out file.')

  close (unit=jpfuncfgout)

end subroutine read_config
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine abortM (cderrmessage)

  ! Description :
  ! -----------
  !    abortM is called after and error arise somewhere in the code.
  !    The message is written in stderr file, and the program is stopped
  !
  !    Calling abort (Fortran intrinsic function) did not work well with the
  !    current calling param file (probably because of a returned non-zero
  !    STATUS value). This might be imporved in the future.

  ! Author :
  ! ------
  !    Josue Bock


  use file_unit, only : &
       jpfunerr
  implicit none
  character (len=*), intent(in) :: cderrmessage
  write (jpfunerr,'(a)') cderrmessage
  stop '  --> stopped by SR abort'

end subroutine abortM
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module config
