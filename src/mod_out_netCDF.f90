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


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MODULE cdf_var_gas

! contains the old common /cdf_var_gas/
! data exchanged between open / write / close gas subroutines

! Author:
! -------
!     Josue Bock

! History:
! --------
!     20/11/2016 -- Josue Bock -- First version

USE gas_common, ONLY : &
&     j1, &
&     j5

IMPLICIT NONE
SAVE

INTEGER :: idgasfile, igascount
INTEGER, ALLOCATABLE :: idvar_gas(:)

END MODULE cdf_var_gas
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
