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


module precision

! Define kind of variables used throughout the model


! Author:
! -------
!     Josue Bock


! Modifications :
! -------------
  ! 28-Oct-2017  Josue Bock  introduce tiny_dp to replace == 0. tests by < tiny_dp

implicit none

save

integer, parameter :: dp = kind(1.d0)
real(kind=dp), parameter :: tiny_dp = tiny(0._dp)

end module precision
