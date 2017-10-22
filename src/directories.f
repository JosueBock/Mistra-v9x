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


      module directories

      implicit none

! input directory: general data files for Mistra
      character (len=*) inpdir
      parameter (inpdir='/local/hmb15gqu/mistra/input/')

      character (len=len(inpdir)+9) inpdir_phot       ! input directory for photolysis data files
      parameter ( inpdir_phot = inpdir//'photolys/')


      end module directories
