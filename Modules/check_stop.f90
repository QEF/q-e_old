!
! Copyright (C) 2002-2004 FPMD-PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! ... This module contains functions to check if the code should
! ... be stopped smootly.
! ... In particular th function check_stop_now return .TRUE. if
! ... either the user has created a given file or if the 
! ... elapsed time is larger than max_seconds
!
!------------------------------------------------------------------------------!
MODULE check_stop
!------------------------------------------------------------------------------!
  !  
  USE kinds
  !
  IMPLICIT NONE
  !
  SAVE
  !
  REAL(dbl)        :: max_seconds = 1.D+7
  LOGICAL, PRIVATE :: tinit = .FALSE.
  !
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_stop_init( val )
       !-----------------------------------------------------------------------
       !
       USE io_global, ONLY : ionode, ionode_id, stdout
       USE io_files,  ONLY : exit_file, stopunit
       !
       IMPLICIT NONE
       !
       REAL(dbl), INTENT(IN) :: val
       LOGICAL               :: tex
       REAL(dbl)             :: seconds
       REAL(dbl)             :: elapsed_seconds
       EXTERNAL                 elapsed_seconds
       !
       IF( tinit ) &
         WRITE( UNIT = stdout, &
                FMT = '("WARNING: check_stop already initialized *** ")' )
       !
       IF ( val > 0.D0 ) max_seconds = val
       !
       IF ( ionode ) THEN
          !
          INQUIRE( FILE = exit_file, EXIST = tex )
          !
          IF ( tex ) THEN
             !
             OPEN( stopunit, FILE = exit_file, STATUS = 'OLD' )
             CLOSE( stopunit, STATUS = 'DELETE' )
             !
          END IF
          !
       END IF
       !
       seconds = elapsed_seconds()
       tinit   = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION check_stop_now()
       !-----------------------------------------------------------------------
       !
       USE mp,        ONLY : mp_bcast
       USE mp_global, ONLY : intra_image_comm 
       USE io_global, ONLY : ionode, ionode_id, stdout
       USE io_files,  ONLY : exit_file, stopunit, iunexit
       !
       IMPLICIT NONE
       !
       LOGICAL   :: check_stop_now, tex
       REAL(dbl) :: seconds
       REAL(dbl) :: elapsed_seconds
       EXTERNAL     elapsed_seconds
       !
       !
       ! ... elapsed_seconds is a C function returning the elapsed solar 
       ! ... time in seconds since the first call to the function itself
       !
       IF( .NOT. tinit ) &
          CALL errore( 'check_stop_now', 'check_stop not initialized', 1 )
       !
       check_stop_now = .FALSE.
       !  
       IF ( ionode ) THEN
          !
          INQUIRE( FILE = TRIM( exit_file ), EXIST = tex )
          !
          IF ( tex ) THEN
             !     
             check_stop_now = .TRUE.
             !
             WRITE( UNIT = stdout, &
                    FMT = '(" *** Program stopped by user request *** ")' )
             !
             OPEN( UNIT = iunexit, FILE = TRIM( exit_file ) )
             CLOSE( UNIT = iunexit, STATUS = 'DELETE' )
             !
          END IF
          !
          seconds = elapsed_seconds()
          !
          IF ( seconds  >  max_seconds ) THEN
             !
             check_stop_now = .TRUE.
             !
             WRITE( UNIT = stdout, &
                    FMT = '(" *** Maximum CPU time exceeded ***")' )
             WRITE( UNIT = stdout, &
                    FMT = '(" *** max_seconds     = ",D10.2," ***")' ) max_seconds
             WRITE( UNIT = stdout, &
                    FMT = '(" *** elapsed seconds = ",D10.2," ***")' ) seconds
             !
          END IF
          !
       END IF
       !
       CALL mp_bcast( check_stop_now, ionode_id, intra_image_comm )
       !
       RETURN
       !
     END FUNCTION check_stop_now
     !
!------------------------------------------------------------------------------!
END MODULE check_stop
!------------------------------------------------------------------------------!
