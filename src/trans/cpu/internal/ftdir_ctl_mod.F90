! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_CTL_MOD

USE PARKIND1, ONLY: JPRB

IMPLICIT NONE

CONTAINS

SUBROUTINE FTDIR_CTL_COMP(PIN,PGTF,KF_FS)

USE PARKIND1,        ONLY: JPIM, JPRB
!USE FOURIER_OUT_MOD, ONLY: FOURIER_OUT
USE FTDIR_MOD,       ONLY: FTDIR
USE TPM_DISTR,       ONLY: D

! Dummy arguments

REAL(KIND=JPRB), INTENT(IN) :: PIN(:,:)
REAL(KIND=JPRB), INTENT(INOUT) :: PGTF(:,:)
!REAL(KIND=JPRB), CONTIGUOUS, INTENT(OUT) :: FOUBUF_IN(:)
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS

INTEGER(KIND=JPIM) :: JGL, IBLEN

!     ------------------------------------------------------------------

! Fourier transform
CALL GSTATS(106,0)

CALL GSTATS(1640, 0)
! If this rank has any Fourier fields, Fourier transform them
IF (KF_FS > 0) THEN
  ! Loop over latitudes
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL)
  DO JGL = 1, D%NDGL_FS
    ! Fourier transform
    CALL FTDIR(PIN,PGTF, KF_FS, JGL)
 enddo
!  DO JGL = 1, D%NDGL_FS
!    ! Save Fourier data in FOUBUF_IN
!    CALL FOURIER_OUT(PGTF, FOUBUF_IN,KF_FS, JGL)
!  ENDDO
  !$OMP END PARALLEL DO
ENDIF
CALL GSTATS(1640, 1)

CALL GSTATS(106,1)

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR_CTL_COMP

END MODULE FTDIR_CTL_MOD



