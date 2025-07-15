! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_CTL_MOD

  CONTAINS

    SUBROUTINE LTDIR_CTL_SEND(SENDID)
      
    USE PARKIND1  ,ONLY : JPIM     ,JPRB
    USE PROGRESS_THREAD
    USE MPI
    
    IMPLICIT NONE
    
    INTEGER(KIND=JPIM),INTENT(IN) :: SENDID
    INTEGER :: STATUS
    
    CALL GSTATS(153,0)
!    CALL TRLTOM_SEND(FOUBUF_IN(IST:IEN),FOUBUF(IST:IEN),2*KF_FS,A2AREQ,ILENS,ILENR,IOFFS,IOFFR)
    print *,'Starting request set ',sendid
    CALL PT_REQSET_START(SENDID,STATUS)
    IF(STATUS .EQ. MPI_ERR_ARG) THEN
       PRINT *,'ERROR STATUS IN LTDIR_CTL_SEND'
    ENDIF
    CALL GSTATS(153,1)
!IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)

  END SUBROUTINE LTDIR_CTL_SEND

    SUBROUTINE LTDIR_CTL(IOFF,KF_FS,KF_UV,KF_SCALARS,  &
 & PSPVOR,PSPDIV,PSPSCALAR, &
 & PSPSC3A,PSPSC3B,PSPSC2, &
 & KFLDPTRUV,KFLDPTRSC)
  
!**** *LTDIR_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL LTDIR_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_FS      - number of fields in Fourier space
!     KF_UV      - local number of spectral u-v fields
!     KF_SCALARS - local number of scalar spectral fields
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
!     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : LALLOPERM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D

USE LTDIR_MOD       ,ONLY : LTDIR
USE TRLTOM_MOD
USE MPI
! USE PROGRESS_THREAD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS,IOFF
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILED2,IST,IEN,IERR

!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

!CALL MPI_WAIT(A2AREQ,MPI_STATUS_IGNORE,IERR)

! Direct Legendre transform

CALL GSTATS(153,0)
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*KF_FS)
CALL GSTATS(153,1)

! if(luse_progress_thread) then
   ! call pause_mpi_helper
! endif

CALL GSTATS(103,0)
CALL GSTATS(1645,0)
IF(KF_FS>0) THEN
   ILED2 = 2*KF_FS
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTDIR(IOFF,IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC)
  ENDDO
!$OMP END PARALLEL DO
ENDIF
CALL GSTATS(1645,1)

!IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)
CALL GSTATS(103,1)

!     -----------------------------------------------------------------

END SUBROUTINE LTDIR_CTL
END MODULE LTDIR_CTL_MOD
