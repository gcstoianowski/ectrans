! (C) Copyright 1987- ECMWF.
! (C) Copyright 1987- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_MOD
CONTAINS
SUBROUTINE LTDIR(IOFF,KM,KMLOC,KF_FS,KF_UV,KF_SCALARS,KLED2,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC)


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
USE TPM_TRANS, ONLY : LATLON
USE TPM_FLT         ,ONLY : S
USE TPM_GEOMETRY    ,ONLY : G

USE PREPSNM_MOD     ,ONLY : PREPSNM
USE PRFI2_MOD       ,ONLY : PRFI2
USE LDFOU2_MOD      ,ONLY : LDFOU2
USE LEDIR_MOD       ,ONLY : LEDIR
USE UVTVD_MOD       ,ONLY : UVTVD      
USE UPDSP_MOD       ,ONLY : UPDSP
USE CDMAP_MOD , ONLY : CDMAP
USE LINKED_LIST_M, ONLY: LINKEDLISTNODE
USE OVERLAP_TYPES_MOD, ONLY: BATCH, BATCHLIST,ACTIVE_BATCHES

!**** *LTDIR* - Control of Direct Legendre transform step

!     Purpose.
!     --------
!        Tranform from Fourier space to spectral space, compute
!        vorticity and divergence.

!**   Interface.
!     ----------
!        *CALL* *LTDIR(...)*

!        Explicit arguments :
!        --------------------  KM     - zonal wavenumber
!                              KMLOC  - local zonal wavenumber

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI2   - prepares the Fourier work arrays for model variables.
!         LDFOU2  - computations in Fourier space
!         LEDIR   - direct Legendre transform
!         UVTVD   -
!         UPDSP   - updating of spectral arrays (fields)

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-24
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
!        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified 94-04-06 R. El khatib Full-POS implementation
!        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
!                             instead of u,v->vor,div
!        MPP Group : 95-10-01 Support for Distributed Memory version
!        K. YESSAD (AUGUST 1996):
!               - Legendre transforms for transmission coefficients.
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        R. El Khatib 12-Jul-2012 LDSPC2 replaced by UVTVD
!     ------------------------------------------------------------------

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS,KLED2,IOFF

REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)

TYPE(LINKEDLISTNODE), POINTER :: IB
TYPE(BATCH), POINTER :: THISBATCH

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IFC, IIFC, IDGLU
INTEGER(KIND=JPIM) :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE
INTEGER(KIND=JPIM) :: ISL, ISLO, MYOFF

!     LOCAL REALS
!REAL(KIND=JPRB) :: ZSIA(KLED2,R%NDGNH),       ZAIA(KLED2,R%NDGNH)
REAL(KIND=JPRB) :: ZEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB) :: ZOA1(KLED2,R%NLED4),         ZOA2(MAX(4*KF_UV,1),R%NLED4)
REAL(KIND=JPRB), ALLOCATABLE :: ZAIA(:,:), ZSIA(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       4.    DIRECT LEGENDRE TRANSFORM.
!              --------------------------
IFC = 2*KF_FS
IIFC = IFC
IF(KM == 0)THEN
  IIFC = IFC/2
ENDIF

IDGLU = MIN(R%NDGNH,G%NDGLU(KM))
ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

ALLOCATE(ZSIA(KLED2,R%NDGNH))
ALLOCATE(ZAIA(KLED2,R%NDGNH))

IF( LATLON.AND.S%LDLL ) THEN

  IF( (S%LSHIFTLL .AND. KM < 2*IDGLU) .OR.&
   & (.NOT.S%LSHIFTLL .AND. KM < 2*(IDGLU-1)) ) THEN

    CALL PREPSNM(KM,KMLOC,ZEPSNM)
    ISLO = S%FA(KMLOC)%ISLD
    ! map from external to internal (gg) roots and split into anti-symmetric / symmetric
    CALL CDMAP(KM,KMLOC,ISL,ISLO,ZEPSNM(R%NTMAX+1),1_JPIM,&
     & R%NDGNH,S%NDGNHD,IFC,ZAIA,ZSIA)

  ENDIF

ELSE

  CALL PRFI2(IOFF,KM,KMLOC,KF_FS,ZAIA,ZSIA)

ENDIF

CALL LDFOU2(KM,KF_UV,ZAIA,ZSIA)

CALL LEDIR(KM,KMLOC,IFC,IIFC,ISL,IDGLU,KLED2,ZAIA,ZSIA,ZOA1,F%RW(1:R%NDGNH))

DEALLOCATE(ZAIA)
DEALLOCATE(ZSIA)

!     ------------------------------------------------------------------

!*       5.    COMPUTE VORTICITY AND DIVERGENCE.
!              ---------------------------------

IF( KF_UV > 0 ) THEN
  CALL PREPSNM(KM,KMLOC,ZEPSNM)

  IB => ACTIVE_BATCHES%HEAD
  DO WHILE (ASSOCIATED(IB))
     SELECT TYPE (THISBATCH => IB%VALUE)
     TYPE IS (BATCH)
        MYOFF = THISBATCH%IOFFGTF * 2 - 1
        IUS = MYOFF
        IUE = IUS + 2 * THISBATCH%NF_UV -1
        IVS = IUE +1
        IVE = IUS + 4 * THISBATCH%NF_UV -1
        IVORS = IUS
        IDIVS = IVS
        IVORE = IUE
        IDIVE = IVE
        CALL UVTVD(KM,THISBATCH%NF_UV,ZEPSNM,ZOA1(IUS:IUE,:),ZOA1(IVS:IVE,:),&
             & ZOA2(IVORS:IVORE,:),ZOA2(IDIVS:IDIVE,:))

!     ------------------------------------------------------------------

!*       6.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

        CALL UPDSP(KM,THISBATCH%NF_UV,THISBATCH%NF_SCALARS,ZOA1(MYOFF:,:), &
             & ZOA2(MYOFF:,:),PSPVOR,PSPDIV,PSPSCALAR,&
             & PSPSC3A,PSPSC3B,PSPSC2 , &
             & THISBATCH%NFLDPTRUV,THISBATCH%NFLDPTRSC)

        IB => IB%NEXT
     END SELECT
  ENDDO
  
ENDIF


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',1,ZHOOK_HANDLE)

END SUBROUTINE LTDIR
END MODULE LTDIR_MOD
