! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! #define DEBUG

MODULE DIR_TRANS_CTL_MOD

USE PARKIND1,          ONLY: JPIM, JPRB
USE OVERLAP_TYPES_MOD, ONLY: BATCH, BATCHLIST, STAT_WAITING, STAT_PENDING, ACTIVE_BATCHES
USE TPM_DISTR       ,ONLY : MYPROC
USE PROGRESS_THREAD
use mpi

IMPLICIT NONE

REAL(KIND=JPRB), ALLOCATABLE :: BGTF(:,:),BIN(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IREQ_RECV(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IREQ_SEND(:,:)
LOGICAL, ALLOCATABLE :: FIRST_PASS(:)
REAL(KIND=JPRB), ALLOCATABLE :: ZCOMBUFR(:,:),ZCOMBUFS(:,:)
INTEGER, PUBLIC :: MAX_COMM(2)
INTEGER, PUBLIC, PARAMETER :: MAX_ACTIVE_BATCHES = 5
INTEGER, PUBLIC, PARAMETER :: STAGE_FINAL = 3
INTEGER :: NCOMM_STARTED(2)
INTEGER :: LAST_SUBMITTED(2)
INTEGER, ALLOCATABLE :: RECV_ID(:),SEND_ID(:,:)
!INTEGER(KIND=JPIM), ALLOCATABLE :: NPTRSPUV(:),NPTRSPSC(:) !NPTRGP(:),
LOGICAL FIRST_TIME

CONTAINS
SUBROUTINE DIR_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)

!**** *DIR_TRANS_CTL* - Control routine for direct spectral transform.

!     Purpose.
!     --------
!        Control routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity
!     PSPDIV(:,:)  - spectral divergence
!     PSPSCALAR(:,:) - spectral scalarvalued fields
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     PGP(:,:,:)  - gridpoint fields

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

USE TPM_GEN,         ONLY: NPROMATR, NOUT
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
USE LINKED_LIST_M,   ONLY: LINKEDLISTNODE
USE TPM_DISTR,       ONLY: D, NPROC, NPRTRNS
USE TPM_TRANS,       ONLY: NGPBLKS,FOUBUF,FOUBUF_IN
USE TRGTOL_MOD,      ONLY: TRGTOL_PROLOG
USE LTDIR_CTL_MOD, ONLY: LTDIR_CTL
USE TRGTOL_MOD,    ONLY: TRGTOL_COMM_RECV
USE FTDIR_CTL_MOD, ONLY: FTDIR_CTL_COMP
USE FOURIER_OUT_MOD, ONLY: FOURIER_OUT
USE PROGRESS_THREAD

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC2(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPUV(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3A(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3B(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP2(:,:,:)

TYPE(BATCH) :: NEW_BATCH
TYPE(BATCH), POINTER :: COMPLETE_COMM_BATCH
INTEGER(KIND=JPIM) :: IBLKS, JBLK
INTEGER(KIND=JPIM) :: NDONE, KNRECV, KRECVCOUNT, KSENDCOUNT, KNSEND
INTEGER(KIND=JPIM), ALLOCATABLE :: KSENDTOT(:), KRECVTOT(:), KSEND(:), KRECV(:), KNDOFF(:),KPTRSPUV(:),KPTRSPSC(:)
INTEGER(KIND=JPIM) :: KGPTRSEND(2,NGPBLKS,NPRTRNS)
INTEGER(KIND=JPIM) :: IVSET(KF_GP)
INTEGER(KIND=JPIM) :: KINDEX(D%NLENGTF)  
TYPE(LINKEDLISTNODE), POINTER :: IB
INTEGER(KIND=JPIM)  :: IOFFSEND, IOFFRECV, IOFFGTF, IOFFGP
INTEGER(KIND=JPIM) :: IST,NACTIVE,IBLEN,IEN,OFFRECV
LOGICAL :: PRODUCTIVE, COMM_COMPL
INTEGER(KIND=JPIM) :: KSENDCOUNT_GLOB
INTEGER I,K,JGL
character(len=25) :: s1,s2,s3,str

!     ------------------------------------------------------------------

! Perform transform

!! WRITE(NOUT,*) "KF_UV_G = ", KF_UV_G
!! WRITE(NOUT,*) "KF_SCALARS_G = ", KF_SCALARS_G
!! WRITE(NOUT,*) "KF_GP = ", KF_GP
!! WRITE(NOUT,*) "KF_FS = ", KF_FS
!! WRITE(NOUT,*) "KF_UV = ", KF_UV
!! WRITE(NOUT,*) "KF_SCALARS = ", KF_SCALARS

IF (NPROMATR > 0) THEN

  ! Determine number of batches
  IBLKS = (KF_GP - 1) / NPROMATR + 1

  IF(.NOT. ALLOCATED(RECV_ID)) THEN
     ALLOCATE(RECV_ID(IBLKS))
     ALLOCATE(SEND_ID(2,IBLKS))
  ENDIF

  ! ================================================================================================
  ! Initialise work arrays shared by all batches
  ! ================================================================================================

  ! Allocate grid-to-Fourier transform buffer
  IF (ALLOCATED(BIN) .AND. SIZE(BIN, 2) /= KF_FS .AND. SIZE(BIN, 1) /= (D%NLENGTF)) THEN
    DEALLOCATE(BIN,BGTF)
  ENDIF
  IF (.NOT. ALLOCATED(BIN)) THEN
    ALLOCATE(BGTF(D%NLENGTF,KF_FS))
    ALLOCATE(BIN(D%NLENGTF,KF_FS))
  ENDIF
  ! Now, force the OS to allocate these arrays right now
  BGTF(1,1)=HUGE(1._JPRB)
  BIN(1,1)=HUGE(1._JPRB)

  ALLOCATE(KSENDTOT(NPROC))
  ALLOCATE(KRECVTOT(NPROC))
  ALLOCATE(KSEND(NPROC))
  ALLOCATE(KRECV(NPROC))
  ALLOCATE(KNDOFF(NPROC))

  ! Create a combined V-set array
  IST = 1
  IF (KF_UV_G > 0) THEN
    IVSET(IST:IST+KF_UV_G-1) = KVSETUV(:)
    IST = IST+KF_UV_G
    IVSET(IST:IST+KF_UV_G-1) = KVSETUV(:)
    IST = IST+KF_UV_G
  ENDIF
  IF (KF_SCALARS_G > 0) THEN
    IVSET(IST:IST+KF_SCALARS_G-1) = KVSETSC(:)
    IST = IST+KF_SCALARS_G
  ENDIF


    ! Call TRGTOL_PROLOG on "global" parameters to determine sizes of communication buffers
  CALL TRGTOL_PROLOG(KF_FS, KF_GP, IVSET, KSENDCOUNT, KRECVCOUNT, KNSEND, KNRECV, KSENDTOT, &
    &                KRECVTOT, KSEND, KRECV, KINDEX, KNDOFF, KGPTRSEND)

  ! Allocate receive request handle array
  KSENDCOUNT_GLOB = SUM(KSENDTOT)

  IF (ALLOCATED(IREQ_RECV) .AND. SIZE(IREQ_RECV,1) /= KNRECV .AND. SIZE(IREQ_RECV,2) /= IBLKS) THEN
    DEALLOCATE(IREQ_RECV)
  ENDIF
  IF (.NOT. ALLOCATED(IREQ_RECV)) THEN
    ALLOCATE(IREQ_RECV(KNRECV,IBLKS))
  ENDIF

  IF (ALLOCATED(IREQ_SEND) .AND. SIZE(IREQ_SEND,1) /= KNSEND .AND. SIZE(IREQ_SEND,2) /= IBLKS) THEN
     DEALLOCATE(IREQ_SEND)
  ENDIF
  IF (.NOT. ALLOCATED(IREQ_SEND)) THEN
    ALLOCATE(IREQ_SEND(KNSEND,IBLKS))
  ENDIF

  ! Allocate communication buffers
  IF (ALLOCATED(ZCOMBUFR) .AND. SIZE(ZCOMBUFR,2) /= KNRECV .AND. SIZE(ZCOMBUFR,1) /= KRECVCOUNT) THEN
    DEALLOCATE(ZCOMBUFR)
  ENDIF
  IF (.NOT. ALLOCATED(ZCOMBUFR)) THEN
    ALLOCATE(ZCOMBUFR(KRECVCOUNT,KNRECV))
  ENDIF
  IF (ALLOCATED(ZCOMBUFS) .AND. SIZE(ZCOMBUFS,2) /= KNSEND .AND. SIZE(ZCOMBUFS,1) /= KSENDCOUNT_GLOB) THEN
    DEALLOCATE(ZCOMBUFS)
  ENDIF
  IF (.NOT. ALLOCATED(ZCOMBUFS)) THEN
    ALLOCATE(ZCOMBUFS(KSENDCOUNT_GLOB,KNSEND))
  ENDIF

  IF(.NOT. ALLOCATED(FIRST_PASS)) THEN
     ALLOCATE(FIRST_PASS(IBLKS))
     DO I=1,IBLKS
        FIRST_PASS(I) = .TRUE.
     ENDDO
     first_time = .true.
  ENDIF
  
  IBLEN = D%NLENGT0B*2*KF_FS
  IF (ALLOCATED(FOUBUF)) THEN
    IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
      DEALLOCATE(FOUBUF)
      ALLOCATE(FOUBUF(MAX(1,IBLEN)))
    ENDIF
  ELSE
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
  IF (ALLOCATED(FOUBUF_IN)) THEN
    IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
      DEALLOCATE(FOUBUF_IN)
      ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
    ENDIF
  ELSE
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  ENDIF

  ! ================================================================================================
  ! Begin overlap loop
  ! ================================================================================================

  IOFFSEND = 1
  IOFFRECV = 1
  IOFFGTF = 1
  IOFFGP = 1

  LAST_SUBMITTED(1) = 0
  LAST_SUBMITTED(2) = 0

  MAX_COMM(1) = 1
  MAX_COMM(2) = 1
  JBLK = 1 ! This keeps track of the last activated batch
  NDONE = 0 ! This keeps track of the number of completed batches
  NCOMM_STARTED(1) = 0 ! This keeps track of the batches in an active communication
  NCOMM_STARTED(2) = 0 ! This keeps track of the batches in an active communication
  NACTIVE = 1 ! This keeps track of the overall number of active batches
  
!  print *,'krecvcount,krecvtot,kf_fs=',krecvcount,krecvtot,kf_fs
!  call gstats(901,0)
  CALL ACTIVATE(1, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, PGP, IOFFSEND, IOFFRECV, &
    &           IOFFGTF, IOFFGP, KSENDCOUNT, KRECVCOUNT)
  IB => ACTIVE_BATCHES%HEAD

    DO JBLK =2,IBLKS+1

!    if(.not. FIRST_PASS(JBLK)) THEN
!       call gstats(902,0)
       !    ENDIF
       IF(JBLK .LE. IBLKS) THEN
          CALL ACTIVATE(JBLK, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, PGP, IOFFSEND, IOFFRECV, &
               &           IOFFGTF, IOFFGP, KSENDCOUNT, KRECVCOUNT)
       ENDIF
!    if(.not. FIRST_PASS(JBLK)) THEN
!       call gstats(902,1)
!    ENDIF
!    print *,'Processing batch ',jblk,' ioffsend=',ioffsend

  ! Check whether any active batches have a completed communication
   SELECT TYPE (THISBATCH => IB%VALUE)
    TYPE IS (BATCH)
!       do while (.not. THISBATCH%COMM_COMPLETE(IREQ_RECV(:,THISBATCH%NBLK)))

!       enddo
    if(luse_progress_thread) then

!       print *,'waiting for request ',THISBATCH%RECV_ID(1)
!       if(.not. FIRST_PASS(JBLK)) THEN
          call gstats(903,0)
!       ENDIF
       CALL PT_REQSET_WAIT(THISBATCH%RECV_ID(1))
!       if(.not. FIRST_PASS(JBLK)) THEN
          call gstats(903,1)
!       endif
    endif
    call gstats(906,0)

    IST = THISBATCH%IOFFGTF
    IEN = THISBATCH%IOFFGTF+THISBATCH%NF_FS-1
    OFFRECV = THISBATCH%MYOFFRECV
    CALL TRGTOL_COMM_RECV(BIN(:,IST:IEN), ZCOMBUFR, &
         & OFFRECV,THISBATCH%NF_FS, THISBATCH%NRECVCOUNT, THISBATCH%NNSEND,THISBATCH%NNRECV, &
         & THISBATCH%NRECVTOT, THISBATCH%NRECV, &
     &                 THISBATCH%NINDEX, THISBATCH%NNDOFF,THISBATCH%IREQ_SEND,THISBATCH%IREQ_RECV)
    
     call gstats(906,1)

     IST = 1+D%NLENGT0B*2*(THISBATCH%IOFFGTF-1)
     IEN = IST + D%NLENGT0B*2*THISBATCH%NF_FS-1

     CALL FTDIR_CTL_COMP(BIN(:,THISBATCH%IOFFGTF:THISBATCH%IOFFGTF+THISBATCH%NF_FS-1), &
          &  BGTF(:,THISBATCH%IOFFGTF:THISBATCH%IOFFGTF+THISBATCH%NF_FS-1), THISBATCH%NF_FS)

  !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JGL)
     DO JGL = 1, D%NDGL_FS
    ! Save Fourier data in FOUBUF_IN
        CALL FOURIER_OUT(BGTF(:,THISBATCH%IOFFGTF:THISBATCH%IOFFGTF+THISBATCH%NF_FS-1), FOUBUF_IN,KF_FS, &
             THISBATCH%NF_FS,JGL,THISBATCH%IOFFGTF)
     ENDDO
  !$OMP END PARALLEL DO

 if(luse_progress_thread) then
       !  call gstats(904,0)
       CALL PT_REQSET_WAIT(THISBATCH%SEND_ID(1))
       !  call gstats(904,1)
    endif
       !   call ACTIVE_BATCHES%REMOVE(IB)
     IB => IB%NEXT
     
  END SELECT

!    CALL ACTIVE_BATCHES%REMOVE(IB)
    

ENDDO
ELSE

  ! No splitting of fields, transform done in one go
  CALL ABORT_TRANS("DIR_TRANS_CTL: NPROMATR = 0 feature disabled for overlap version")

ENDIF



!call gstats(901,1)

!    CALL TRGTOL_COMM_RECV(BIN, ZCOMBUFR, 1, &
!         &                 KF_FS, KRECVCOUNT, KNRECV, KRECVTOT, KRECV, &
!     &                 KINDEX, KNDOFF)

!CALL FTDIR_CTL_COMP(BIN,BGTF,FOUBUF_IN,KF_FS)


 
     IF(FIRST_TIME) THEN

#ifdef DEBUG
s1 = 'bin'
write(s2,9) iblks,jblk,myproc-1
9 format(i0,'.',i0,'.',i0)
s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')
do k=1,kf_fs
   do i=1,D%NLENGTF
      write(11,10) i,k,bin(i,k)
   enddo
enddo

10 format(i8,i8,F22.8)
   close(11)
#endif

     
#ifdef DEBUG
s1 = 'bgtf'
write(s2,9) iblks,jblk,myproc-1
s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')
do i=1,D%NLENGTF
   do k=1,kf_fs
      write(11,10) k,i,bgtf(k,i)
   enddo
enddo

   close(11)

s1 = 'foubuf_in'
write(s2,9) iblks,jblk,myproc-1
s3 = trim(s2)
str = trim(s1) // s3
open(11,file=str,form='formatted',status='unknown',action='write')
do k=1,kf_fs * D%NLENGT0B*2
   write(11,12) k,foubuf_in(k)
   enddo

12 format(i8,F22.8)
   close(11)
#endif

endif

FIRST_TIME = .FALSE.

CALL LTDIR_CTL(1, KF_FS, KF_UV, KF_SCALARS, &
         &     PSPVOR=PSPVOR, PSPDIV=PSPDIV, PSPSCALAR=PSPSCALAR)

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_CTL

SUBROUTINE ACTIVATE(N, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, PGP, IOFFSEND, IOFFRECV, &
  &                 IOFFGTF, IOFFGP, SENDCNTMAX, RECVCNTMAX)

  INTEGER,                      INTENT(IN)    :: N
  INTEGER(KIND=JPIM),           INTENT(IN)    :: KF_GP
  INTEGER(KIND=JPIM),           INTENT(IN)    :: KF_SCALARS_G
  INTEGER(KIND=JPIM),           INTENT(IN)    :: KF_UV_G
  INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KVSETUV(:)
  INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)    :: KVSETSC(:)
  REAL(KIND=JPRB),    OPTIONAL, INTENT(IN)    :: PGP(:,:,:)
  INTEGER(KIND=JPIM),           INTENT(INOUT) :: IOFFSEND
  INTEGER(KIND=JPIM),           INTENT(INOUT) :: IOFFRECV
  INTEGER(KIND=JPIM),           INTENT(INOUT) :: IOFFGTF, IOFFGP
  INTEGER(KIND=JPIM),           INTENT(IN)    :: SENDCNTMAX
  INTEGER(KIND=JPIM),           INTENT(IN)    :: RECVCNTMAX
!  INTEGER(KIND=JPIM),           INTENT(INOUT)    ::  NPTRSPUV(:),NPTRSPSC(:) !NPTRGP(:),
  INTEGER :: STATUS
  CLASS(BATCH), POINTER :: NEW_BATCH

  ! Add a new batch to the list
  call gstats(902,0)
  CALL ACTIVE_BATCHES%APPEND(BATCH(N, KF_GP, KF_SCALARS_G, KF_UV_G, KVSETUV, KVSETSC, IOFFSEND, &
       &                     IOFFRECV, IOFFGTF, IOFFGP, SENDCNTMAX, RECVCNTMAX,  &
       &                     IREQ_SEND(:,N),IREQ_RECV(:,N)))
  call gstats(902,1)

     SELECT TYPE (NEW_BATCH => ACTIVE_BATCHES%TAIL%VALUE)
  TYPE IS (BATCH)

     IF(FIRST_PASS(N)) THEN
        CALL NEW_BATCH%INIT_RECVS(ZCOMBUFR,RECV_ID(N))
        CALL NEW_BATCH%INIT_SENDS(ZCOMBUFS,SEND_ID(:,N))
        !        IF(STATUS .EQ. MPI_ERR_ARG) THEN
        !           PRINT *,'ERROR IN STATUS, SENDS IN ACTIVATE'
        !        ENDIF
        FIRST_PASS(N) = .FALSE.
     ENDIF
     if(luse_progress_thread) then
        NEW_BATCH%RECV_ID(1) = RECV_ID(N)
        NEW_BATCH%RECV_ID(2) = SEND_ID(2,N)
        NEW_BATCH%SEND_ID(1) = SEND_ID(1,N)
        NEW_BATCH%SEND_ID(2) = SEND_ID(2,N)
!        print *,'Starting receive request set ',NEW_BATCH%RECV_ID(1)
        call gstats(904,0)
        CALL PT_REQSET_START(NEW_BATCH%RECV_ID(1),STATUS)
        call gstats(904,1)
        IF(STATUS .EQ. MPI_ERR_ARG) THEN
           PRINT *,'ERROR IN STATUS, RECVS IN ACTIVATE'
        ENDIF
     endif
!     IF (NCOMM_STARTED(1) < MAX_COMM(1)) THEN
     call gstats(905,0)
     CALL NEW_BATCH%START_COMM(PGP,ZCOMBUFS,ZCOMBUFR,BIN)
     call gstats(905,1)
     !      NCOMM_STARTED(1) = NCOMM_STARTED(1) + 1
!      LAST_SUBMITTED(1) = LAST_SUBMITTED(1) + 1
!    ENDIF
  END SELECT

END SUBROUTINE ACTIVATE

END MODULE DIR_TRANS_CTL_MOD
