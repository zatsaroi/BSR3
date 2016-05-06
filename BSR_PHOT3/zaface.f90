!======================================================================
      Subroutine ZAFACE (IBUG,AC,ION,km,R,DRI,nch,nopen,LCH,ECH,CF, &
                         iauto, mfgi, info, F,G,FP,GP )
!======================================================================
!      INTERFACE TO CREES ASSYMPTOTIC PROGRAM
!----------------------------------------------------------------------
!      INPUT: 
!
!      IBUG    = DEBUG PRINT
!      AC      = ACCURACY
!      ION     = RESIDUAL CHARGE ON TARGET (Z-NELC)
!      km      = MAX.ORDER OF LONG RANGE POTENTIAL CF
!      R       = BOUNDARY RADIUS
!      DRI     = FOX-GOODWIN STEP
!      NCH     = TOTAL NUMBER OF CHANNELS
!      NOPEN   = NUMBER OF OPEN CHANNELS
!      CF      = LONG RANGE POTENTIAL COEFFICIENTS, UPPER TRIANGLE FORM
!      LCH     = CHANNEL ANGULAR MOMENTA
!      ECH     = CHANNEL ENERGY (RYDBERS)
!      iauto   = mode for automatic choose og MFG
!      mfgi    = suggested number of points in outer region (MFG)
!
!      OUTPUT:
!
!      F,G     = sin/cos solutions ON BOUNDARY
!      FP,GP   = their derivatives
!
!      info    = if not equal 0 - skip the energy!
!-----------------------------------------------------------------------
      IMPLICIT REAL(8) (A-H,O-Z)
 
! ... here the parameters from ASYPCK:

      PARAMETER (mchf=200, MLMX = 5, MXMFG=540)
 
      DIMENSION LCH(nch),CF(nch,nch,km),ECH(nch), &
                F(nch,nch), G(nch,nch), FP(nch,nch), GP(nch,nch)
 
      DIMENSION FA(MCHF,MCHF,2), FAD(MCHF,MCHF,2)
 
! ... CREES COMMON BLOCK ...
 
      COMMON /SAVE/ EL(MCHF), EN(MCHF), BLAM(MCHF,MCHF,MLMX),          &
             F1(2,MCHF,MCHF),F2(2,MCHF,MCHF),B1,B2,ECON,EROR,SQK,DELNU,&
             DEL1,SZ,R1,R2,R0,RA,RB,DR,KJ(MCHF),METHOD,MFG,INTER,METH, &
             ISOLN,KMETH,IPUNC,NCHF,NCHOP,NCHCL,LAMMX,MT6,NBUG,IFILE,  &
             JFILE,KFILE,MFGP1,IALT,JALT,KALT
       
      if(nch.gt.MCHF) Stop ' ZAFACE:  nch > MCHF'
      if(mfgi+3.gt.MXMFG) Stop ' ZAFACE:  mfgi > MXMFG'
 
      info = 0; if(nopen.le.0) then;  info = 1; Return; end if
      info = 0; if(nch.le.0) then;  info = 2; Return; end if

      MT6 = 6;         ! print unit
 
      MFG = mfgi       ! NUMBER OF STEPS WITH FOX-GOODWIN
 
! ... IAUTO = 0 FOR USING SPECIFIED VALUE OF MFG,
! ...       = 1 FOR CONDITIONAL AUTOMATIC INCREASE MFG  
! ...       = 2 FOR UNCONDITIONAL AUTOMATIC INGREASE MFG
 
! ... NBUG  = DEBUG PARAMETER (=0 FOR NO L/P OUTPUT FROM ASYPCK)
 
      NBUG = IBUG    !  =0,1,2
 
! ... EROR  = ACCURACY PARAMETER FOR ASYMPTOTIC SOLUTIONS
 
      EROR = AC         !  from 0.01 to 0.0001
 
      R0 = R            !  radial point for soulutions

! ... DR   = STEP LENGTH FOR FOX-GOODWIN INTEGRATION; < 0.2
 
      DR = DRI
 
! ... first, solutions are obtained in the R1=R-DR and R2=R+DR,
! ... then by interpolation - in R. So, for interpolation it is 
! ... better small DR, but for FOX-GOODWIN - larger DR:
! ... solution by expantion will be seeked in R + DR*MFG
!
! ... Begin with 0.1 ! 

      IF (IBUG.GT.0) THEN
        WRITE (MT6,'(/a/)') &
        ' * ASYMPTOTIC PROGRAM FROM MARTIN CREES, 1980 *'
        WRITE (MT6,'(a,i3)') ' MFG   =',MFG
        WRITE (MT6,'(a,i3)') ' IAUTO =',IAUTO
        WRITE (MT6,'(a,i3)') ' NBUG  =',NBUG
        WRITE (MT6,'(a,F14.7)') ' EROR  =',EROR
        WRITE (MT6,'(a,F14.7)') ' DR    =',DR
      ENDIF
 
      ISOLN = 1; IALT = 0; JALT = 0;  KALT = 0

! ... ISOLN = 1 means that the solutions and their derivatives
! ... will be calculated at RO, not only solutions in R1 and R2,
! ... and then the parameters IALT,JALT,KALT are ignored.


! ... no output for radial function on some mesh: 

      IFILE = 0;  JFILE = 0;  KFILE = 0

! ... ion chaege and asymptotic multipoles:

      SZ = ION;   LAMMX = MAX(km,1);  if(LAMMX.gt.MLMX) LAMMX=MLMX
 
      DO K = 1,LAMMX
        kphase= (-1)**k                                     !  ???
        BLAM(1:NCH,1:NCH,K) = CF(1:NCH,1:NCH,K)  !* kphase  !  ???
      END DO
 
! ... orbital momentums and channel energies:

      NCHF = NCH;  NCHOP = NOPEN

      DO I = 1,NCH
        BLAM(I,I,1) = BLAM(I,I,1) - (LCH(I)* (LCH(I)+1))
        EL(I) = LCH(I)
        EN(I) = ECH(I)
        if(EN(i).eq.0.d0) then; info=3; Return; end if 
      END DO
 
! ... the following parameters are important only for the case
! ... when all channels are closed: 

      DEL1 = 0.0D0;  ECON = 0.0D0;  METH = 1; 
 
      IF(IBUG.GT.1) THEN
        WRITE (MT6,'(a,i2)') ' NCHAN =',NCHF
        WRITE (MT6,'(a,i2)') ' NOPEN =',NOPEN
        WRITE (MT6,'(a/(6F12.5))') ' CHANNEL ENERGIES ...', &
                                    (EN(I),I=1,NCH)
      END IF

! ... CALL ASYMPTOTIC PACKAGE:
 
        CALL ASYPCK(IAUTO,FA,FAD)
 
! ... RENORMALISE OPEN CHANNEL SOLUTIONS:
 
      if(NCHOP.gt.0) then

       if(SZ.le.0.d0) SZ = 1.0;    A1 = SQRT(SZ)
       FA (1:NCHF,1:NCHOP,1:2) = FA (1:NCHF,1:NCHOP,1:2) / A1
       FAD(1:NCHF,1:NCHOP,1:2) = FAD(1:NCHF,1:NCHOP,1:2) / A1
      
      end if


! ... find output:

      F (1:nch,1:nchop) = FA (1:NCHF,1:NCHOP,1)  ! sine solution      
      FP(1:nch,1:nchop) = FAD(1:NCHF,1:NCHOP,1)  ! its derivative      

      G (1:nch,1:nchop) = FA (1:NCHF,1:NCHOP,2)  ! cosine solution      
      GP(1:nch,1:nchop) = FAD(1:NCHF,1:NCHOP,2)  ! its derivative      

      
      if(nchop.lt.nch) then
       m = nchop + 1
       G (1:nch,m:nch) = FA (1:NCHF,m:NCHF,1) 
       GP(1:nch,m:nch) = FAD(1:NCHF,m:NCHF,1)       
       F (1:nch,m:nch) = 0.d0 
       FP(1:nch,m:nch) = 0.d0       
      end if

      END SUBROUTINE ZAFACE
 


 
