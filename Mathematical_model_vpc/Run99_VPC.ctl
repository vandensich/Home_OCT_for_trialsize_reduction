$PROBLEM Simulation for PC-VPC

$INPUT C IDORIG ID TIME DV MDV CMT AMT TAD AGE CENTER EYE DRUG DIA DOSE DROP DROP DROP DROP DROP DROP
$DATA NOTAL_NM_20230421_JH.csv

IGNORE=@
IGNORE=(CMT.GT.2)
IGNORE=(DV.EQ.-1)

IGNORE=(ID.EQ.10021)
IGNORE=(ID.EQ.10022)
IGNORE=(ID.EQ.10031)
IGNORE=(ID.EQ.10042)
IGNORE=(ID.EQ.10052)
IGNORE=(ID.EQ.10051)
IGNORE=(ID.EQ.10102)
IGNORE=(ID.EQ.20011)
IGNORE=(ID.EQ.20031)

$SUBROUTINE ADVAN13 TOL=8

$MODEL
      NCOMP=2
      COMP=(PK)   ;1
      COMP=(CST)  ;2

$PK
  IF (DRUG.EQ.3) THEN ; Eylea
    T12 = 9
  ENDIF
  IF (DRUG.EQ.4) THEN ; Lucentis
    T12 = 6.5
  ENDIF
  KEL = LOG(2)/T12
  
MU_1  = THETA(1);                   
CSTDEL = EXP(MU_1 + ETA(1)) 

MU_2  = THETA(2);                   
KIN   = EXP(MU_2 + ETA(2)) 

MU_3  = THETA(3);                   
EC50_Eylea = EXP(MU_3 + ETA(3)) 

MU_4  = THETA(4);                   
EC50_Lucentis = EXP(MU_4 + ETA(3)) 

MU_5  = THETA(5);                   
CSTMIN = EXP(MU_5 + ETA(4)) 

MU_6  = THETA(6);                   
CST0  = EXP(MU_6 + ETA(5)) 

MU_7  = THETA(7);                   
PK0   = EXP(MU_7 + ETA(6)) 

MU_8  = THETA(8);                   
HILL  = EXP(MU_8 ) 

IF (DRUG.EQ.3) THEN ; Eylea
  EC50 = EC50_Eylea
ENDIF
IF (DRUG.EQ.4) THEN ; Lucentis
  EC50 = EC50_Lucentis
ENDIF

S1 = 1
S2 = 1

EMAX = CSTDEL / CSTMIN
KOUT = KIN / (CSTMIN + CSTDEL)

A_0(1) = PK0
A_0(2) = CST0



$THETA  3.55   ; 1 CSTDEL  (update with final estimate)
$THETA  4.01      ; 2 KIN     (update with final estimate)
$THETA -1.67      ; 3 EC50_Eylea
$THETA -3.82      ; 4 EC50_Lucentis
$THETA 5.42      ; 5 CSTMIN
$THETA 5.54      ; 6 CST0
$THETA -1.2      ; 7 PK0
$THETA  0.779      ; 8 HILL
$THETA  5.49      ; 9 additive error CST

 
$OMEGA  1.02      ; ETA(1) for CSTDEL
$OMEGA  0.00 FIX     ; ETA(2) for KIN (fixed)
$OMEGA  2.82      ; ETA(3) for EC50 parameters
$OMEGA  3.34E-2      ; ETA(4) for CSTMIN
$OMEGA  3.83E-2      ; ETA(5) for CST0
$OMEGA  4.54      ; ETA(6) for PK0

$SIGMA
  1.0 FIX  

$DES
DADT(1) = - KEL*A(1)
DADT(2) =   KIN - KOUT*(1 + EMAX*((A(1) + 1E-08)**HILL)/(EC50**HILL+(A(1) + 1E-08)**HILL))*A(2)

$ERROR
$ERROR
PKC    = A(1)
CST    = A(2)

ADD_ERR = THETA(9) ; Residual error for CST

IPRED = CST
W     = SQRT((ADD_ERR)**2)
Y     = IPRED + W*ERR(1) 

IRES = IPRED - DV
IWRES = IRES / W



$EST  MAXEVAL=9999 NSIG=6 SIGL=9 PRINT=5 METHOD=1 INTER NOABORT

