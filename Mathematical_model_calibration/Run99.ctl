$PROBLEM Indirect response model

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
  T12 = 1 ; for those patients without treatment; not used as no dose given
  ; https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.9b01191

  IF (DRUG.EQ.3) THEN ; Eylea
    T12 = 9
  ENDIF
  IF (DRUG.EQ.4) THEN ; Lucentis
    T12 = 6.5
  ENDIF
  KEL = LOG(2)/T12

MU_1  = THETA(1);                   
CSTDEL = EXP(MU_1 + ETA(1)) ; CSTDEL

MU_2  = THETA(2);                   
KIN   = EXP(MU_2 + ETA(2)) ; KIN

MU_3  = THETA(3);                   
EC50_Eylea = EXP(MU_3 + ETA(3)) ; EC50_Eylea ; Sharing ETA(3) with EC50_Lucentis

MU_4  = THETA(4);                   
EC50_Lucentis = EXP(MU_4 + ETA(3)) ; EC50_Lucentis ; Sharing ETA(3) with EC50_Eylea

MU_5  = THETA(5);                   
CSTMIN = EXP(MU_5 + ETA(4)) ; CSTMIN

MU_6  = THETA(6);                   
CST0  = EXP(MU_6 + ETA(5)) ; CST0

MU_7  = THETA(7);                   
PK0   = EXP(MU_7 + ETA(6)) ; PK0

MU_8  = THETA(8);                   
HILL  = EXP(MU_8 ) ; HILL

EC50=1

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

$DES
DADT(1) = - KEL*A(1)
DADT(2) =   KIN - KOUT*(1 + EMAX*((A(1) + 1E-08)**HILL)/(EC50**HILL+(A(1) + 1E-08)**HILL))*A(2)

$ERROR
PKC    = A(1)
CST    = A(2)

ADD_ERR = THETA(9) ; Residual error for CST

IPRED = CST
W     = SQRT((ADD_ERR)**2)
Y     = IPRED + W*ERR(1) 

IRES = IPRED - DV
IWRES = IRES / W

$THETA
 3        ;  1 CSTDEL
 4          ;  2 KIN
-1.5         ;  3 EC50_Eylea
-5          ;  4 EC50_Lucentis
(5, 5.5)    ;  5 CSTMIN
 5.5        ;  6 CST0
-2          ;  7 PK0
(0,1)       ;  8 HILL
;=========================================================================   
(0, 5)      ; 9 additive error CST
;=========================================================================  
$OMEGA 
1           ; ETA(1) for CSTDEL
0 FIX       ; ETA(2) for KIN
3          ; ETA(3) for EC50 parameters
0.2           ; ETA(4) for CSTMIN
0.2         ; ETA(5) for CST0
4           ; ETA(6) for PK0

;=========================================================================  
$SIGMA
1 FIX 
;========================================================================= 

$EST  MAXEVAL=9999 NSIG=6 SIGL=9 PRINT=5 METHOD=1 INTER NOABORT
NOTHETABOUNDTEST NOOMEGABOUNDTEST NOSIGMABOUNDTEST

$COV UNCONDITIONAL
$TABLE  C IDORIG ID TIME DV MDV CMT AMT TAD AGE CENTER EYE DRUG DIA DOSE PKC CST MU_1 MU_2 MU_3 MU_4 MU_5 MU_6 MU_7 MU_8 ADD_ERR ETAS(1:6) PRED IPRED IRES IWRES NPDE CWRES
NOPRINT ONEHEADER NOAPPEND
FILE = run1.fit
