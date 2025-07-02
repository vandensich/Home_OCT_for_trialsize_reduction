$PROBLEM Indirect response model

$INPUT TIME DV ID CMT AMT DOSE MDV AGE CENTER DRUG DIA C EYE TAD
$DATA Data_154_Notal_longer_handpicked.csv

IGNORE=@
IGNORE=(CMT.GT.2)
IGNORE=(DV.EQ.-1)

IGNORE=(TIME.GT.168)

$SUBROUTINE ADVAN13 TOL=9

$MODEL
      NCOMP=2
      COMP=(PK)   ;1
      COMP=(CST)  ;2

$PK
  T12 = 1 ; for those patients without treatment; not used as no dose given
  ; https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.9b01191

  
  T12 = 9
 
  KEL = LOG(2)/T12

MU_1  = THETA(1);                   
CSTDEL = EXP(MU_1 + ETA(1)) ; CSTDEL

MU_2  = THETA(2);                   
KIN   = EXP(MU_2 + ETA(2)) ; KIN

MU_3  = THETA(3);                   
EC50 = EXP(MU_3 + ETA(3)) ; EC50_Trialdrug ; Sharing ETA(3) with EC50_Lucentis



MU_4  = THETA(4);                   
CSTMIN = EXP(MU_4 + ETA(4)) ; CSTMIN

MU_5  = THETA(5);                   
CST0  = EXP(MU_5 + ETA(5)) ; CST0

MU_6  = THETA(6);                   
PK0   = EXP(MU_6 + ETA(6)) ; PK0
            
MU_7  = THETA(7);                   
HILL  = EXP(MU_7 ) ; HILL


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

ADD_ERR = THETA(8) ; Residual error for CST

IPRED = CST
W     = SQRT((ADD_ERR)**2)
Y     = IPRED + W*ERR(1) 

IRES = IPRED - DV
IWRES = IRES / W

$THETA
 3.5        ;  1 CSTDEL
 4          ;  2 KIN
-2          ;  3 EC50
(5, 5.4)    ;  4 CSTMIN
 5.5        ;  5 CST0
-3.61 FIX  ;  6 PK0
(0,1)       ;  7 HILL
;=========================================================================   
(0, 5)      ; 8 additive error CST
;=========================================================================  
$OMEGA 
1           ; ETA(1) for CSTDEL
0 FIX       ; ETA(2) for KIN
3           ; ETA(3) for EC50
0.2         ; ETA(4) for CSTMIN
0.1         ; ETA(5) for CST0
0 FIX       ; ETA(6) for PK0

$SIGMA
1 FIX 
;========================================================================= 
$SIML (1234567) BOOTSTRAP=20 NOREPLACE SUBP=200

$EST  MAXEVAL=9999 NSIG=5 SIGL=9 PRINT=5 METHOD=1 INTER NOABORT
NOTHETABOUNDTEST NOOMEGABOUNDTEST NOSIGMABOUNDTEST
$COV UNCONDITIONAL
$TABLE  C ID TIME DV MDV CMT AMT TAD AGE CENTER EYE DRUG DIA DOSE PKC CST MU_1 MU_2 MU_3 MU_4 MU_5 MU_6 MU_7 ADD_ERR ETAS(1:6) PRED IPRED IRES IWRES NPDE CWRES
NOPRINT ONEHEADER NOAPPEND
FILE = run1.fit
