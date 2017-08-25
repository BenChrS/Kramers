#!/bin/bash

tau=2
a=3.1
chi=6
potw=1
  #  D=666.6667 #25.5692 #2    #35.77709
  # temperature=1.0
B_eff=0
npTen=0
npTwo=3					
avOpt=1
avNum=1
transition=0
npoints=10
pdaBool=0
vdaBool=0
corrFuncNr=3
potNr=0   # 0:kein Potenzial 1: harmonischer Oszillator 4:Mexican-hat 5: abschnitsweise definiertes Potenzial 6: "" mit modifizierter
noiseNr=2 # Frequenz
sdeSolverNr=3
Kramers=0
initCondNr=0 #0: x,v fest, 1: v thermisch verteilt, 2: x th. vert., 3: x,v th. ver 
testBool=0
a_Corr=9
b_Corr=10
pota=1
potMy=1.99999866
potNy=0

cat > input.txt <<EOF
$tau
$a
$chi
$potw
$B_eff
$npTen
$npTwo
$avOpt
$avNum
$transition
$npoints
$pdaBool
$vdaBool
$corrFuncNr
$potNr
$noiseNr
$sdeSolverNr
$Kramers
$initCondNr
$testBool
$a_Corr
$b_Corr
$pota
$potMy
$potNy
EOF


./colnoise

