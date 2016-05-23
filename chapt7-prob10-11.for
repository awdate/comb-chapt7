          program react
          DIMENSION TD(20),HPR(20)
        OPEN(UNIT=8,FILE='OO')
C THIS IS CONST VOL COMBUSTION - PROBs 10 and 11 - CHAPTER 7
c See FUNCTION ROUTINES FOR C2H6 AND C8h18
C FUEL CmHn + RST/FI O2 + n_{H2O} H2O ( 2 % ) ----> 2 CO2 + (3+n)H2O 
C   T0 = 300 K ,  P0 = 1 ATM
C COMPRESSION RATIO ( CR ) = 10 ( PV^{1.4} = C ) - CPM=1200 MMIX=29
C SIMULATION OF EXPLOSION AT TDC
C PROGRAM IS WRITTEN FOR 2 FUELS - C2H6 AND C8H18
C STUDT EFFECT OF FI 
          WRITE(*,*)' READ CR  '
          READ(*,*)CR
          FI=1.0
          TI=300*CR**0.4
          PI=101325*CR**1.4
          XH2O=0.02
C For Ethane C2H6
c          AM=2
c          AN=6
c          xx=0.1
c          yy=1.65
c          afu=6.1857e9
c          efu=15098
c For iso-octane

          AM=8
          AN=18
          xx=0.25
          yy=1.5
          afu=4.049e10
          efu=20131

          ASTOIC=AM+AN/4.0
          RU=8314
          WFU=AM*12+AN
          WO2=32
          WH2O=18
          WCO2=44
          WN2=28
C DETERMINE ADIABATIC TEMPERATURE
         ANH2O=XH2O/(1.0-XH2O)*(1.0+ASTOIC/FI*4.76)
         HREACT=HFU(TI)+ASTOIC/FI*(HO2(TI)+3.76*HN2(TI))
     1         +ANH2O*HH2O(TI)
          write(8,*)' This is Ex-Prob10, chapter 7 '
          write(8,*)' TI = ',Ti,' Hfu = ', HFU(TI),' HO2 = ',HO2(TI)
     1   ,  ' HN2 = ',HN2(TI),' HH2O = ',HH2O(TI),' ANH2O = ',ANH2O
     1   , ' Astoic = ',astoic,' HREACT = ',HREACT
          TD(1)=TI+500
          i=0
21        i=i+1
          TAD=TD(I)
          HPR(I)=AM*HCO2(TAD)
     1    +( AN/2.0 + ANH2O)*HH2O(TAD)
     1   + (ASTOIC/FI-ASTOIC)*HO2(TAD)
     1   + ASTOIC/FI*3.76*HN2(TAD)
          WRITE(8,*)I,TD(I),HPR(I)
          TD(I+1)=TD(I)+100
          ilast=i+1
          IF(TD(I+1).LT.2800.0)GO TO 21
          write(8,*)' ilast = ',ilast
c
          DO 22 I=1,ilast-1
c          IF(HPR(I).LT.HREACT.AND.HPR(I+1).GT.HREACT)THEN
          IF(HPR(I).LT.0.0.AND.HPR(I+1).GT.0.0)THEN
          RATIO=(HREACT-HPR(I))/(HPR(I+1)-HPR(I))
          TREQ=(TD(I+1)-TD(I))*RATIO+TD(I)
          ENDIF
22        CONTINUE          
          TAD=TREQ
          WRITE(8,*)TAD,HREACT
C DETERMINE H0F
         HPROD=AM*HCO2(TI)
     1   +( AN/2.0 + ANH2O)*HH2O(TI)
     1   + (ASTOIC/FI-ASTOIC)*HO2(TI)
     1   + ASTOIC/FI*3.76*HN2(TI)
C
          H0F=(HREACT-HPROD)*1000.0/WFU
          write(8,*)HCO2(TI),HO2(TI),HN2(TI),HH2O(TI),HPROD,H0F
C DETERMINE CP 
          CPMIX=1200
c extinction temperature
          TERM=WFU+ASTOIC/FI*(WO2+3.76*WN2)+ANH2O*WH2O
          OFUI=WFU/TERM
          TEXT=TI+OFUI*H0F/CVMIX
          OO2I=ASTOIC/FI*WO2/TERM
          ON2I=ASTOIC/FI*3.76*WN2/TERM
          OH2OI=ANH2O*WH2O/TERM
          OCO2I=0.0
          WRITE(8,*) ' ************************************* '
          WRITE(8,*)' FI = ',FI,' PI = ', PI,' TI = ',TI,' Tad = ',TAD
          WRITE(8,*)' OO2I = ',OO2I,' OFUI = ',OFUI,' ON2I = ',ON2I
          WRITE(8,*)' OH2OI = ',OH2OI,' TEXT = ',TEXT
          write(8,*)' ANH2O = ',ANH2O,' H0f = ',H0f,' CPMIX = ',CPMIX
C SOLVE FOR T ( DTREF ms )
          DTREF=0.1
          TIME=0.0
          DELT=DTREF
          WRITE(8,*)' DELT = ',DELT,' milli SECONDS'
          WRITE(8,*)' TIME ( ms )    TSTAR    TT     OFU     OO2    OCO2P 
     1    OH2O     P ( bars )    DPDT    CHECK'

          WRITE(8,*) ' ************************************* '

          TOLD=TI
          OFUO=OFUI
          OO2O=OO2I
          OCO2O=OCO2I
          OH2OO=OH2OI
          POLD=PI
C BEGIN A STEP
          ISTEP=0
1000      ISTEP=ISTEP+1
c TIME IS IN ms
          TIME=TIME+DELT
C CALCULATE RFU          
          term=ofuo/wfu+oo2o/wo2+oco2o/wco2+oh2oo/wh2o+on2i/wn2
          AMMIX=1./term
          RMIX=RU/AMMIX
          CVMIX=CPMIX-RMIX
          RHO=POLD/(RMIX*TOLD)
          AKF=AFU*WFU
          RFU=AKF*EXP(-EFU/TOLD)*OFUO**XX*OO2O**YY
          RFU=RFU*RHO**(XX+YY)/(WFU**XX*WO2**YY)
          SOURCE=ABS(H0F/CVMIX/RHO*RFU)
c delt/1000 = dt ( s )
          DT=DELT/1000.0
          TNEW=TOLD+DT*SOURCE
          TSTAR=(TNEW-TI)/(TAD-TI)
C
          WRITE(*,*)' ISTEP = ',ISTEP,' TIME = ',TIME,' mili SEC'
          WRITE(*,*)' TNEW = ',TNEW
          OFU=OFUI-CVMIX/H0F*(TNEW-TI)
          OO2=ASTOIC*WO2/WFU*(OFU-OFUI)+OO2I
          OCO2=AM*WCO2/WFU*(OFUI-OFU)+OCO2I
          OH2O=(AN/2.0)*WH2O/WFU*(OFUI-OFU)+OH2OI
          P=PI/TI*TNEW
          DPDT=(P-POLD)/DELT/1E5
          TT=TNEW/1000
          PP=P/1E5
          OCO2P=OCO2*10
          CHECK=ON2I+OFU+OO2+OH2O+OCO2
          WRITE(8,500)TIME,TSTAR,TT,OFU,OO2,OCO2P,OH2O,PP,DPDT,CHECK
500     FORMAT(E10.5,9(1x,E10.3))         
          WRITE(*,*)' Ofu = ',OFU,' OO2 = ',OO2,' OCO2 = ',OCO2
     1      ,' OH2O = ',OH2O,' PNEW = ',P,' Mmix = ',AMMIX
          WRITE(*,*)' RFU = ' ,RFU,' RHO = ',RHO
C SET OLD VALUES
          TOLD=TNEW
          OFUO=OFU
          OO2O=OO2          
          OCO2O=OCO2
          OH2OO=OH2O
          POLD=P
C fOR C2H6
c          IF(TSTAR.GT.0.1)DELT=0.1*DTREF
c          IF(TSTAR.GT.0.15)DELT=0.05*DTREF
c          IF(TSTAR.GT.0.20)DELT=0.01*DTREF
c          IF(TSTAR.GT.0.25)DELT=0.005*DTREF
c          IF(TSTAR.GT.0.35)DELT=0.003*DTREF
C FOR C8H18
          IF(TSTAR.GT.0.1)DELT=0.5*DTREF
          IF(TSTAR.GT.0.15)DELT=0.2*DTREF
          IF(TSTAR.GT.0.20)DELT=0.1*DTREF
          IF(TSTAR.GT.0.25)DELT=0.05*DTREF
          IF(TSTAR.GT.0.35)DELT=0.01*DTREF
c          IF(TSTAR.LT.1.0)GO TO 1000
          IF(OFU.GT.0.0)GO TO 1000
          STOP
          END                     
C**************************************************
          FUNCTION HFU(TT)
C**************************************************
          T=TT
C COEFFICIENTS OF C8H18
         a1=12.524584
         a2=-1.0101883e-2
         a3=2.219926e-4
         a4=-2.8486372e-7
         a5=1.1241014e-10
         b1=-2.98434406e4
         b2=-19.710999
         IF(T.GT.1000)THEN
         a1=20.943071
         a2=4.4169102e-2
         a3=-1.5326163e-5
         a4=2.305448e-9
         a5=-1.2976573e-13
         b1=-3.5575509e4
         b2=-81.063773
         ENDIF
         TERM=A1*T+A2/2.*T**2+A3/3.*T**3+A4/4.*T**4+A5/5.*T**5
         TERM=TERM+B1
         HFU=TERM*8.314
         RETURN
         END

C**************************************************
          FUNCTION HFU1(TT)
C**************************************************
          T=TT
C COEFFICIENTS OF C2H6
         a1=4.2914257
         a2=5.5015495e-3        
         a3=5.9943846e-5
         a4=-7.0846647e-8
         a5=2.686858e-11
         b1=-1.1522206e4
         b2=2.6667899
         IF(T.GT.1000)THEN
         a1=4.0466641
         a2=1.535388e-2
         a3=-5.4703949e-6
         a4=8.7782654e-10
         a5=-5.2316753e-14
         b1=-1.244735e4
         b2=-0.96869831
         ENDIF
         TERM=A1*T+A2/2.*T**2+A3/3.*T**3+A4/4.*T**4+A5/5.*T**5
         TERM=TERM+B1
         HFU=TERM*8.314
         RETURN
         END
C**************************************************
          FUNCTION HCO2(TT)
C**************************************************
          T=TT
C COEFFICIENTS OF CO2
         a=24.99735
         b=55.18696
         c=-33.69137
         d=7.948387
         e=-0.136638
         f=-403.6075
         g=228.2431
         IF(T.GT.1200)THEN
         a=58.16639
         b=2.720074
         c=-0.492289
         d=0.038844
         e=-6.447293
         f=-425.6075
         g=263.6125
         ENDIF
         H0F=-393520
         TH=T/1000.0
         TERM=A*TH+B/2.*TH**2+C/3.*TH**3+D/4.*TH**4
     1    -E/TH+F-H0F/1000.0
         DHSCO2=1000.0*TERM
         HCO2=H0F+DHSCO2
         RETURN
         END
C**************************************************
          FUNCTION HO2(TT)
C**************************************************
          T=TT
C COEFFICIENTS OF O2
         a=29.659
         b=6.137261
         c=-1.186521
         d=0.09578
         e=-0.219663
         f=-9.861391
         g=237.948
         H0F=0.0
         TH=T/1000.0
         TERM=A*TH+B/2.*TH**2+C/3.*TH**3+D/4.*TH**4
     1    -E/TH+F
         DHSO2=1000*(TERM-HOF/1000.0)
         HO2=H0F+DHSO2
         RETURN
         END
C**************************************************
          FUNCTION HN2(TT)
C**************************************************
         T=TT
C COEFFICIENTS OF N2
         a=26.092
         b=8.218801
         c=-1.976141
         d=0.159274
         e=0.044434
         f=-7.98923
         g=221.02
         H0F=0.0
         TH=T/1000.0
         TERM=A*TH+B/2.*TH**2+C/3.*TH**3+D/4.*TH**4
     1    -E/TH+F
         DHSN2=1000*(TERM-HOF/1000.0)
         HN2=H0F+DHSN2
         RETURN
         END

C**************************************************
          FUNCTION HH2O(TT)
C**************************************************
         T=TT
C COEFFICIENTS OF H2O
         a=30.092
         b=6.832514
         c=6.793435
         d=-2.53448
         e=0.0823139
         f=-250.881
         g=223.3967
         IF(T.GT.1700)THEN
         a=41.96426
         b=8.622053
         c=-1.49978
         d=0.098119
         e=-11.15764
         f=-272.1787
         g=219.7809
         ENDIF
         H0F=-241830
         TH=T/1000.0
         TERM=A*TH+B/2.*TH**2+C/3.*TH**3+D/4.*TH**4
     1    -E/TH+F-H0F/1000.0
         DHSH2O=1000*TERM
         HH2O=H0F+DHSH2O
         RETURN
         END

