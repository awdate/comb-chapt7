
          program react
          DIMENSION TD(200),HPR(200)
          OPEN(UNIT=8,FILE='OO')
          OPEN(UNIT=9,FILE='O1')
C THIS IS VARIABLE VOL COMBUSTION - PROB 12 - CHAPTER 7 - Expansion stroke
c See new source term S4 
C  CmHn + RST/FI O2 + n{H2O} H2O ( 2 % ) ----> 2 CO2 + (n/2+n_H2O)H_{2}O 
C   T0 = 300 K ,  P0 = 1 ATM
C COMPRESSION RATIO ( CR ) = 10.5 ( PV^{1.4} = C ) RPM=1200, AMMIX=29
C SIMULATION OF EXPLOSION AT tdc
C STUDY EFFECT OF FI 
          WRITE(*,*)' READ TWALL,RPM  '
          READ(*,*)TWALL,RPM
          FI=1
          PI=22.0/7.0
          CR=10.5
          P0=101325
          T0=300
          TIN=T0
          PIN=P0
C          
          XH2O=0.02
          GAMMA=1.4
          BORE=0.075
          STROKE=0.075
          CRANK=0.0375
          CROD=0.30
          AC=PI*BORE**2/4.0
          OMEGA=2*PI*RPM/60.0
C CLEARANCE
          X0=STROKE/(CR-1.0)
          VCV0=PI*BORE**2/4.0*X0
C FUEL - ISO OCTANE C8H18 OR C2H6
C          AM=2
C          AN=6
          AM=8.0
          AN=18.0
          AF=4.049E10
          EBYRU=20131
          XX=0.25
          YY=1.5
          ASTOIC=AM+AN/4.0
          RU=8314
          WFU=AM*12+AN
          WO2=32
          WH2O=18
          WCO2=44
          WN2=28
         write(8,*)' TIN = ',TIN,' pin = ',pin
         write(*,*)' TIN = ',TIN,' pin = ',pin
C DETERMINE ADIABATIC TEMPERATURE
          hfutin=hfu(tin)
          ho2tin=ho2(tin)
          hn2tin=hn2(tin)
          hh2otin=hh2o(tin)
         ANH2O=XH2O/(1.0-XH2O)*(1.0+ASTOIC/FI*4.76)
         HREACT=HFUTIN+ASTOIC/FI*(HO2TIN+3.76*HN2TIN)
     1         +ANH2O*HH2OTIN
          
          write(8,*)' This is Ex-Prob12, chapter 7 '
          write(8,*)' TIN = ',TIN,' Hfu = ', HFUTIN
          write(8,*)' HO2 = ',HO2TIN
          write(8,*)' HN2 = ',HN2TIN,' HH2O = ',HH2OTIN
          write(8,*)' ANH2O = ',ANH2O
          write(8,*) ' Astoic = ',astoic,' HREACT = ',HREACT
          TD(1)=TIN+500
          i=0
21        i=i+1
          TAD=TD(I)
          HPR(I)=AM*HCO2(TAD)
     1    +( AN/2.0 + ANH2O)*HH2O(TAD)
     1   + (ASTOIC/FI-ASTOIC)*HO2(TAD)
     1   + ASTOIC/FI*3.76*HN2(TAD)
c          WRITE(8,*)I,TD(I),HPR(I)
          TD(I+1)=TD(I)+100
          ilast=i+1
          IF(TD(I+1).LT.2800.0)GO TO 21
          write(8,*)' ilast = ',ilast
c
          DO 22 I=1,ilast-1
          IF(HPR(I).LT.HREACT.AND.HPR(I+1).GT.HREACT)THEN
c          IF(HPR(I).LT.0.0.AND.HPR(I+1).GT.0.0)THEN
          RATIO=(HREACT-HPR(I))/(HPR(I+1)-HPR(I))
          TREQ=(TD(I+1)-TD(I))*RATIO+TD(I)
          ENDIF
22        CONTINUE          
          TAD=TREQ
c          WRITE(8,*)TAD,HREACT
C DETERMINE H0F
         HPROD=AM*HCO2(TIN)
     1   +( AN/2.0 + ANH2O)*HH2O(TIN)
     1   + (ASTOIC/FI-ASTOIC)*HO2(TIN)
     1   + ASTOIC/FI*3.76*HN2(TIN)
C
          H0F=(HREACT-HPROD)*1000.0/WFU
c       write(8,*)HCO2(TIN),HO2(TIN),HN2(TIN),HH2O(TIN),HPROD,H0F
C DETERMINE CP AND CV
          CPMIX=1200
          AMMIX=29
          RMIX=8314/AMMIX
          CVMIX=CPMIX-RMIX
c extinction temperature
          TERM=WFU+ASTOIC/FI*(WO2+3.76*WN2)+ANH2O*WH2O
          OFUI=WFU/TERM
          TEXT=TI+OFUI*H0F/CVMIX
          RHOI=PIN*AMMIX/(RU*TIN)
**************************************************************
C TSTART AND PSTART AT T = 1000 FROM THE  PROBLEM 10
***************************************************************
          TOLD=1000
          POLD=35.6E5
          OFUI=0.0567
          OO2I=0.199
          OH2OI=0.0188
          OCO2I=0.0149
          ON2I=1-ofui-oo2i-oco2i-oh2oi
          WRITE(8,*) ' ************************************* '
          WRITE(8,*)' FI = ',FI,' PSTART = ', POLD,' TSTART = ',TOLD
     1           ,' Tad = ',TAD,' RPM = ',RPM, ' OMEGA = ',OMEGA 
          WRITE(8,*)' OO2I = ',OO2I,' OFUI = ',OFUI,' ON2I = ',ON2I
          WRITE(8,*)' OCO2I = ',OCO2I,' OH2Oi = ',OH2OI
          WRITE(8,*)' TEXTINCTION = ',TEXT, ' TWALL = ',TWALL
          WRITE(8,*)' AMMIX = ',AMMIX,' CVMIX = ',CVMIX
          WRITE(8,*)' RMIX = ',RMIX,' CPMIX = ',CPMIX          
          write(8,*)' ANH2O = ',ANH2O,' H0f = ',H0f,' rhoi = ',rhoi
          WRITE(8,*) ' ************************************* '
          write(8,*)' ***** solution begins ' 
C SOLVE FOR T -CYCLE BEGINS - DTREF ( s )
          DTREF=0.00001
          TIME=0.0
          DELT=DTREF
          WRITE(8,*)' DELT = ',DELT,' SECONDS'
          WRITE(8,*)' TIME(ms)    Theta      Tnew       OFU        OO2
     1        P        RHON         QW     alfa'

          WRITE(9,*)' TIME      THETA       alfa1       alfa2       pad              
     1     ug      rfu      sor'
C
          OFUO=OFUI
          OO2O=OO2I
          OCO2O=OCO2I
          OH2OO=OH2OI
          RHOO=POLD*AMMIX/(RU*TOLD)
          VCVOLD=VCV0
C BEGIN A STEP
          ISTEP=0
          THETA=0.0
1000      ISTEP=ISTEP+1
c time is in seconds
          TIME=TIME+DELT
          THETA=THETA+OMEGA*DELT
          TERM=SQRT(CROD**2-(CRANK*SIN(THETA))**2)
          XP=CRANK*(1.-COS(THETA))+CROD-TERM
          UP=CRANK*OMEGA*SIN(THETA)*(1+CRANK*COS(THETA)/TERM)
C DETREMINE RHO - MASS CONS
          RHON=RHOO/(1.+UP*DELT/(X0+XP))
C ENERGY EQN
          AKF=AF*WFU
          RFU=AKF*EXP(-EBYRU/TOLD)*(OFUO)**XX*OO2O**YY
          RFU=RFU*RHON*(XX+YY)/(WFU**XX*WO2**YY)
          AW=PI*BORE*(XP+X0)
          VCV=PI*BORE**2/4.*(XP+X0)
C HEAT TRANSFER COEFFICIENT - BORMAN 
          PAD=PIN*(VCV0/VCV)**GAMMA
          P=RHON*RMIX*TOLD
          IF(P.LT.PIN)P=PIN
          UPAVE=2*STROKE*RPM/60.0
          UG=2.28*UPAVE+0.00324*VCV*(P-PAD)*TIN/(VCV0*PIN)
c kg=0.08 w/m-K -- 
          akg=0.08
          amug=608e-7
          anug=amug/rhon
c worchni correlation
          alfa1=0.035*akg/bore*(ug*bore/anug)**0.8
c ragland and borman correlation
          ALFA2=820*(UG*P/1e6)**0.8/TOLD**0.53/BORE**0.2
          alfa=(alfa1+alfa2)/2.0
          QW=ALFA*(TWALL-TOLD)
C HEAT GEN
          S2=ABS(H0F*RFU)
C HEAT LOSS
          S3=QW*AW/VCV
C WORKDONE
         s4=pold*up/(x0+xp)
C ENERGY EQN
          SOURCE=S2+S3-S4
          TNEW=TOLD+DELT*SOURCE/(RHON*CVMIX)
          IF(TNEW.LT.TIN)TNEW=TIN
          IF(TNEW.GT.TAD)TNEW=TAD
          TSTAR=(TNEW-TIN)/(TAD-TIN)
C  OFU BECAUSE QW IS FINITE
          OFU=OFUO-DELT*RFU/RHON
C
          WRITE(*,*)' ISTEP = ',ISTEP,' TIME = ',TIME,' SEC'
          WRITE(*,*)' S2 = ',S2,' S3 = ',S3,' s4 = ',s4
          write(*,*)' source = ',source
          WRITE(*,*)' alfa1 = ',alfa1,' alfa2 = ', alfa2
          OO2=ASTOIC*WO2/WFU*(OFU-OFUI)+OO2I
          OCO2=AM*WCO2/WFU*(OFUI-OFU)+OCO2I
          OH2O=(AN/2.0)*WH2O/WFU*(OFUI-OFU)+OH2OI
          DPDT=(P-POLD)/DELT/1E5
          PP=P/1E5
          OCO2P=OCO2*10
          TMS=TIME*1000
          ANG=theta*180.0/PI
          WRITE(9,500)TMS,ANG,alfa1,alfa2,pad,ug,rfu,source
          WRITE(8,500)TMS,ANG,Tnew,OFU,OO2,PP,RHON,QW,alfa
500     FORMAT(E10.5,9(1x,E10.3))         
          WRITE(*,*)' Ofu = ',OFU,' OO2 = ',OO2,' OCO2 = ',OCO2
     1              ,' OH2O = ',OH2O,' PNEW = ',P          
c cylinder mass - must be const
          cmass=rhon*vcv
          WRITE(*,*)' RFU = ' ,RFU,' RHON = ',RHON,' cmass = ',cmass
C SET OLD VALUES
          IF(OO2.LT.0.0)go to 2000
          IF(Ofu.LT.0.0)go to 2000
          if(source.eq.0.0)go to 2000
          TOLD=TNEW
          OFUO=OFU
          OO2O=OO2          
          OCO2O=OCO2
          OH2OO=OH2O
          POLD=P
          RHOO=RHON
          VCVOLD=VCV
          IF(THETA.LT.PI)GO TO 1000
2000      write(8,*)'******************* End of Burning *************'
          write(9,*)'******************* End of Burning *************'
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

