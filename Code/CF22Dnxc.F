Subroutine CF22Dx(F,D1F,RA,RB,D1RA,D1RB,TA,TB,NGrid,ijzy)
**************************************************************
*                                                            *
*  CF22Dnxc evaluates the nonseparable exchange-correlation  *
*      part of the CF22D functional on the density grid.     *
*                                                            *
*  OUTPUT:                                                   *
*     F   - Functional values                                *
*     D1F - First derivatives with respect to RA, RB, GA, GB *
*              TA, TB                                        *
*                                                            *
*  INPUT:                                                    *
*       ijzy - 1 CF22D                                       *
*                                                            *
*     RA,B   - Spin densities                                *
*     D1RA,B - Spin density gradients                        *
*     TA,B   - Spin kinetic energy densities                 *
*     NGrid  - number of grid points                         *
*                                                            *
*  RP (09/12)                                                *
*                                                            *
**************************************************************
      Implicit Real*8(A-H,O-Z)
c
      INTEGER NGrid
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC      
      REAL*8  F(NGrid),D1F(NGrid,7),RA(NGrid),RB(NGrid),
     $        D1RA(NGrid,3),D1RB(NGrid,3),TA(NGrid),TB(NGrid)
      Save One,Two,Three,Four,Five,Six,Seven,Eight,Nine
      Data One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Five/5.0d0/,
     $  Six/6.0d0/,Seven/7.0d0/,Eight/8.0d0/,Nine/9.0d0/
c
        dRA = 1
        dRB = 2
        dGA = 3
        dGB = 4
        dGC = 5
        dTA = 6
        dTB = 7
c
        G   = 0.004d+0
        ome = 2.5d+0
      If (ijzy.eq.1) then
c
! By Yiwei Liu
c       CF22D
c      
        CC000=0.244161168d0
        CC001=-0.389728151d0
        CC002=-1.829675858d0
        CC003=1.396044771d0
        CC004=2.315047133d0
        CC005=0.397552547d0
        CC010=1.082144406d0
        CC011=-7.894560034d0
        CC012=-3.656253030d0
        CC013=2.574496508d0
        CC014=4.031038406d0
        CC020=-3.931389433d0
        CC021=0.333519075d0
        CC022=-3.032270318d0
        CC023=3.673752289d0
        CC030=3.005997956d0
        CC031=-6.463733874d0
        CC032=-4.596755225d0
        CC100=0.964839180d0
        CC101=0.363791944d0
        CC102=1.646506623d0
        CC103=-3.504641550d0
        CC104=-3.922228074d0
        CC110=0.843718076d0
        CC111=10.779373313d0
        CC112=2.293612669d0
        CC113=7.088363286d0
        CC120=2.598770741d0
        CC121=-0.088522116d0
        CC122=7.18080903d0
        CC200=-1.017514009d0
        CC201=1.735020310d0
        CC202=3.499241561d0
        CC203=0.922224945d0
        CC210=-2.21290392d0
        CC211=0.243080429d0
        CC212=17.306321840d0
        CC300=0.311402396d0
        CC301=-3.257126009d0
        CC302=-3.372399742d0

      End If
      Tiny = 1.0d-8
      F12   = Two * Six
      F24   = Four * Six
      F28   = Four * Seven
      F2o3  = Two / Three
      F3o2  = Three / Two
      F1o3  = One / Three
      F4o3  = Four / Three
      F7o3  = Seven / Three
      F8o3  = Eight / Three
      F10o3 = F2o3 * Five
      F28o9 = F28 / Nine
      PI    = Four * ATan(One)
C
C     Local Spin Density factor.
C
      Ax = -F3o2*(F4o3*PI)**(-F1o3)
      Do 10 iG = 1, NGrid
c
c alpha component
c
        If(RA(i).gt.Tiny.and.TA(i).gt.Tiny) then
          F3o5 = Three/Five
          F6   = Six
          F5o3 = Five/Three
          pX    = RA(i)
          GamX2 = D1RA(i)
          S2    = GamX2*pX**(-F8o3)
          U     = G*S2/(One+G*S2)
          E = Ax*pX**F4o3
          Tau    = TA(i)
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3) 
          Tsig =TauUEG/Tau
          Wsig =(Tsig - One)/(Tsig + One)
c          
          FU = U
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FT = Wsig          
c
          FCF22D = CC000 + FT*CC001 + FT**2*CC002 + FT**3*CC003 + 
     $     FT**4*CC004 + FT**5*CC005 + 
     $    FU*CC010 + FT*FU*CC011 + FT**2*FU*CC012 + FT**3*FU*CC013 + 
     $     FT**4*FU*CC014 + FU**2*CC020 + 
     $    FT*FU**2*CC021 + FT**2*FU**2*CC022 + FT**3*FU**2*CC023 + 
     $     FU**3*CC030 + FT*FU**3*CC031 + 
     $    FT**2*FU**3*CC032 + FV*CC100 + FV*FT*CC101 + 
     $     FV*FT**2*CC102 + FV*FT**3*CC103 + 
     $    FV*FT**4*CC104 + FV*FU*CC110 + FV*FT*FU*CC111 + 
     $     FV*FT**2*FU*CC112 + FV*FT**3*FU*CC113 + 
     $    FV*FU**2*CC120 + FV*FT*FU**2*CC121 + FV*FT**2*FU**2*CC122 + 
     $     FV**2*CC200 + FV**2*FT*CC201 + 
     $    FV**2*FT**2*CC202 + FV**2*FT**3*CC203 + FV**2*FU*CC210 + 
     $     FV**2*FT*FU*CC211 + 
     $    FV**2*FT**2*FU*CC212 + FV**3*CC300 + FV**3*FT*CC301 + 
     $     FV**3*FT**2*CC302
c
          F(i) = F(i) + E*FCF22D
c
          If(MaxDrF.ge.1) then
            ER = F4o3*E/pX
            S    = Sqrt(S2)
            GamX = Sqrt(GamX2)
            SR   = -F4o3*S/pX
            SG   = S/GamX
            US   = Two*G*S/((One+G*S*S)**2)
c         
            dFVdR =  (ome/(Three*pX**F2o3))
     $                *(One+ome*pX**F1o3)**(-Two)
c
            dWdT = Two/((One + Tsig)**2)
            dTdR = ((F6*PI*PI)**F2o3)*(pX**F2o3)/Tau
            dTdTau = -TauUEG/Tau**2
c
            dFTdR = dWdT*dTdR
            dFTdTau=dWdT*dTdTau
c     
          dFCF22DdFV = CC100 + FT*CC101 + FT**2*CC102 + FT**3*CC103 + 
     $     FT**4*CC104 + FU*CC110 + 
     $    FT*FU*CC111 + FT**2*FU*CC112 + FT**3*FU*CC113 + 
     $     FU**2*CC120 + FT*FU**2*CC121 + 
     $    FT**2*FU**2*CC122 + Two*FV*CC200 + Two*FV*FT*CC201 + 
     $     Two*FV*FT**2*CC202 + Two*FV*FT**3*CC203 + 
     $    Two*FV*FU*CC210 + Two*FV*FT*FU*CC211+Two*FV*FT**2*FU*CC212+ 
     $     Three*FV**2*CC300 + Three*FV**2*FT*CC301 + 
     $    Three*FV**2*FT**2*CC302
c
          dFCF22DdFU = CC010 + FT*CC011 + FT**2*CC012 + FT**3*CC013 + 
     $     FT**4*CC014 + Two*FU*CC020 + 
     $    Two*FT*FU*CC021 + Two*FT**2*FU*CC022 + Two*FT**3*FU*CC023 + 
     $     Three*FU**2*CC030 + Three*FT*FU**2*CC031 + 
     $    Three*FT**2*FU**2*CC032 + FV*CC110 + FV*FT*CC111 + 
     $     FV*FT**2*CC112 + FV*FT**3*CC113 + 
     $    Two*FV*FU*CC120 + Two*FV*FT*FU*CC121+Two*FV*FT**2*FU*CC122+ 
     $     FV**2*CC210 + FV**2*FT*CC211 + 
     $    FV**2*FT**2*CC212
c    
          dFCF22DdFT = CC001 + Two*FT*CC002 + Three*FT**2*CC003 + 
     $     Four*FT**3*CC004 + Five*FT**4*CC005 + FU*CC011 + 
     $    Two*FT*FU*CC012 + Three*FT**2*FU*CC013+Four*FT**3*FU*CC014+ 
     $     FU**2*CC021 + Two*FT*FU**2*CC022 + 
     $    Three*FT**2*FU**2*CC023 + FU**3*CC031 + Two*FT*FU**3*CC032 + 
     $     FV*CC101 + Two*FV*FT*CC102 + 
     $    Three*FV*FT**2*CC103 + Four*FV*FT**3*CC104 + FV*FU*CC111 + 
     $     Two*FV*FT*FU*CC112 + Three*FV*FT**2*FU*CC113 + 
     $    FV*FU**2*CC121 + Two*FV*FT*FU**2*CC122 + FV**2*CC201 + 
     $     Two*FV**2*FT*CC202 + 
     $    Three*FV**2*FT**2*CC203 + FV**2*FU*CC211 + 
     $     Two*FV**2*FT*FU*CC212 + FV**3*CC301 + Two*FV**3*FT*CC302
c
        dFCF22DdR = dFCF22DdFV*dFVdR+dFCF22DdFU*US*SR+dFCF22DdFT*dFTdR
        dFCF22DdG = dFCF22DdFU*US*SG
        dFCF22DdT = dFCF22DdFT*dFTdTau
c
            D1F(i,dRA) = D1F(i,dRA) + ER*FCF22D
     $                      + E*dFCF22DdR
            D1F(i,dGA) = D1F(i,dGA) + E*dFCF22DdG/(Two*GamX)
            D1F(i,dTA) = D1F(i,dTA) + E*dFCF22DdT
            endIf
          endIf
c
c beta component
c
        If(RB(i).gt.Tiny.and.TB(i).gt.Tiny) then
          F3o5 = Three/Five
          F6   = Six
          F5o3 = Five/Three
          pX    = RB(i)
          GamX2 = D1RB(i)
          S2    = GamX2*pX**(-F8o3)
          U     = G*S2/(One+G*S2)
          E = Ax*pX**F4o3
          Tau    = TB(i)
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3) 
          Tsig =TauUEG/Tau
          Wsig =(Tsig - One)/(Tsig + One)
c          
          FU = U
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FT = Wsig          
c
          FCF22D = CC000 + FT*CC001 + FT**2*CC002 + FT**3*CC003 + 
     $     FT**4*CC004 + FT**5*CC005 + 
     $    FU*CC010 + FT*FU*CC011 + FT**2*FU*CC012 + FT**3*FU*CC013 + 
     $     FT**4*FU*CC014 + FU**2*CC020 + 
     $    FT*FU**2*CC021 + FT**2*FU**2*CC022 + FT**3*FU**2*CC023 + 
     $     FU**3*CC030 + FT*FU**3*CC031 + 
     $    FT**2*FU**3*CC032 + FV*CC100 + FV*FT*CC101 + 
     $     FV*FT**2*CC102 + FV*FT**3*CC103 + 
     $    FV*FT**4*CC104 + FV*FU*CC110 + FV*FT*FU*CC111 + 
     $     FV*FT**2*FU*CC112 + FV*FT**3*FU*CC113 + 
     $    FV*FU**2*CC120 + FV*FT*FU**2*CC121 + FV*FT**2*FU**2*CC122 + 
     $     FV**2*CC200 + FV**2*FT*CC201 + 
     $    FV**2*FT**2*CC202 + FV**2*FT**3*CC203 + FV**2*FU*CC210 + 
     $     FV**2*FT*FU*CC211 + 
     $    FV**2*FT**2*FU*CC212 + FV**3*CC300 + FV**3*FT*CC301 + 
     $     FV**3*FT**2*CC302
c
          F(i) = F(i) + E*FCF22D
c
          If(MaxDrF.ge.1) then
            ER = F4o3*E/pX
            S    = Sqrt(S2)
            GamX = Sqrt(GamX2)
            SR   = -F4o3*S/pX
            SG   = S/GamX
            US   = Two*G*S/((One+G*S*S)**2)
c         
            dFVdR =  (ome/(Three*pX**F2o3))
     $                *(One+ome*pX**F1o3)**(-Two)
c
            dWdT = Two/((One + Tsig)**2)
            dTdR = ((F6*PI*PI)**F2o3)*(pX**F2o3)/Tau
            dTdTau = -TauUEG/Tau**2
c
            dFTdR = dWdT*dTdR
            dFTdTau=dWdT*dTdTau
c     
          dFCF22DdFV = CC100 + FT*CC101 + FT**2*CC102 + FT**3*CC103 + 
     $     FT**4*CC104 + FU*CC110 + 
     $    FT*FU*CC111 + FT**2*FU*CC112 + FT**3*FU*CC113 + 
     $     FU**2*CC120 + FT*FU**2*CC121 + 
     $    FT**2*FU**2*CC122 + Two*FV*CC200 + Two*FV*FT*CC201 + 
     $     Two*FV*FT**2*CC202 + Two*FV*FT**3*CC203 + 
     $    Two*FV*FU*CC210 + Two*FV*FT*FU*CC211+Two*FV*FT**2*FU*CC212+ 
     $     Three*FV**2*CC300 + Three*FV**2*FT*CC301 + 
     $    Three*FV**2*FT**2*CC302
c
          dFCF22DdFU = CC010 + FT*CC011 + FT**2*CC012 + FT**3*CC013 + 
     $     FT**4*CC014 + Two*FU*CC020 + 
     $    Two*FT*FU*CC021 + Two*FT**2*FU*CC022 + Two*FT**3*FU*CC023 + 
     $     Three*FU**2*CC030 + Three*FT*FU**2*CC031 + 
     $    Three*FT**2*FU**2*CC032 + FV*CC110 + FV*FT*CC111 + 
     $     FV*FT**2*CC112 + FV*FT**3*CC113 + 
     $    Two*FV*FU*CC120 + Two*FV*FT*FU*CC121+Two*FV*FT**2*FU*CC122+ 
     $     FV**2*CC210 + FV**2*FT*CC211 + 
     $    FV**2*FT**2*CC212
c    
          dFCF22DdFT = CC001 + Two*FT*CC002 + Three*FT**2*CC003 + 
     $     Four*FT**3*CC004 + Five*FT**4*CC005 + FU*CC011 + 
     $    Two*FT*FU*CC012 + Three*FT**2*FU*CC013+Four*FT**3*FU*CC014+ 
     $     FU**2*CC021 + Two*FT*FU**2*CC022 + 
     $    Three*FT**2*FU**2*CC023 + FU**3*CC031 + Two*FT*FU**3*CC032 + 
     $     FV*CC101 + Two*FV*FT*CC102 + 
     $    Three*FV*FT**2*CC103 + Four*FV*FT**3*CC104 + FV*FU*CC111 + 
     $     Two*FV*FT*FU*CC112 + Three*FV*FT**2*FU*CC113 + 
     $    FV*FU**2*CC121 + Two*FV*FT*FU**2*CC122 + FV**2*CC201 + 
     $     Two*FV**2*FT*CC202 + 
     $    Three*FV**2*FT**2*CC203 + FV**2*FU*CC211 + 
     $     Two*FV**2*FT*FU*CC212 + FV**3*CC301 + Two*FV**3*FT*CC302
c
        dFCF22DdR = dFCF22DdFV*dFVdR+dFCF22DdFU*US*SR+dFCF22DdFT*dFTdR
        dFCF22DdG = dFCF22DdFU*US*SG
        dFCF22DdT = dFCF22DdFT*dFTdTau
c
            D1F(i,dRB) = D1F(i,dRB) + ER*FCF22D
     $                      + E*dFCF22DdR
            D1F(i,dGB) = D1F(i,dGB) + E*dFCF22DdG/(Two*GamX)
            D1F(i,dTB) = D1F(i,dTB) + E*dFCF22DdT
            endIf
          endIf
   10   Continue
      Return
      End
