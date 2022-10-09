      Subroutine CF22Dc(F,D1F,RA,RB,D1RA,D1RB,TA,TB,NGrid,ijzy)
************************************************************************
*                                                                      *
*  CF22Dc evaluates the correlation part of the CF22D functional on the*
*  the grid.                                                           *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  OUTPUT:                                                             *
*     F      - Functional values                                       *
*     D1F    - First derivatives with respect to RA, RB, GA, GB        *
*              TA, TB                                                  *
*                                                                      *
*  INPUT:                                                              *
*                                                                      *
*       ijzy - 1 CF22D                                                 *
*                                                                      *
*     RA,B   - Spin densities                                          *
*     D1RA,B - Spin density gradients                                  *
*     TA,B   - Spin kinetic energy densities                           *
*     NGrid  - number of grid points                                   *
*                                                                      *
*  RP (09/12), YZ (12/08)                                              *
*                                                                      *
************************************************************************
      Implicit Real*8(A-H,O-Z)
      Real*8 LSDA
      INTEGER NGrid
      REAL*8  F(NGrid),D1F(NGrid,7),RA(NGrid),RB(NGrid),
     $        D1RA(NGrid,3),D1RB(NGrid,3),TA(NGrid),TB(NGrid)
      Integer dRA, dRB, dTA, dTB, dGA, dGB, dGC
      Save F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11
      DATA F1/1.0D+00/,  F2/2.0D+00/,  F3/3.0D+00/,
     $     F4/4.0D+00/,  F5/5.0D+00/,  F6/6.0D+00/,
     $     F7/7.0D+00/,  F8/8.0D+00/,  F9/9.0D+00/, 
     $     F10/10.0D+00/,F11/11.0D+00/
      parameter( pi = 3.1415926535897932384626433832795d0 )
      dRA = 1
      dRB = 2
      dGA = 3
      dGB = 4
      dGc = 5
      dTA = 6
      dTB = 7
C
      DTol   = 1.0d-8
C
      F1o3 = F1/F3
      F2o3 = F2/F3
      F5o3 = F5/F3
C
      Pi34   = F3/(F4*Pi)
C
      if (ijzy.eq.1) then
! By Yiwei Liu
C     Parameters for CF22D
        at00=0.873863376D+00
        at01=0.078066142D+00
        at02=6.576550257D+00
        at03=-1.126030147D+00
        at04=-3.244797887D+00
        at05=-2.186090839D+00
        at06=-3.489135041D+00
        at07=3.090689716D+00
        at08=3.866592474D+00
        at09=  0.0000000D+00
        at10=  0.0000000D+00
        at11=  0.0000000D+00
C           
        bt00=0.828203832D+00
        bt01=-2.518707202D+00
        bt02=10.436806314D+00
        bt03=3.588267084D+00
        bt04=-5.789404145D+00
        bt05=3.353560215D+00
        bt06=-2.432384384D+00
        bt07=-1.147183331D+00
        bt08=2.991316045D+00
        bt09=  0.0000000D+00
        bt10=  0.0000000D+00
        bt11=  0.0000000D+00
C
       endif
      
      DO i = 1,NGrid
       RhoA = RA(i)
       RhoB = RB(i)
       Rho = RhoA + RhoB
       TauA = TA(i)/F2
       TauB = TB(i)/F2
       Tau = TauA + TauB

       If(Rho.gt.DTol.and.Tau.gt.DTol) then
        RS = (Pi34/Rho)**F1o3
        Zeta = (RhoA-RhoB)/Rho
        TauUEG=F3*(F3*Pi*Pi)**(F2o3)*Rho**(F5o3)/F10
        Tsig =TauUEG/Tau
        Wsig =(Tsig - F1)/(Tsig + F1)
        Fsig1=(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))

        Fsig2=(bt0 + Wsig*(bt1 + Wsig*(bt2 + Wsig*(bt3 + Wsig*(
     &            bt4 + Wsig*(bt5 + Wsig*(bt6 + Wsig*(bt7 + Wsig*(
     &            bt8 + Wsig*(bt9 + Wsig*(bt10+Wsig*bt11)))))))))))

        Y = (D1RA(i,1) + D1RB(i,1))**F2
     $      + (D1RA(i,2) + D1RB(i,2))**F2
     $      + (D1RA(i,3) + D1RB(i,3))**F2
        GRho = Sqrt(Y)
c       
c      lsdac is a subroutine to evaluate the Perdew-Wang-91 correlation functional 
c      local spin density approximation (LSDA) to the correlation energy of a uniform 
c      electron gas. (Phys. Rev. B 45, 13244 (1992)). Users should provid their own
c      for this LSDA correlation functional or they may find this routine on Kieron 
c      Burke's Web site at http://www.chem.uci.edu/~kieron/dftold2/pubs/PBE.asc
c
        Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ)
        LSDA = Rho*PotLC
c
c      PBEH0 is a subroutine to evaluate the H0 term in the PBE correlation functional
c      (Phys. Rev. Lett. 77, 3865 - 3868 (1996)) Users should provid their own
c      for this H0 subroutine or they may find this routine on Kieron
c      Burke's Web site at http://www.chem.uci.edu/~kieron/dftold2/pubs/PBE.asc
c
        Call PBEH0(Rho,GRho,RS,Zeta,PotLC,dLdS,dLdZ,H,dHdR,dHdG,dHdZ)
        GGA = Rho*H 
        E1 = LSDA*Fsig1
        E2 = GGA*Fsig2
        F(i) = F(i)+ E1 +E2
        
c
c     functional derivatives
c
         RSP = -RS/(F3*Rho)
         dZdA = (F1-Zeta)/Rho
         dZdB = (-F1-Zeta)/Rho
         dLdRA = dLdS*RSP + dLdZ*dZdA
         dLdRB = dLdS*RSP + dLdZ*dZdB
         dF1dW=( at1 + Wsig*(F2  *at2 + Wsig*(F3*at3 + Wsig*(
     &            F4 *at4 + Wsig*(F5 *at5 + Wsig*(F6  *at6 + Wsig*(
     &            F7*at7 + Wsig*(F8*at8 + Wsig*(F9 *at9 + Wsig*(
     &            F10  *at10+ Wsig*F11*at11))))))))))
         dF2dW=( bt1 + Wsig*(F2  *bt2 + Wsig*(F3*bt3 + Wsig*(
     &            F4 *bt4 + Wsig*(F5 *bt5 + Wsig*(F6  *bt6 + Wsig*(
     &            F7*bt7 + Wsig*(F8*bt8 + Wsig*(F9 *bt9 + Wsig*(
     &            F10  *bt10+ Wsig*F11*bt11))))))))))
         dWdT = F2/((F1 + Tsig)**F2)
         dTdR = Tsig*F5/(F3*Rho) 
         dTdTau = -Tsig/Tau
         dF1dR = dF1dW*dWdT*dTdR
         dF1dTau=dF1dW*dWdT*dTdTau
         dF2dR = dF2dW*dWdT*dTdR
         dF2dTau=dF2dW*dWdT*dTdTau
         dLDdRA = PotLC + Rho*dLdRA
         dLDdRB = PotLC + Rho*dLdRB
         dHdRA = dHdR + dHdZ*dZdA
         dHdRB = dHdR + dHdZ*dZdB
         dGRhodY = F1/(F2*GRho)
         dHdY = dHdG * dGRhodY
         dHdGA = dHdY
         dHdGC = dHdY*F2  
         dGGAdRA = H + Rho*dHdRA
         dGGAdRB = H + Rho*dHdRB
         dGGAdGA = Rho*dHdGA
         dGGAdGB = dGGAdGA
         dGGAdGC = Rho*dHdGC
C
         dE1dRA = dLDdRA*Fsig1 + LSDA*dF1dR
         dE1dRB = dLDdRB*Fsig1 + LSDA*dF1dR
         dE1dKA = LSDA*dF1dTau
         dE1dKB = dE1dKA
C
         dE2dRA = dGGAdRA*Fsig2 + GGA*dF2dR
         dE2dRB = dGGAdRB*Fsig2 + GGA*dF2dR 
         dE2dKA = GGA*dF2dTau
         dE2dKB = dE2dKA
         dE2dGA = dGGAdGA*Fsig2
         dE2dGB = dGGAdGB*Fsig2
         dE2dGC = dGGAdGC*Fsig2   

         D1F(i,dRA)=D1F(i,dRA) + dE1dRA + dE2dRA 
         D1F(i,dRB)=D1F(i,dRB) + dE1dRB + dE2dRB  
         D1F(i,dTA)=D1F(i,dTA) + (dE1dKA + dE2dKA)/F2
         D1F(i,dTB)=D1F(i,dTB) + (dE1dKB + dE2dKB)/F2
         D1F(i,dGA)=D1F(i,dGA) + dE2dGA
         D1F(i,dGB)=D1F(i,dGB) + dE2dGB
C     GC is the dot product of the vectors D1RA and D1RB
         D1F(i,dGC)=D1F(i,dGC) + dE2dGC 
       Endif
      Enddo
      Return
      End
