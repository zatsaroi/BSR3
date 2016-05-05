!======================================================================
      Subroutine J_case
!======================================================================
!     calculate f-values in LSJ mode
!----------------------------------------------------------------------
      Use bsr_dmat; Use dmatrix; Use term_LS

      Implicit none
      Character(80) :: AS
      Character(64) :: Label1,Label2
      Character( 5) :: AA
      Integer :: i,j,ii,jj,j1,j2,kp,met1,met2, KLT1,KLT2
      Real(8) :: g1,g2,E1,E2, de, SL,SV, S,C, alfaL,alfaV, alfL,alfV
      Real(8) :: FL,FV,WL,WV,GFL,GFV,GWL,GWV,dmatL,dmatV, ANGS, ANGSA
      Integer, external :: ITRI
      Real(8), external :: ANGS_AIR,DJ_fact,DJM_fact

!----------------------------------------------------------------------
      if(jmode.eq.2) then
       met1=jot1; met2=jot2
      else
       met1=0;    met2=0
      end if

      Deallocate(C1); Allocate(C1(kdm1))
      Deallocate(C2); Allocate(C2(kdm2))
!----------------------------------------------------------------------
!... loop over initial set:

      rewind(inb1);  read(inb1,'(a)') AS
      Do j1=1,mstate1
       Call Read_sol(ctype1,inb1,kdm1,C1,Label1,E1,jot1)

       if(istate1.gt.0.and.istate1.ne.j1) Cycle

       alfaL=0.0; alfaV=0.0

!... loop over final set:

      rewind(inb2);  read(inb2,'(a)') AS
      Do j2=1,mstate2
       Call Read_sol(ctype2,inb2,kdm2,C2,Label2,E2,jot2)

       if(istate2.gt.0.and.istate2.ne.j2) Cycle

!----------------------------------------------------------------------
! ... calculate J-factor when jmode <> 2:

      if(ITRI(jot1,jot2,kpol+kpol+1).eq.0) Cycle

      if(kpol.gt.0.and.(met1.ne.jot1.or.met2.ne.jot2)) then

       Do i=1,nterm1; Do j=1,nterm2
        DJ(i,j)=DJ_fact(ILterm1(i),ILterm2(j),ISterm1(i),ISterm2(j),JOT1,JOT2,kpol)
       End do; End do
       if(ktype.eq.'M') then
        Do i=1,nterm1; Do j=1,nterm2
         DJM(i,j)= &
          DJM_fact(ILterm1(i),ILterm2(j),ISterm1(i),ISterm2(j),JOT1,JOT2,kpol)
        End do; End do
       else
	      DJM = DJ
       end if
       met1=jot1; met2=jot2

      end if
!----------------------------------------------------------------------
      de = abs(E2 - E1); if(de.lt.2.d-8) Cycle

      SL=0.d0; SV=0.d0
      if(jmode.eq.1) then

       Do i=1,kdm1; ii=LSP1(i) 
       Do j=1,kdm2; jj=LSP2(j)
        SL = SL + C1(i)*DL(i,j)*C2(j)*DJ (ii,jj)
        SV = SV + C1(i)*DV(i,j)*C2(j)*DJM(ii,jj)
       End do; End do
      
      else

       SL=0.d0; SV=0.d0
       Do j=1,kdm2
        C = SUM(DL(1:kdm1,j)*C1(1:kdm1));  SL = SL + C2(j)*C
        C = SUM(DV(1:kdm1,j)*C1(1:kdm1));  SV = SV + C2(j)*C
       End do

      end if

      dmatL = SL
      dmatV = SV/de

      if(kpol.eq.0) then
       write(nur,'(/i4,f14.8,2x,a)') jot1-1,E1,Label1
       write(nur,'( i4,f14.8,2x,a)') jot2-1,E2,Label2
       write(nur,'(a,2E12.5,f10.2)') 'S = ',SL,SL*SL,SL*SL*100
       Cycle
      end if

      if(ktype.eq.'E') then
       SL = SL*SL; SV = SV*SV/(de*de)
       if(kpol.ne.1) SV=0.d0
      else
       SL = kpol*(kpol+kpol-1)*(SL/(kpol+1)+(ge/2.d0)*SV)**2; SV=0.d0
      end if

      if(SL.eq.0.0) Cycle

      kp = kpol+kpol+1
      S = 1.d0;  Do i = 1,kp,2;  S = S * i;  End do
      S = 2*kp*(kpol+1)/(S*S)/kpol 
      S = S * (de/c_au)**kp

      GWL = S*SL; GFL = c_au**3/2*GWL/(de*de)
      GWV = S*SV; GFV = c_au**3/2*GWV/(de*de)

      if(ktype.eq.'M') then
       GFL = GFL / c_au**2
       GWL = GWL / c_au**2
       SL = SL * 4    ! due to defition of corr. atomic unit
      end if

      g1 = jot1; g2 = jot2

      alfL = 2*dmatL*dmatL/(E2-E1)/(kp*g1); alfaL = alfaL + alfL   
      alfV = 2*dmatV*dmatV/(E2-E1)/(kp*g1); alfaV = alfaV + alfV

      KLT1 = IPT1*jot1 
      KLT2 = IPT2*jot2 

      if(E2.gt.E1) then
       FL = GFL/g1; WL = GWL/g2/time_au
       FV = GFV/g1; WV = GWV/g2/time_au
       write(nur,'(/i4,f14.8,2x,a)') KLT1,E1,trim(Label1) 
       write(nur,'( i4,f14.8,2x,a)') KLT2,E2,trim(Label2) 
      else
       FL = GFL/g2; WL = GWL/g1/time_au
       FV = GFV/g2; WV = GWV/g1/time_au
       write(nur,'(/i4,f14.8,2x,a)') KLT2,E2,trim(Label2) 
       write(nur,'( i4,f14.8,2x,a)') KLT1,E1,trim(Label1) 
      end if

      de = abs(E1-E2)*au_cm; ANGS=1.d+8/de; ANGSA=ANGS_AIR(ANGS)
      write(nur,'(f11.2,a5,2(3x,f10.2,a10))') &
                DE,' CM-1',ANGS,' ANGS(VAC)',ANGSA,' ANGS(AIR)'

      AA='FIK= '
      if(GF.ne.'f') then; FL=GFL; FV=GFV; AA='GF = '; end if

      if(ialfa.eq.0) then
       write(nur,'(1x,a1,i1,2x,a4,1PD12.5,3x,a5,D12.5,3x,a6,D12.5)') &
                   ktype,kpol,'S = ',SL,AA,FL,'AKI = ',WL
       if(SV.ne.0.d0) & 
       write(nur,'(9x,1PD12.5,8x,D12.5,9x,D12.5)') SV,FV,WV
      else
       write(nur,'(1x,a1,i1,2x,a4,1PD12.5,3x,a5,D12.5,3x,a6,D12.5, &
                   a8,d16.8,0P2f10.5)') &
                   ktype,kpol,'S = ',SL,AA,FL,'AKI = ',WL, &
                   '   RME =',dmatL, alfL,alfaL
       if(SV.ne.0.d0) & 
       write(nur,'(9x,1PD12.5,8x,D12.5,9x,D12.5,8x,d16.8,0P2f10.5)') &
                   SV,FV,WV, dmatV, alfV,alfaV 
      end if      

      End do; End do   !  over solutions

      End Subroutine J_case
 