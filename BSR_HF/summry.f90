!====================================================================== 
      Subroutine SUMMRY 
!====================================================================== 
!     The results of a calculation are summarized.   These include
!     the following for each electron:
!          E(NL)   - diagonal energy parameter
!          I(NL)   - -(1/2)<nl|L|nl>
!          KE      - I(NL) + Z <r>
!          REL     - Relativistic shift (mass-velocity, Darwin term,
!                    spin-spin contact term)
!          SIGMA   - screening parameter.
!          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0
!          1/R**3  - expected value of <1/r**3>
!          1/R     - expected value of <1/r>
!          R       - expected mean radius
!          R**2    - expected value of <r**2>
!
!     These results are followed by:
!
!          KINETIC ENERGY (EK)
!          POTENTIAL ENERGY (EP) = ETOTAL - EN
!          RATIO                 = EP/EN
!          TOTAL ENERGY (ETOTAL)
!
!----------------------------------------------------------------------
      Use bsr_hf
      Use hf_orbitals 

      Implicit none
      Integer :: i
      Real(8) :: ekinp, en, rom3, rp2, rz, epotl,ratio, rh
      Real(8) :: ro1(nbf),hli(nbf),rom1(nbf),AZ(nbf),S(nbf)
      Real(8), external ::  azl, bhl_hf, quadr_hf, rk

      write(log,'(//80("-"))')
      write(log,'(/24X,a,3x,a,5x,a,3x,a)') 'ATOM',ATOM,'TERM',TERM 

! ... COMPUTE AND PRiNT ONE-ELECTRON PARAMETERS

      write(log,'(//3x,a,7X,a,8X,a,8X,a,7X,a,5X,a,4X,a)') &
               'nl','E(nl)','I(nl)','KE(nl)','AZ(nl)','S(nl)','ns'

      EN = 0.d0
      Do i = 1,nbf 
       RO1(i)  = QUADR_hf(i,i,1) 
       HLi(i) = -0.5*bhl_hf(i,i)
       ROM1(i) = QUADR_hf(i,i,-1) 
       EKiNP  = HLi(i) + Z*ROM1(i) 
       EN     = EN + qsum(i)*EKiNP 
       RH     = 3*nbs(i)*nbs(i) - lbs(i)*(lbs(i)+1) 
       S(i)   = Z - 0.5*RH/RO1(i) 
       AZ(i)  = azl(z,h,ks,lbs(i)+1) * p(lbs(i)+2,i)
       write(log,'(1X,A4,4F13.6,F10.3,i6)') & 
             ebs(i),E(i,i),HLi(i),EKiNP,AZ(i),S(i),mbs(i) 
      End do 

! ...  Compute Moments

      write(log,'(//3x,a,7X,a,8X,a,11X,a,10X,a,6X,a,2x,a)') &
           'nl','1/R**3','1/R','R','R**2','Delta(r)','max_R'
      Do i = 1,nbf 
       ROM3 = 0.d0;  if (lbs(i) /= 0) ROM3 = QUADR_hf(i,i,-3) 
       RP2 = QUADR_hf(i,i,2) 
       RZ = 0.D0;   if (lbs(i) == 0) RZ = AZ(i)**2/(4.*Pi)
       write(log,'(1X,A4,4F13.6,F10.3,F10.3)') &
             ebs(i), ROM3, ROM1(i), RO1(i), RP2, RZ, t(mbs(i)+ks)  
      End do 
      
      EPOTL = ETOTAL - EN 
      RATiO = EPOTL/EN 
      Call core_energy
      write(log,'(//5X,a,5X,F16.8)') 'CORE ENERGY (a.u.)', ECORE
      write(log,'(//5X,a,5X,F16.8)') 'TOTAL ENERGY (a.u.)', ETOTAL

      rel = .true.
      lh=-2;  meth=''; krk=-1
      Call update_int(0)
      Call energy                                                      
      write(log,'(5X,a,5X,F16.8)')   'REL ENERGY         ', ETOTAL

      write(log,'(T20,a,F16.8)') 'Kinetic   ', EN
      write(log,'(T20,a,F16.8)') 'Potential ', EPOTL
      write(log,'(T20,a,F16.8)') 'Ratio     ', RATiO 

      End Subroutine SUMMRY 

