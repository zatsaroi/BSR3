!====================================================================== 
      Subroutine Summry
!====================================================================== 
!     The results of a calculation are summarized in "SUMMRY" file 
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer :: i, j, ipt(nlevels)
      Real(8) :: ar1,ar2,arm1, EN,EK,EPOTL,RATIO
      Real(8), external :: quadr, bhl

      write(log,'(/80(''-'')/)')
      write(log,'(a,a,a,a)') 'ATOM:  ',ATOM,'     TERM:  ',TERM 

      write(log,'(/a)') 'Convergence (latest difference):'
      write(log,'(a,T40,1Pd10.2)') 'Orbital diff. =', orb_diff
      write(log,'(a,T40,1Pd10.2)') 'SCF diff. =', scf_diff

! ... Compute Moments

      write(log,'(/2x,a,7x,a,8x,a,7x,a,7x,a,7x,a,7x,a,5x,a/)') &
           'nl','E(nl)','1/R','R','R**2','dmp','ns','max_R'

      EN = 0.d0
      Do i = 1, nbf 
       ar1  = quadr(i,i,1) 
       ar2  = quadr(i,i,2)
       arm1 = quadr(i,i,-1) 
       write(log,'(a5,f15.8,3f9.3,1PE12.2,i6,0Pf10.2)') &
        ebs(i), e(i), arm1, ar1, ar2, dpm(i), mbs(i), t(mbs(i)+ks)  
       EK =  -0.5*bhl(i,i) + Z*quadr(i,i,-1)
       EN = EN + qsum(i)*EK
      End do 

      EPOTL = ETOTAL - EN 
      RATIO = EPOTL/EN 
      
      write(log,'(/a,T20,f20.12)') 'Total energy', Etotal
      write(log,'( a,T20,f20.12)') 'Kinetic     ', EN
      write(log,'( a,T20,f20.12)') 'Potential   ', EPOTL
      write(log,'( a,T20,f20.12)') 'Ratio       ', RATIO 

      write(log,'(/a)') 'Optimized states:'

      write(log,'(/a/)') ' Level     Energy      Leading CSFs'

      Call SortR(nlevels,elevel,ipt)
      Do j = 1,nlevels; i = ipt(j)
       Call print_level(i,ncfg)
      End do

      Call Write_cc(ipt(1))

      End Subroutine Summry 


!======================================================================
      Subroutine print_level(i,nc)
!======================================================================
! ... print major contributors to ASF:
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer, intent(in) :: i,nc
      Integer :: ip(nc), j,ii
      Real(8) :: CP(nc)
     
      CP = -abs(coefs(ip_level(i)+1:ip_level(i)+nc)) 
      Call Sortr(nc,CP,ip)
      CP = coefs(ip_level(i)+1:ip_level(i)+nc) 
      ii = min(5,nc)      
      write(log,'(i4,f16.8,5(i4,f8.4))') &
        level(i),elevel(i),(ip(j),CP(ip(j)),j=1,ii) 

      End Subroutine print_level


!======================================================================
      Subroutine Write_cc(i)
!======================================================================
! ... print major contributors to ASF:
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer, intent(in) :: i
      Integer :: ip(ncfg), ic
      Real(8) :: CP(ncfg)
     
      CP = -abs(coefs(ip_level(i)+1:ip_level(i)+ncfg)) 
      Call Sortr(ncfg,CP,ip)
      CP = coefs(ip_level(i)+1:ip_level(i)+ncfg) 

      AF = trim(name)//'.cc'
      Open(nuc,file=AF)
      write(nuc,'(4x,a2,4x,a2,1x,F18.8)') atom,term,elevel(i)
      write(nuc,'(20a4)') (ebs(ic),ic=1,ncore)
      Do ic = 1,ncfg
       Call Pri_conf (nuc,ip(ic),CP(ip(ic)))
      End do

      End Subroutine Write_cc

