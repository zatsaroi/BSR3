!======================================================================
      Subroutine Check_mat(met)
!======================================================================
!     check the overlap matrix for big overlaps and assign new
!     orthogonality conditions if needed
!----------------------------------------------------------------------
      Use bsr_mat,     only: S_ovl, S_pert, pri, nuc
      Use target,      only: BFT
      Use conf_LS,     only: WC
      Use channel    
      Use phys_orb_LS
      Use bsr_matrix
      Use spline_param
      Use spline_orbitals

      Implicit none
      Integer :: met, i,j,ij,i1,i2,j1,j2,ich,jch,ii,jj,it,jt,is,js,kp
      Real(8) :: S, v(ns)

      met = 0
      write(pri,'(/a,f6.3)') &
       'Checking the overlap matrix for overlaps > S_ovl =', S_ovl
!----------------------------------------------------------------------      
! ... analize the overlap matrix:

      Do ich=1,nch; it = iptar(ich) 
       i1=1; if(it.gt.1) i1=ip_tar(it-1)+1; i2=ip_tar(it)
      Do jch=1,ich; jt = iptar(jch)
       j1=1; if(jt.gt.1) j1=ip_tar(jt-1)+1; j2=ip_tar(jt)
       
       if(ich.eq.jch) Cycle
       
       ij=ich*(ich-1)/2+jch

       Do ii=i1,i2; i=ip_phy(ii); is=ip_sub(ii)  
       Do jj=j1,j2; j=ip_phy(jj); js=ip_sub(jj)

        if(lch(ich).ne.lbs(j)) Cycle
        if(lch(jch).ne.lbs(i)) Cycle

        v(1:ns) = matmul (hcc(1:ns,1:ns,ij),pbs(1:ns,i))
        S = SUM(v(:)*pbs(:,j)); S = abs(S)

        i = ipch(ich)
        if(S.gt.S_ovl.and.IBORT(i,js).ne.0) then
         write(pri,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,a6,2f10.3)') &
         '<',elc(ich),'|',ebs(js),'>=0', &
         BFT(it),elc(ich), BFT(jt),elc(jch),S,S_ovl
         write(nuc,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,a6,2f10.3)') &
         '<',elc(ich),'|',ebs(js),'>=0', &
         BFT(it),elc(ich), BFT(jt),elc(jch),S,S_ovl
         met = met + 1
         IBORT(i,js) = 0; IBORT(js,i)=0
        end if

        End do; End do

       End do ! over jch

       if(ncp.eq.0) Cycle

       Do kp=1,npert; v(1:ns)=hcb(1:ns,ich,kp)
        
        Do ii=1,nphys_sub; is=jp_sub(ii)  
         if(lch(ich).ne.lbs(is)) Cycle
         S = SUM(v(:)*pbs(:,is)); S = abs(S)
         if(S.gt.S_ovl) &
         write(pri,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,i6,f10.3)') &
         '<',elc(ich),'|',ebs(is),'>=0', &
         BFT(it),elc(ich), 'perturber   ',kp,S

         i=ipch(ich)
         if(S.gt.S_ovl.and.IBORT(i,is).ne.0) then
          write(nuc,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,i6,f10.3)') &
          '<',elc(ich),'|',ebs(is),'>=0', &
          BFT(it),elc(ich), 'perturber   ',kp,S
          IBORT(i,is) = 0; IBORT(is,i)=0
          met = met + 1
         end if
        End do ! over phys.orb

       End do ! over perturbers

      End do ! over ich

! ... check perturbers: 

      if(ncp.eq.0) Return

       Do i=1,npert; Do j=1,i; ij=i*(i-1)/2+j
        S = abs(hbb(ij))
        if(S.gt.S_pert.and.i.ne.j) then
         write(pri,'(f10.5,2i5,a)') &
           S, i,j , ' - suspicious perturber overlap '
         is = ippert(i)-ippert(i-1)          
         js = ippert(j)-ippert(j-1)          
         ii=is; if(js.lt.is) ii=js
         i1=ippert(ii-1)+1+ipconf(nch)
         i2=ippert(ii)+ipconf(nch)         
         WC(i1:i2) = 0.d0
         write(pri,'(a,i5,a)') 'pertuber',ii,'  was removed !!!'
         met=met+1
        end if
        if(S.lt.0.999.and.i.eq.j) then
         write(pri,'(f10.5,2i5,a)') &
          S, i,j , ' - suspicious perturber normalization '
        end if
       End do; End do

      End Subroutine Check_mat
