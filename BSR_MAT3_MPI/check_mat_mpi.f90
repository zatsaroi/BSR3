!======================================================================
      Subroutine Check_mat(met)
!======================================================================
          
      Use MPI

      USE bsr_mat
      USE bsr_matrix
      USE channel    
      USE target, only: BFT
      USE orb_LS, only: ief
      Use conf_LS
      USE phys_orb_LS
      USE spline_param
      USE spline_orbitals

      Implicit none
      
      Integer :: met, i,j,ij,i1,i2,j1,j2,ich,jch,ii,jj,it,jt,is,js,kp,io,jo
      Real(8) :: S, v(ns)
      Real(8), Allocatable :: Rarr(:,:)
      Integer, Allocatable :: Iarr(:,:)

      met = 0

      if(pri.gt.0) write(pri,'(/a,f8.3)') &
         'Check overlap matrix for big overlaps > s_ovl =',s_ovl

!----------------------------------------------------------------------      
! ... analize the overlap matrix:

      Do ich=1,nch; it = iptar(ich); io = ipch(ich) 
       i1=1; if(it.gt.1) i1=ip_tar(it-1)+1; i2=ip_tar(it)
      Do jch=1,ich; jt = iptar(jch); jo = ipch(jch) 
       j1=1; if(jt.gt.1) j1=ip_tar(jt-1)+1; j2=ip_tar(jt)
       if(ich.eq.jch) Cycle
       ij=icc(ich,jch); if(ij.eq.0) Cycle
       Do ii=i1,i2; i=ip_phy(ii); is=ip_sub(ii)  
       Do jj=j1,j2; j=ip_phy(jj); js=ip_sub(jj)

        if(lch(ich).ne.lbs(j)) Cycle
        if(lch(jch).ne.lbs(i)) Cycle
        v(1:ns) = matmul (hcc(1:ns,1:ns,ij),pbs(1:ns,i))
        S = SUM(v(:)*pbs(:,j)); S = abs(S)
        if(S.gt.S_ovl.and.IBORT(i,js).ne.0) then
         if(pri.gt.0) &
         write(pri,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,a6,2f10.3)') &
         '<',elc(ich),'|',ebs(js),'>=0', &
         BFT(it),elc(ich), BFT(jt),elc(jch),S,S_ovl
         met = met + 1
         IBORT(io,js)=-jch; IBORT(js,io)=-jch 
         OBS(io,js)=S; OBS(io,js)=S 
        end if

        End do; End do

       End do ! over jch

       if(ncp.eq.0) Cycle

       Do kp=1,npert
       
        i = icb(ich,kp); if(i.eq.0) Cycle; v(1:ns)=hcb(1:ns,i)
        
        Do ii=1,nphys_sub; is=jp_sub(ii)  
         if(lch(ich).ne.lbs(is)) Cycle
         S = SUM(v(:)*pbs(:,is)); S = abs(S)
         if(S.gt.S_ovl.and.pri.gt.0) &
         write(pri,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,i6,2f10.3)') &
         '<',elc(ich),'|',ebs(is),'>=0', &
         BFT(it),elc(ich), 'perturber   ',kp,S,S_ovl
         if(S.lt.S_ovl.or.IBORT(io,is).eq.0) Cycle
         j=nch+kp
         IBORT(io,is)=-j; IBORT(is,io)=-j; OBS(io,is)=S; OBS(io,is)=S 
         met = met + 1
        End do ! over phys.orb
       End do ! over perturbers

      End do ! over ich

!----------------------------------------------------------------------

      if(ncp.gt.0) then

        Do i=1,npert; Do j=1,i; ij=ibb(i,j); if(ij.eq.0) Cycle
         S = abs(hbb(ij))
         if(S.gt.S_pert.and.i.ne.j) then
          if(pri.gt.0) write(pri,'(f10.5,2i5,a)') &
          S, i,j , ' - suspicious perturber overlap '
          is = ippert(i)-ippert(i-1)          
          js = ippert(j)-ippert(j-1)          
          ii=is; if(js.lt.is) ii=js
          i1=ippert(ii-1)+1+ipconf(nch)
          i2=ippert(ii)+ipconf(nch)         
          WC(i1:i2) = 0.d0
          write(prj,'(a,i5,a)') 'pertuber',ii,'  was removed !!!'
          met=met+1
         end if
         if(S.lt.0.999.and.i.eq.j) then
          write(prj,'(f10.5,2i5,a)') &
          S, i,j , ' - suspicious perturber normalization '
         end if
        End do; End do

       end if

!----------------------------------------------------------------------

      Call MPI_REDUCE(met,i,1,MPI_integer,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) met=i

      Call MPI_BCAST(met,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(met.eq.0) Return

! ... collect IBORT:

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if(myid.eq.0) Allocate(Iarr(mbf,mbf),Rarr(mbf,mbf))

      Call MPI_REDUCE(IBORT,Iarr,mbf*mbf,MPI_INTEGER,MPI_MIN, &
                      0,MPI_COMM_WORLD,ierr)

      Call MPI_REDUCE(OBS,Rarr,mbf*mbf,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

! ... output new orthogonal constraints:

      if(myid.eq.0) then

      i=LEN_TRIM(AF_cfg); AF_cfg(i-2:i)=ALSP
      Open(nuc,file=AF_cfg,position='APPEND')

      Do i = 1,nbf; ich=ief(i); if(ich.le.0) Cycle; it=iptar(ich)
       Do j = 1,nbf; jch=-Iarr(i,j); if(jch.le.0) Cycle;
         if(Iarr(i,j).ge.0) Cycle
         S = Rarr(i,j)
         if(jch.le.nch) then; jt=iptar(jch)
           write(nuc,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,a6,2f10.3)') &
           '<',ebs(i),'|',ebs(j),'>=0', &
           BFT(it),elc(ich), BFT(jt),elc(jch),S,S_ovl
         else
          kp=jch-nch
          write(nuc,'(a1,a4,a1,a4,a3,5x,a12,a6,5x,a12,i6,2f10.3)') &
          '<',ebs(i),'|',ebs(j),'>=0', &
          BFT(it),elc(ich), 'perturber',kp,S,S_ovl
         end if
       End do
      End do

      Close(nuc)
             
! ... clean IBORT and OBS:

      Do i = 1,nbf
       Do j = 1,nbf
        if(Iarr(i,j).ge.0) Cycle
        IBORT(i,j)=0; IBORT(j,i)=0; OBS(i,j)=0.d0; OBS(i,j)=0.d0 
       End do
      End do
      Deallocate(Rarr,Iarr)
      end if

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_BCAST(OBS,mbf*mbf,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(IBORT,mbf*mbf,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      End Subroutine Check_mat
