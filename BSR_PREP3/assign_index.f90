!======================================================================
      Subroutine Assign_index(m)
!======================================================================
!     assign set index for orbital m
!----------------------------------------------------------------------
      Use bsr_prep

      Implicit real(8) (A-H,O-Z)
      Character(4) :: ELW
      Character(4), external :: ELF4

       kbs(m)=-1; ELW = EBS(m); S = 0.d0
       knew = New_index(lbs(m),ksmax,nbf,lbs,kbs)

       Do k = 1,knew-1
        S=0.d0           
        Do i = 1,nbf
         if(lbs(m).ne.lbs(i).or.k.ne.kbs(i)) Cycle
         S=max(S,abs(OBS(i,m)))
        End do
        if(S.lt.eps_ovl) then; kbs(m)=k; Exit; end if
       End do  

       if(kbs(m).eq.-1) then  ! the orbital belongs to new set  

        kbs(m) = knew
        EBS(m)=ELF4(nbs(m),lbs(m),kbs(m))
        write(pri,'(a,a,a,a,f15.8)') &
              elw,' --> ',EBS(m),'    new orbitals and new set index',S

       else                   ! check the same label for diff.orbitals            

        EBS(m)=ELF4(nbs(m),lbs(m),kbs(m))
        Do i = 1,nbf
         if(EBS(m).ne.EBS(i)) Cycle
         kbs(m) = knew
         EBS(m)=ELF4(nbs(m),lbs(m),kbs(m))
         write(pri,'(a,a,a,a,f15.8)') &
              elw,' --> ',EBS(m),'    new orbitals and new set index',S
        End do
        if(kbs(m).ne.knew) write(pri,'(a,a,a,a,f15.8)') &
              elw,' --> ',EBS(m),'    new orbitals but old set index',S
       
       end if

      End Subroutine Assign_index
