!======================================================================
      Subroutine Gen_matrix(jtype,jpol)
!======================================================================
!     merge data and generate the interaction matrix for current itype
!----------------------------------------------------------------------
      Use bsr_mat
      Use cmdata, nc => ncdata

      Implicit none
      Integer, intent(in) :: jtype,jpol
      Integer :: i,j,k,n
      Real(8) :: t1,t2

      Call CPU_time(t1)

! ... prepare the data:
       
      n = nblk(jpol,jtype)          ! number of blocks
      i = kblk(jpol,jtype,1)
      if(n.eq.1.and.jpblk(i).le.0) n = 0

      if(n.eq.0) then          ! nothing to do
       Return  

      elseif(n.eq.1) then      ! simple case
       nc = jpblk(i)-ipblk(i) + 1
       j = ipblk(i)-1
       Do i = 1,nc;  IPT(i)=j+i;  End do

      else                      ! need merge data from 
                                ! different blocks
       Do i = 1,n               
        j = kblk(jpol,jtype,i)
        ipi(i) = ipblk(j)  
        ipj(i) = jpblk(j)
       End do
       Call Merge_cdata(n, ipi, ipj, nc, EPS_c)

! ... release the blocks:       

       Do i=1,n; jpblk(kblk(jpol,jtype,i))=-1;  End do

      end if

! ... reassign the first block for current itype:

      i = kblk(jpol,jtype,1); j=jtype; k=jpol
      iblk(k,j)=i; nblk(k,j)=1; jpblk(i)=0

! ... generate matrix:

      Select case (icase)
       Case(1,2,3,4,5,8,9,10); Call I_data(jtype,jpol) 
       Case(6);                Call L_data(jtype,jpol)  
       Case(7);                Call Z_data(jtype,jpol)  
       Case(11);               Call O_data
      End Select

      Call CPU_time(t2)

      if(debug.gt.1) &
      write(pri,'(a,3(a,i2),a,i3,a,i6,a,f8.2,a)') 'Gen_matrix:', &
        '  icase=',icase, '  itype=',jtype, '  kpol=',jpol, &
        '  n =',n , '  nc=',nc,'  time=',t2-t1,' sec'

      End Subroutine Gen_matrix
