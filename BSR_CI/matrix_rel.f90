!======================================================================
      Subroutine matrix_rel
!======================================================================
!     contribution of rel.integrals
!======================================================================
      Use  bsr_ci
      Use  c_data
      Use  z_core

      Implicit real(8) (A-H,O-Z)

      open(nur,form='UNFORMATTED',status='SCRATCH')

! ... s-o interaction:

      if(iso.gt.0) then

      icase = 7

      Call Gen_Zval

      Call Alloc_c_data(1,0,0,mblock,nblock,kblock,eps_c)
      Call Read_int_bnk(icase)
      Call Add_matrix(1,0)

      Do j=1,ncdata; i=IPT(j)
       ic = k3(i);  jc = k4(i)
       ii = k1(i);  jj = k2(i)
        write(nur) cdata(i)*Z_int(ii,jj), ic, jc, icase
      End do

      end if  ! over iso

! ... s-o-o inreraction:

      if(isoo.gt.0) then
       
      Call Alloc_c_data(1,0,kmax,mblock,nblock,kblock,eps_c)

      icase = 8
      Call Read_int_bnk(icase)
      Do k=0,kmax
       Call Add_matrix(1,k)
       Call I_data(icase,k) 
      End do
                                                                                            
      icase = 9
      Call Read_int_bnk(icase)
      Do k=0,kmax
       Call Add_matrix(1,icase)
       Call I_data(icase,k) 
      End do

      end if  ! over isoo

! ... s-s inreraction:

      if(iss.gt.0) then
       
      Call Alloc_c_data(1,0,kmax,mblock,nblock,kblock,eps_c)

      icase = 10
      Call Read_int_bnk(icase)
      Do k=0,kmax
       Call Add_matrix(1,k)
       Call I_data(icase,k) 
      End do

      end if  ! over isoo

      End Subroutine matrix_rel