!======================================================================
      Subroutine Record_is_done 
!======================================================================
      Use bsr_breit
      Use term_exp

      rewind(nur)
      write(nur) ic_case
      write(nur) IS_NEED

      End Subroutine Record_is_done 
           
!======================================================================
      Subroutine Read_is_done 
!======================================================================
      Use bsr_breit
      Use term_exp

      Implicit none
      Integer :: i,j, ij
      Integer, external :: DEF_ij

      ij=DEF_ij(ic_case,ic_case)
      if(allocated(IS_NEED)) Deallocate(IS_NEED)
      Allocate(IS_NEED(ij)); IS_NEED = 0

      if(allocated(JS_NEED)) Deallocate(JS_NEED)
      Allocate(JS_NEED(ic_case)); JS_NEED=0

      rewind(nur)
      read(nur) i
      if(i.ne.ic_case) Stop 'different ic_case in det_done file'

      read(nur) IS_NEED

      JS_NEED=0
      Do i=1,ic_case
       Do j=1,i
        ij = DEF_ij(i,j)
        if(IS_NEED(ij).eq.0) Cycle 
        JS_NEED(i)=1; Exit
       End do
      End do

      End Subroutine Read_is_done 

           
!======================================================================
      Subroutine Define_is_done 
!======================================================================
      Use bsr_breit
      Use term_exp
      Use conf_LS,      only: ne
      Use symc_list_LS, only: JC_need

      Implicit none
      Integer :: i,j,ic,jc, ij,ijc, kt, kdt
      Integer, allocatable :: icase(:)
      Integer, external :: DEF_ij
      Real(8) :: C
 
      Allocate(icase(ic_case))
      rewind(nud)
      Do jc = 1,ic_case  
       read(nud) ic
       read(nud) j
       read(nud) C
       read(nud) j
       read(nud) j
       read(nud) j
       read(nud) j
       icase(jc) = ic
      End do

      ij=Def_ij(ic_case,ic_case)
      if(allocated(IS_NEED)) Deallocate(IS_NEED)
      Allocate(IS_NEED(ij)); IS_NEED = 0

      Do i=1,ic_case; ic=icase(i)
       Do j=1,i; jc=icase(j)
        ijc = DEF_ij(ic,jc)
        ij = DEF_ij(i,j)
        if(JC_NEED(ijc).ne.0) IS_NEED(ij) = 1 
       End do
      End do
      Deallocate(icase)

      if(allocated(JS_NEED)) Deallocate(JS_NEED)
      Allocate(JS_NEED(ic_case)); JS_NEED=0

      Do i=1,ic_case
       Do j=1,ic_case
        ij = DEF_ij(i,j)
        if(IS_NEED(ij).eq.0) Cycle 
        JS_NEED(i)=1; Exit
       End do
      End do

      End Subroutine Define_is_done

