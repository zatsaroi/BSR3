!======================================================================
      Subroutine Read_arg(nu)
!======================================================================
!     read arguments from file unit 'nu'
!----------------------------------------------------------------------
      Use bsr_phot
    
      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,i1,i2
      Character(80) :: AS

      rewind(nu); Elow=0.d0; Ehigh=0.d0; Estep=0.d0 
    1 read(nu,'(a)',end=2) AS
      i=INDEX(AS,'='); if(i.lt.2) go to 1
      i1=i+1; i2=INDEX(AS,'!')-1; if(i2.lt.i1) i2=LEN_TRIM(AS)
      Select Case(AS(1:i-1))
       Case('klsp'  );  read(AS(i1:i2),*) klsp
       Case('iauto' );  read(AS(i1:i2),*) iauto
       Case('mfgi ' );  read(AS(i1:i2),*) mfgi
       Case('AC'    );  read(AS(i1:i2),*) AC
       Case('awt'   );  read(AS(i1:i2),*) AWT
       Case('e_exp' );  read(AS(i1:i2),*) e_exp
       Case('DR'    );  read(AS(i1:i2),*) DR
       Case('ibug'  );  read(AS(i1:i2),*) ibug
       Case('elow'  );  read(AS(i1:i2),*,end=1) ELOW
       Case('estep' );  read(AS(i1:i2),*,end=1) ESTEP
       Case('ehigh' );  read(AS(i1:i2),*,end=1) EHIGH
       Case('nwt'   );  read(AS(i1:i2),*) nwt
       Case('ikm'   );  read(AS(i1:i2),*) ikm
       Case('ath'   );  read(AS(i1:i2),*) athreshold
       Case('bth'   );  read(AS(i1:i2),*) bthreshold
      End Select
      go to 1
    2 Continue   

      Call Read_iarg('klsp',klsp )

      End Subroutine Read_arg

