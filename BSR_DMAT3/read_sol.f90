!======================================================================
      Subroutine Read_sol(ctype,nu,nc,C,Label,E,jot)
!======================================================================
!     read one solution from c-,l-,j- or b-files 
!
!     nu    -  file unit
!     ctype -  type of file (c,l,j,b)
!     nc    -  dimention of solution
!     C     -  solution
!     Label -  spectroscopic notation
!     E     -  energy 
!     jot   -  (2J+1) for j-type calculations
!
!     the first line in files 'l,j,b' is supposed to have already 
!     been read !
!----------------------------------------------------------------------
      Implicit none

      Character(1), intent(in) :: ctype
      Character(64), intent(out) :: Label
      Character(80) :: AS
      Integer, intent(in) :: nu,nc 
      Integer, intent(inout) :: jot 
      Integer :: i
      Real(8), intent(out) :: E
      Real(8), intent(out) :: C(nc)
 
      Select Case(ctype)

      Case('b')

       read(nu,'(7x,a64)') Label
       read(nu,'(D20.10)') E  
       read(nu,'(5D15.8)') C

      Case('c')
	  
       Call R_expn(nu,nc,C)
       Call Idef_LabelC(nu,0,0,Label)        
       rewind(nu); read(nu,'(a)') AS; read(AS,'(15x,f16.8)') E
       jot=0
       i=INDEX(AS,'=')
       if(i.gt.0) then; read(AS(i+1:),*) jot; jot=jot+1; end if

      Case('j') 

      1 read(nu,'(a80)',end=2) AS
       if(AS(14:14).eq.'.') then
         read(AS,'(6x,F16.8,2x,a64)') E, Label 
         read(nu,'(7F11.8)') C
        elseif(AS(5:5).eq.'J') then
         read(AS(8:),*) jot; jot=jot+1; go to 1
        else
         go to 1
        end if
        Return
      2 write(*,*) 'Read_sol: read the file end '
        Stop

       Case('l') 

      3 read(nu,'(a80)',end=4) AS
        if(AS(14:14).eq.'.') then
         read(AS,'(6x,F16.8,2x,a64)') E, Label
         read(nu,'(7F11.8)') C
        else
         go to 3
        end if
        Return
      4 write(*,*) 'Read_sol: read the end of file'
        Stop

      Case Default
	    
       Stop ' Read_sol: unknown case '

      End Select
 
      End Subroutine Read_sol


