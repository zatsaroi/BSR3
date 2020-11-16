!======================================================================
      Subroutine get_case
!======================================================================
!     This routine obtains information about the problem to be solved
!     and reads main parameters
!----------------------------------------------------------------------
      Use bsr_mchf, nint_my => nint

      Implicit none
      Integer :: anumber = 0, nlevel1=0, nlevel2=0, i
      Character(200) :: A_core, A_conf
      Character(20) :: AA

! ... name of case:
                                                                                                  
      Call Read_name(name)

      if(LEN_TRIM(name).eq.0.or.name.eq.'?') then
       write(*,*)
       write(*,*) 'BSR_MCHF is a name-driving procedure.'
       write(*,*) 
       write(*,*) 'Each case (with spesific name) requires four file:'
       write(*,*)
       write(*,*) 'name.c - list of configurations (states)'
       write(*,*)
       write(*,*) 'name.bnk  | int.bnk - angular coefficients (after DBSR_BREIT)'
       write(*,*)
       write(*,*) 'name.knot | knor.dat - B-splines and nuclear parameters, optional'
       write(*,*)
       write(*,*) 'name.inp - running parameters, optional'
       write(*,*)
       write(*,*) 'If files are absent - I will create examples to work with.'
       write(*,*)
       write(*,*) 'Call as:   bsr_mchf  name  parameter=value  ... '
       write(*,*)
       write(*,*) 'with any paramters which are different from default values'
       write(*,*) 'or absent in the name.inp'
       write(*,*)
       Stop
      end if

! ... check the input-data file:

      z = 0.d0        !  nulify z in modules

      AF_dat = trim(name)//'.inp'
      if(Icheck_file(AF_dat).ne.0) then

       Open(inp,file=AF_dat)

       Call Read_apar(inp,'atom'   ,atom   )
       Call Read_rpar(inp,'z'      ,z      )
       Call Read_rpar(inp,'Z'      ,z      )
       Call Read_rpar(inp,'atw'    ,atw    )

       Call Read_apar(inp,'varied' ,avaried)     
       Call Read_rpar(inp,'scf_tol',scf_tol)
       Call Read_rpar(inp,'orb_tol',orb_tol)
       Call Read_rpar(inp,'end_tol',end_tol)
       Call Read_ipar(inp,'max_it' ,max_it )

       Call Read_ipar(inp,'method' ,method )
       Call Read_ipar(inp,'all'    ,all    )
       Call Read_ipar(inp,'irhs'   ,irhs   )
       Call Read_ipar(inp,'newton' ,newton )
       Call Read_ipar(inp,'rotate' ,rotate )
       Call Read_ipar(inp,'debug'  ,debug  )
       Call Read_ipar(inp,'n_corr' ,n_corr )
                                           
       Call Read_ipar(inp,'ilzero' ,ilzero )
       Call Read_ipar(inp,'ibzero' ,ibzero )

       Call Read_ipar(inp,'nlevels',nlevel1)

      end if

! ... check the command-line arguments:

      Call Read_aarg('atom'   ,atom   )
      Call Read_rarg('z'      ,z      )
      Call Read_rarg('Z'      ,z      )
      Call Read_rarg('atw'    ,atw    )

      Call Read_aarg('varied' ,avaried)     

      Call Read_rarg('scf_tol',scf_tol)
      Call Read_rarg('orb_tol',orb_tol)
      Call Read_rarg('end_tol',end_tol)
      Call Read_iarg('max_it' ,max_it )

      Call Read_iarg('method' ,method )
      Call Read_iarg('all'    ,all    )
      Call Read_iarg('irhs'   ,irhs   )
      Call Read_iarg('newton' ,newton )
      Call Read_iarg('rotate' ,rotate )
      Call Read_iarg('debug'  ,debug  )
      Call Read_iarg('n_corr' ,n_corr )
	
      Call Read_iarg('ilzero' ,ilzero )
      Call Read_iarg('ibzero' ,ibzero )

      Call Read_iarg('nlevels',nlevel2)

      Call Read_rarg('aweight',aweight)
      Call Read_rarg('bweight',bweight)
      Call Read_iarg('ac',ac)

!----------------------------------------------------------------------------------
! ... define levels for optimization:     at the moment only equal weights  ???

      if(nlevel2.gt.0) then
       nlevels=nlevel2
       Allocate(level(nlevels),weight(nlevels),elevel(nlevels), &
                ip_level(nlevels),labeln(nlevels))
       level = 0; weight = 1.d0;  weight=weight /sqrt(SUM(weight*weight))
       Call Read_iarr('level',nlevels,level)
      elseif(nlevel1.gt.0) then 
       nlevels=nlevel1
       Allocate(level(nlevels),weight(nlevels),elevel(nlevels), &
                ip_level(nlevels),labeln(nlevels))
       level = 1; weight = 1.d0
       if(nlevels.gt.1) then
        Do i=1,nlevels
         read(inp,*)  AA,level(i),AA,weight(i)
        End do
       end if
       weight=weight /sqrt(SUM(weight*weight))
      else
       nlevels=1
       Allocate(level(nlevels),weight(nlevels),elevel(nlevels), &
                ip_level(nlevels),labeln(nlevels))
       level = 1; weight = 1.d0
      end if

!------------------------------------------------------------------------------
! ... check the atom if available:

      anumber = NINT(Z)

      if(len_trim(atom).gt.0.or.anumber.gt.0) then
       Call Def_atom_LS(anumber,atom,A_core,A_conf)
       z = anumber
      else
       Stop 'Cannot define atom' 
      end if

! ... check c-file

      AF_cfg = TRIM(name)//'.c'
      Call Read_aarg('c',AF_cfg)
      if(Icheck_file(AF_cfg).eq.0) then        
       write(*,*) 'Absent c-file: ', trim(AF_cfg)
       Stop
      end if
      Open(nuc,file=AF_cfg)

! ... check bnk-file 

      AF_bnk = TRIM(name)//'.bnk'
      Call Read_aarg('bnk',AF_bnk)
      if(Icheck_file(AF_bnk).eq.0) then        
       AF_bnk = 'int_bnk'
       if(Icheck_file(AF_bnk).eq.0) then
        write(*,*) 'Absent bnk-file: ', trim(AF_bnk)
        Stop
       end if
      end if
      Open(nub,file=AF_bnk,form='UNFORMATTED')

!--------------------------------------------------------------------------------------

      AF_log = TRIM(name)//'.log'
      Open(log,file=AF_log)

      write(log,'(80("=")/T20,a/80("="))') &
         'MULTICONFIGURATION B-SPLINE HARTREE-FOCK'
 
      write(scr,'(79("=")/T15,a/79("="))') &
         'MULTICONFIGURATION B-SPLINE HARTREE-FOCK'

      write(scr,'(/a,T20,a/)') 'Name of case:',trim(name)
      write(log,'(/a,T20,a )') 'Name of case:',trim(name)

      if(Icheck_file(AF_dat).eq.0) Call Write_inp

      End Subroutine get_case


!======================================================================
      Subroutine Write_inp
!======================================================================
!     This routine prepare default input file 
!----------------------------------------------------------------------
      Use bsr_mchf

      Implicit none
      Integer :: i

      Open(inp,file=AF_dat)
      rewind(inp)

      write(inp,'(a/)') 'Main parameters:'

      write(inp,'(a,1x,a,T40,a)') 'atom    = ',trim(atom),'- atomic symbol'

      write(inp,'(a,f4.1,T40,a)') 'z       = ',z,         '- nuclear number'

      write(inp,'(/a,a)')         'varied  =  ',trim(adjustl(avaried))

      write(inp,'(/a,i2,T40,a)')  'eol     = ',eol, '- optimzation mode'

      write(inp,'(/a,i2,T40,a)')  'nlevels = ',nlevels, '- levels fo optimzation'

      if(nlevels.gt.1) then
       Do i=1,nlevels
        write(inp,'(a,i5,5x,a,f10.5)')  'level',level(i),'weight',weight(i)
       End do
      end if

      if(len_trim(physical).gt.0) &
      write(inp,'(/a,a)') 'physical=  ',trim(physical)

      Call Write_run_par(inp)

      write(inp,'(/a,a)')  'All parameters from input files ', &
                            'can be replaced from command line as:'
      write(inp,'(/10x,a/)') 'dbsr_mchf [name] par1=value par2=value par3=value ... '

      write(inp,'(80("-"))')                     
      write(inp,'(/a/)') 'Name-driven fine-name and key-words for their re-definition:' 

      write(inp,'(a,T20,a,T40,a)') 'name.inp',  ' ',       '- input parameters'
      write(inp,'(a,T20,a,T40,a)') 'name.log',  ' ',       '- run output and summry'
      write(inp,'(a,T20,a,T40,a)') 'name.c',    'c=...',   '- input configurations'
      write(inp,'(a,T20,a,T40,a)') 'name.bnk',  'bnk=...', '- angular coefficients'
      write(inp,'(a,T20,a,T40,a)') ' ',         'int_bnk', '- generic name'
      write(inp,'(a,T20,a,T40,a)') 'name.knot', 'knot=...','- B-spline parameters'
      write(inp,'(a,T20,a,T40,a)') ' ',         'knot.dat','- generic name'
      write(inp,'(a,T20,a,T40,a)') 'name.bsw',  'inp=...' ,'- input w.f. if any'
      write(inp,'(a,T20,a,T40,a)') 'name.bsw',  'out=...', '- output w.f.'
      write(inp,'(a,T20,a,T40,a)') 'name.l',    'l=...',   '- expansion coef.s in the DBSR format' 

      write(inp,'(/80("-")/)')                     
      write(inp,'(a)') ' Additional information for input parameters:'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' varied   - possible options -> all, none, list of nl, =last, n=..., n>...'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' eol      - indicates the mode for the state weights:'
      write(inp,'(a)') '            =1 - equally weighted' 
      write(inp,'(a)') '            =5 - statistically weighed, default' 
      write(inp,'(a)') '            =9 - defined by user'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' relevant level-block parameters:'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' if eol = 1 or 5, repeat for each J-block: '
      write(inp,'(a)') '            '
      write(inp,'(a)') ' levels = 1,1,2 - block, list of levels to be optimized' 
      write(inp,'(a)') '            '
      write(inp,'(a)') ' if eol = 9, repeat for each level to be optimized'
      write(inp,'(a)') '            '
      write(inp,'(a)') ' weights = 1,2,0.5 - block, level, weight' 
      write(inp,'(a)') '            '
      write(inp,'(a)') ' if level information is absent - program will optimized the first '
      write(inp,'(a)') ' level in each block'          
      write(inp,'(a)') '            '
      write(inp,'(a)') ' ilzero = 0 means  l+1 zero B-splines in the beginning for large component'

      write(inp,'(80("-"))')

      End Subroutine write_inp


!======================================================================
      Subroutine Write_run_par(nu)
!======================================================================
!     print running parameters 
!----------------------------------------------------------------------
      Use bsr_mchf
      Integer, intent(in) :: nu

      write(nu,'(/80("-")/)')
      write(nu,'(a)') 'Running parameters:'
      write(nu,*)

      write(nu,'(a,1PE9.2,T40,a)') 'scf_tol = ',scf_tol, &
                '- tolerance for energy convergence'
      write(nu,'(a,1PE9.2,T40,a)') 'orb_tol = ',orb_tol, &
                '- tolerance for orbital convergence'
      write(nu,'(a,1PE9.2,T40,a)') 'end_tol = ',end_tol, &
                '- tolerance for ending zeros'
      write(nu,'(a,i3,T40,a)') 'max_it  = ',max_it, &
                '- max. number of iterations'
      write(nu,'(a)') 
      write(nu,'(a,i2,T40,a)') 'ilzero  = ',ilzero, &
                '- initial zeros for larger component' 
      write(nu,'(a,i2,T40,a)') 'ibzero  = ',ibzero, &
                '- final zeros for larger component'

      write(nu,'(/80("-")/)')
      write(nu,'(a)') 'Additonal options (not applied if = 0)'
      write(nu,*)
 
      write(nu,'(a,i2,T40,a)') 'method  = ',method,  &
                '- method for solving MCHF equation'
      write(nu,'(a,i2,T40,a)') 'all     = ',all,    &
                '- collective optimization'  
      write(nu,'(a,i2,T40,a)') 'irhs    = ',irhs,   &
                '- convert right-hand-side to main matrix'  
      write(nu,'(a,i2,T40,a)') 'newton  = ',newton, &
                '- use Newton-Rapson method'  
      write(nu,'(a,i2,T40,a)') 'rotate  = ',rotate, &
                '- use rotations'  
      write(nu,'(a,i2,T40,a)') 'debug   = ',debug,  &
                '- additional debug output'

      End Subroutine Write_run_par
