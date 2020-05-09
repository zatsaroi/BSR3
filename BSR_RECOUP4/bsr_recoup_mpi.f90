!======================================================================
!     PROGRAM       B S R _ R E C O U P
!
!                   C O P Y R I G H T -- 2017
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!                        
!======================================================================
!     Make recoupling the LS Hamiltionian/Overlap matrixes 
!     to the JK coupling based on the target expansions over LS-states
!======================================================================
!
!     INPUT FILES:
!
!     bsr_par        -  parameters of calculations
!     target         -  description of target states and partial waves                        
!                       in JK-coupling
!     target_LS      -  description of target states and partial waves
!                          in LS-coupling                                
!     bsr_mat_LS.nnn -  LS Hamiltonian matrix for nnn-partion wawe 
!                    
!     target_exp     -  LSJ target expansions over LS-states
!
!     OUTPUT FILES:  
!                    
!     bsr_recoup.log -  running information
!
!     bsr_mat_JK.nnn -  LS Hamiltonian matrix for nnn-partion wawe
!                    
!     INPUT ARGUMENT:  see subroutine read_arg.f90
!
!----------------------------------------------------------------------
      Use MPI
      Use bsr_recoup

      Implicit real(8) (A-H,O-Z)
      Character(80) :: AS

      Call MPI_INIT(ierr)
      Call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

! ... read arguments:
if(myid.eq.0) write(*,*) 'nprocs',nprocs
      if(myid.eq.0)  Call Read_arg
      Call br_arg
if(myid.eq.0) write(*,*) 'klsp1,klsp2',klsp1,klsp2
!---------------------------------------------------------------------- 
! ... read target information:

      if(myid.eq.0) then

       Call Check_file(AF_tar)
       Open(nut,file=AF_tar,status='OLD',action='READ')
       Call R_target(nut)
if(myid.eq.0) write(*,*) AF_tar
       Call Check_file(BF_tar)
       Open(mut,file=BF_tar,status='OLD',action='READ')
       Call R_target_ion(mut)
       Call R_channels_ion(mut)
if(myid.eq.0) write(*,*) BF_tar

      end if

      Call br_target
      Call br_target_ion
      Call br_channels_ion

! ... read recoupling coefficients:                                                                        
                                                                                                            
      if(myid.eq.0) then                                                                                   
       Call Check_file(AF_expn)                                                                              
       Open(nue,file=AF_expn,status='OLD',action='READ')                                                                                
       Call Read_ipar(nue,'mpert',mpert)                                                                     
       Allocate(ip_expn(0:ntarg),it_expn(mpert),c_expn(mpert))                                               
       ip_expn = 0                                                                                           
       rewind(nue)                                                                                           
       k = 0                                                                                                 
       Do it = 1,ntarg                                                                                       
        read(nue,*)  AS,jt,AS,AS,npert                                                                       
        if(it.ne.jt) Stop 'it <> jt  in target_LS_expn'                                                      
        Do i = 1,npert                                                                                       
         k = k + 1                                                                                           
         read(nue,*) it_expn(k),c_expn(k)                                                                    
        End do                                                                                               
        ip_expn(it) = k      
       End do
       if(k.ne.mpert) Stop 'mpert ??? on target_LS_expn'
      end if

if(myid.eq.0) write(*,*) AF_expn

! ... broadcast recoupling coefficients: 

      Call MPI_BCAST(mpert, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(myid.ne.0)  Allocate(ip_expn(0:ntarg),it_expn(mpert),c_expn(mpert))                                               
      Call MPI_BCAST(ip_expn(0:ntarg),ntarg+1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(it_expn,mpert,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      Call MPI_BCAST(c_expn,mpert,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!----------------------------------------------------------------------
! ... loop over partial waves:

      Do klsp = klsp1,klsp2
 
       i=len_trim(AF_p)
       write(AF_p(i-2:i),'(i3.3)') klsp
       open(pri,file=AF_p,action='WRITE')
       if(myid.le.debug) then
        write(AF,'(a,a,i4.4)') trim(AF_p),'.',myid
        if(myid.eq.0) AF = AF_p
        open(pri,file=AF)
       else
        pri=0 
       end if

       if(pri.gt.0) write(pri,'(a,i5/)') 'MPI: nprocs = ',nprocs

       if(pri.gt.0) &
       write(pri,'(a,i3/)') 'BSR_RECOUP:  partial wave klsp =', klsp

       Call CPU_TIME(t1);  Call SUB1;  Call CPU_TIME(t2)

       if(pri.gt.0) &
       write(pri,'(/a,5x,T20,f8.2,a/)') 'Total time:', (t2-t1)/60, ' min'

      End do  ! over klsp

      Call MPI_FINALIZE(ierr)

      End  ! program bsr_mat
