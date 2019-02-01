# ifdef VISIT
!!!====================================================================
!!!====================================================================
!!!====================================================================
!!!
!!!         BEGIN VISIT SUBROUTINES
!!!
!!!====================================================================
!!!====================================================================
!!!====================================================================

!NOTES:
! Variables not declared at the header of each subroutine come from
! the other modules in FVCOM.


SUBROUTINE INIT_VISIT
!  USE MOD_UTILS
  USE CONTROL, only : fvcom_version, casetitle, PAR
  USE LIMS, only : MYID
  USE MOD_VISIT
  USE MOD_UTILS
  implicit none
  include "visitfortransiminterface.inc"

  TYPE(VisitMeshType),pointer  :: VMS
! local variables
  integer err,ind
  character(len=100) :: SimName, timestamp
  integer :: LSimName, LTitle

# ifdef VISIT_DEBUG

  character(LEN=3) :: ch3
  character(len=100) :: debugname

  
  if(.not. MSR) then


     write(ch3,'(i3)') myid
     if(myid<10 .and. myid>=1) then
     debugname="FVCOM_VISIT_DEBUG." // trim(adjustl(ch3)) // ".log"
     else if(myid<100 .and. myid>=10) then
     debugname="FVCOM_VISIT_DEBUG." // trim(adjustl(ch3)) // ".log"
     else if(myid<1000 .and. myid>=100) then
     debugname="FVCOM_VISIT_DEBUG." // trim(adjustl(ch3)) // ".log"
     end if

     VDB=164 + myid
     CALL FOPEN(VDB, TRIM(debugname) ,"ofr")
     
  else
     VDB=IPT
     
  end if


  WRITE(VDB,*) "BEGIN INIT_VISIT: PID=",VisitParRank

  SLAVECALLBACKCOUNT=0
  BROADCASTSTRCOUNT=0
  BROADCASTINTCOUNT=0

  if (MSR) then
     WRITE(VDB,*)"=================================="
     WRITE(VDB,*)"       VISIT PARAMETERS"
     WRITE(VDB,*)"=================================="
     WRITE(VDB,*)"VISIT_TWODMESH=",VISIT_TWODMESH
     WRITE(VDB,*)"VISIT_SSHMESH=",VISIT_SSHMESH
     WRITE(VDB,*)"VISIT_BATHYMESH=",VISIT_BATHYMESH
     WRITE(VDB,*)"VISIT_LAYERMESH=",VISIT_LAYERMESH
     WRITE(VDB,*)"VISIT_LEVELMESH=",VISIT_LEVELMESH

     WRITE(VDB,*)"VISIT_COMMAND_PROCESS=",VISIT_COMMAND_PROCESS
     WRITE(VDB,*)"VISIT_COMMAND_SUCCESS=",VISIT_COMMAND_SUCCESS
     WRITE(VDB,*)"VISIT_COMMAND_FAILURE=",VISIT_COMMAND_FAILURE

  end if

# else

  if(MSR) VDB = IPT

# endif

  ! SET MESH VARIABLE TYPES
  VISIT_MESH(VISIT_TWODMESH)%MESH=  VISIT_TWODMESH

  VISIT_MESH(VISIT_SSHMESH)%MESH=  VISIT_SSHMESH

  VISIT_MESH(VISIT_BATHYMESH)%MESH=  VISIT_BATHYMESH

  VISIT_MESH(VISIT_LAYERMESH)%MESH=  VISIT_LAYERMESH

  VISIT_MESH(VISIT_LEVELMESH)%MESH=  VISIT_LEVELMESH

# ifdef VISIT_DEBUG
  if (MSR) then
     VMS =>  VISIT_MESH(VISIT_TWODMESH)
     if(associated(VMS))  WRITE(VDB,*)"associated"
     Call PRINT_MESH(VMS)
  end if
# endif


! INITIALIZE SOME VISIT VARIABLES
  VisitRunFlag=1 ! Set to zero to wait for viewer at the first time step
  !MYID is result of MPI_comm_rank PLUS one!
  VisitParRank=MYID-1

  VisitOneTStep=.FALSE.
  VisitHalt=.false.
  VisitAnimate=.FALSE.
  VISIT_CMD_DUMP=.FALSE.

  VisitStep=VISIT_STEP_INT


  err = visitsetupenv()

#  if defined (MULTIPROCESSOR)
  err = visitsetparallel(PAR)
  err= visitsetparallelrank(VisitParRank)

# ifdef VISIT_DEBUG
  WRITE(VDB,*)"visitsetparallel(PAR)=",err,": PAR=",PAR,": PID=",VisitParRank
  WRITE(VDB,*)"visitsetparallelrank=",err,": PID=",VisitParRank
# endif
  
# endif


  if(MSR) then
!     call get_timestamp(timestamp)
     LTitle=Len_Trim(casetitle)
     SimName= trim(fvcom_version) ! could include fvcom timestamp BUT
     ! NO COLONS!
     LSimName=Len_Trim(SimName)
!     WRITE(VDB,*)trim(SimName),",",LSimName
!     WRITE(VDB,*)trim(casetitle),',',Ltitle

     err = visitinitializesim(trim(simname),lsimname, trim(casetitle), ltitle,&
          & "/no/useful/path", 15, visit_f77nullstring,&
          & visit_f77nullstringlen, visit_f77nullstring, visit_f77nullstringlen)

!     WRITE(VDB,*) "MASTER WROTE VISIT SIM FILE: ",trim(simname)
  end if


# ifdef VISIT_DEBUG
  WRITE(VDB,*) "END INIT_VISIT: PID",VisitParRank
# endif
  
end SUBROUTINE INIT_VISIT

!!====================================================================
!!====================================================================
!!
!! VISITCHECK decides when check the viewer for connections and commands
!!
!!====================================================================
!!====================================================================


SUBROUTINE VisitCheck
  USE MOD_VISIT
  USE CONTROL, only : IINT, IEXT
  implicit none
  
  
  if (VisitStep .EQ. VISIT_STEP_EXT) then
     Call VisitSimLoop
     
  elseif (VisitStep .EQ. VISIT_STEP_INT) then
     if (IEXT .EQ. 1) Call VisitSimLoop
     
  elseif (VisitStep .EQ. VISIT_STEP_10XINT) then
     if(IEXT .EQ. 1) then
        VisitStepCount=VisitStepCount+1
        if (mod(VisitStepCount,10) .EQ. 0) Call VisitSimLoop
        
     end if
     
  elseif (VisitStep .EQ. VISIT_STEP_100XINT) then
     if(IEXT .EQ. 1) then
        VisitStepCount=VisitStepCount+1
        if (mod(VisitStepCount,100) .EQ. 0) Call VisitSimLoop
        
     end if
     
  end if
  
END SUBROUTINE VisitCheck



!!====================================================================
!!====================================================================
!!
!! VisitSimLoop
!!
!! This is the loop that does all the work of talking to visit
!!====================================================================
!!====================================================================


subroutine VisitSimLoop
  USE MOD_VISIT
# if defined (MULTIPROCESSOR)
  USE MPI
# endif
  implicit none
  include "visitfortransiminterface.inc"
  integer :: visitstate, result, blocking, loopcount,IERR, ind
  integer, external :: processvisitcommand
  integer, external :: visitbroadcastintfunction

  TYPE(VISITMESHTYPE), POINTER :: vmp ! VISIT MESH POINTER

  TYPE(VISITLAG), POINTER :: lmp ! VISIT LAGRANGIAN POINTER

  TYPE(VisitSphereVel), pointer :: VSV


! VisitRunFlag == 1       Simulation is Running

! VisitRunFlag == 0       Simulation is Stopped


# ifdef VISIT_DEBUG
  WRITE(VDB,*)"BEGIN DoVisitSim: PID", VisitParRank
# endif


  IF (VISITANIMATE) then
     IF (MSR) WRITE(VDB,*)"Updating Plots:"
     result = visittimestepchanged()
     result = visitupdateplots()
  END IF
  
  if (VisitOneTStep) then
# ifdef VISIT_DEBUG
  WRITE(VDB,*)"Hit VisitOneTStep at beginning of visitsimloop"
# endif
     VisitOneTStep=.false.
     VisitHalt=.true.
  end if

  DO loopcount=1,10000 ! make sure loop is not infinite!

     ! LOGIC TO DO ONE TIME STEP ONLY
     
# ifdef VISIT_DEBUG
     WRITE(VDB,*)"Start VISIT SIM LOOP COUNT=", loopcount,": PID=",VisitParRank
     WRITE(VDB,*)"visitRunFlag=",visitrunflag
# endif
     
     if(VisitRunFlag.eq.1) then
        blocking = 0 ! Simulation running
     else !
        blocking = 1 ! Simulation stopped, waiting for commands
     endif
 
!     !Allow the one step call back to over ride visitrunflag
!     if (VisitOneTStep .eq. .TRUE.) then
!        blocking =0
!        VisitOneTStep=.False.
!        VisitRunFlag=0  ! Make sure FVCOM only goes one step!
!     end if

# ifdef VISIT_DEBUG
     WRITE(VDB,*)"blocking=",blocking
# endif

     if(MSR) visitstate = visitdetectinput(blocking, -1)
     
# ifdef VISIT_DEBUG
     IF(MSR) WRITE(VDB,*)"Sending visitstate from master, VISITSTATE=",VISITSTATE
     IF(.not.MSR) WRITE(VDB,*)"Recieving visitstate from master"
# endif

     IERR= visitbroadcastintfunction(visitstate, 0)

# ifdef VISIT_DEBUG
     WRITE(VDB,*) "Sent: VISITSTATE= ",VISITSTATE,": PID=",VisitParRank
# endif

     if (visitstate.lt.0) then
        
        if(MSR)then ! equivlent to VisitParRank=0
           WRITE(VDB,*)'VISITSTATE error codes: ',VISITSTATE
           WRITE(VDB,*)"-5: Logic error (fell through all cases)"
           WRITE(VDB,*)'-4: Logic error (no descriptors but blocking)'
           WRITE(VDB,*)'-3: Logic error (a socket was selected but not one we set)'
           WRITE(VDB,*)"-2: Unknown error in select"
           WRITE(VDB,*)"-1: Interrupted by EINTR in select"
        end if
        
        ! Stop the simulation!
        CALL PSTOP
        
     elseif (visitstate.eq.0) then

# ifdef VISIT_DEBUG
        WRITE(VDB,*) "VisitRunFlag= ", VisitRunFlag,": PID=",VisitParRank
# endif

        if (VisitHalt) then
# ifdef VISIT_DEBUG
           WRITE(VDB,*)"Hit VisitHalt while exiting loop with visitsta&
                &te == 0 "
# endif
           VisitRunFlag=0
           VisitHalt=.false.
        else
           exit ! return to US_FVCOM Main loop
        end if

        !!call simulate_one_timestep()
     elseif (visitstate.eq.1) then
        VisitRunFlag = 0
        result = visitattemptconnection()

        if (MSR) then
           if (result.eq.1) then
              WRITE(VDB,*) 'VisIt connected!'
           else
              WRITE(VDB,*)'VisIt did not connect!'
           endif
        end if


        if (result.eq.1) VisitAnimate=.TRUE.


     elseif (visitstate.eq.2) then

# ifdef VISIT_DEBUG
        if(MSR) WRITE(VDB,*)'Run ProcessVisitCommand'
# endif
        IERR = ProcessVisitCommand()
        
        if(IERR.eq.0) then
           if(MSR) WRITE(VDB,*)'ProcessVisitCommand returned disconnect PID=', visitParRank
           result = visitdisconnect()
           VisitRunFlag = 1
           VisitAnimate=.FALSE.
           VisitOneTStep=.false.
           VisitHalt=.false.
           VISIT_CMD_DUMP=.FALSE.

!           VisitStep= VISIT_STEP_10XINT
           
           VMP => VISIT_MESH(VISIT_TWODMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_SSHMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_BATHYMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_LAYERMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           VMP => VISIT_MESH(VISIT_LEVELMESH)
           Call DeAllocate_Mesh(VMP,IERR)

           NULLIFY(VMP)

#  if defined (NETCDF_IO)
           LMP => VISIT_LAGDATA
           Call DeAllocate_LAG(LMP,IERR)

           NULLIFY(LMP)
# endif

# if defined (SPHERICAL)
           VSV =>VisitSphericalVel
           Call DeAllocate_SphereVel(VSV,IERR)

           VSV =>VisitSphericalAVel
           Call DeAllocate_SphereVel(VSV,IERR)

           VSV =>VisitSphericalWindVel
           Call DeAllocate_SphereVel(VSV,IERR)

# if defined (ICE)
           VSV =>VisitSphericalIceVel
           Call DeAllocate_SphereVel(VSV,IERR)
# endif           

           NULLIFY(VSV)
# endif


           exit ! Return to simulation loop
        endif

     endif
     
# ifdef VISIT_DEBUG
     WRITE(VDB,*)"END VISIT SIM LOOP COUNT=", loopcount,": PID=",VisitParRank
# endif

  end DO
  

# ifdef VISIT_DEBUG
  WRITE(VDB,*)"END DoVisitSim: PID=",VisitParRank
# endif

end subroutine VisitSimLoop

!!=============================================================
!!=============================================================
!!=============================================================

! PARRALLEL HELPER COMMAND FOR DOVISITSIM
integer function ProcessVisitCommand()
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  integer :: IERR, processcount
  integer :: command, success
  
# ifdef VISIT_DEBUG
  WRITE(VDB,*)"BEGIN ProcessVisitCommand: PID=",VisitParRank
# endif

  IERR=0

  if(MSR) then

# ifdef VISIT_DEBUG
     WRITE(VDB,*)"Master Processing Engine command:"
# endif
     success =  visitprocessenginecommand()
     
     if(success.eq.1) then

# ifdef VISIT_DEBUG
        WRITE(VDB,*)"Master Success Processing Engine command: BROADCAST2SLAVES"
#endif
        command = VISIT_COMMAND_SUCCESS
        IERR=1
     else
# ifdef VISIT_DEBUG
        WRITE(VDB,*)"Master Failure Processing Engine command: BROADCAST2SLAVES"
# endif
        command = VISIT_COMMAND_FAILURE
        IERR=0
     end if
     
     call BroadcastSlaveCommand(command)
     
  else
     
     processcount=0
     DO

# ifdef VISIT_DEBUG
        WRITE(VDB,*)"Slave asking for process engine command: processcount=&
             &",processcount
# endif
        processcount=processcount+1
        call BroadcastSlaveCommand(command)
        SELECT CASE (command)
           
        CASE (VISIT_COMMAND_PROCESS)

# ifdef VISIT_DEBUG
           WRITE(VDB,*)"Slave: process engine command"
# endif

           success=visitprocessenginecommand()
           IERR=1

        CASE (VISIT_COMMAND_SUCCESS)
# ifdef VISIT_DEBUG
           WRITE(VDB,*)"Slave recieved: visit command success"
# endif
           IERR =1
           exit
        CASE (VISIT_COMMAND_FAILURE)
# ifdef VISIT_DEBUG
           WRITE(VDB,*)"Slave recieved: visit command failure"
# endif
           IERR=0
           exit
           
        CASE DEFAULT
# ifdef VISIT_DEBUG
           WRITE(VDB,*)"Slave recieved: !Unrecognized visit command!"
           WRITE(VDB,*)"RANK: ",VisitParRank,"; Bad Command, Calling PSTOP:"
# endif
           call PSTOP
        END SELECT
     END DO

  end if

# ifdef VISIT_DEBUG
  WRITE(VDB,*)"END ProcessVisitCommand: PID=",VisitParRank
# endif

  ProcessVisitCommand = IERR
end function ProcessVisitCommand
!! END PARALLEL HELPER FOR DoVisitSim


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c These functions must be defined to satisfy the visitfortransiminterface lib.
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c---------------------------------------------------------------------------
!c visitcommandcallback
!c---------------------------------------------------------------------------
subroutine visitcommandcallback (cmd, lcmd, intdata, floatdata,&
     & stringdata, lstringdata)
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(in) :: cmd, stringdata
  integer, intent(in) ::lcmd, lstringdata, intdata
  real, intent(in) :: floatdata
  integer :: result

# ifdef VISIT_DEBUG
  WRITE(VDB,*) "BEGIN VISITCOMMANDCALLBACK: PID=",VisitParRank
# endif
  
  if(visitstrcmp(cmd, lcmd, "Halt/Step", 9).eq.0) then
     if (VisitRunFlag .EQ. 1) then
        ! VisitRunFlag = 1      ! if running, halt
        VisitHalt=.TRUE.     ! catch at next entry to visitsimloop

     else  ! visitRunFlag .eq. 0  - stopped
        visitRunFlag =1       ! if stopped go one time step
        VisitOneTStep=.TRUE.
     end if

  elseif(visitstrcmp(cmd, lcmd, "Run", 3).eq.0) then
     VisitRunFlag = 1

  elseif(visitstrcmp(cmd, lcmd, "Dump Restart", 12).eq.0) then
     VISIT_CMD_DUMP = .true.

  elseif(visitstrcmp(cmd, lcmd, "Toggle Step", 11).eq.0) then


     SELECT CASE(VisitStep)

     CASE(VISIT_STEP_100XINT)
        VisitStep=VISIT_STEP_EXT
        if(MSR) then
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"     VISIT TIME STEP CHANGED TO EVERY EXT"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
        end if
     CASE(VISIT_STEP_EXT)
        VisitStep=VISIT_STEP_INT
        if (MSR) then
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"     VISIT TIME STEP CHANGED TO ONE IINT"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
        end if
     CASE(VISIT_STEP_INT)
        VisitStep=VISIT_STEP_10XINT
        if(MSR) then
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"     VISIT TIME STEP CHANGED TO 10 X IINT"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
        end if
     CASE(VISIT_STEP_10XINT)
        VisitStep=VISIT_STEP_100XINT
        if(MSR) then
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"     VISIT TIME STEP CHANGED TO 100 X IINT"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
        end if
     CASE DEFAULT
        VisitStep=VISIT_STEP_INT
        if (MSR) then
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)" VISIT TIME STEP CHANGED TO ONE IINT: DEFAULT"
           WRITE(VDB,*)"================================================"
           WRITE(VDB,*)"================================================"
        end if
     END SELECT


!  elseif(visitstrcmp(cmd, lcmd, "ext", 3).eq.0) then
!     VisitStep=VISIT_STEP_EXT
!
!  elseif(visitstrcmp(cmd, lcmd, "int", 3).eq.0) then
!     VisitStep=VISIT_STEP_INT
!
!  elseif(visitstrcmp(cmd, lcmd, "10Xint", 6).eq.0) then
!     VisitStep=VISIT_STEP_10XINT
!     VisitStepCount=0
!
!  elseif(visitstrcmp(cmd, lcmd, "100Xint", 7).eq.0) then
!     VisitStep=VISIT_STEP_100XINT
!     VisitStepCount=0
!

  endif


  
# ifdef VISIT_DEBUG
  WRITE(VDB,*) "END VISITCOMMANDCALLBACK: PID=",VisitParRank
# endif

end subroutine visitcommandcallback

!c---------------------------------------------------------------------------
!c visitbroadcastintfunction
!c---------------------------------------------------------------------------
integer function visitbroadcastintfunction(value, sender)
# if defined (MULTIPROCESSOR)
  USE MPI
  USE MOD_VISIT
# endif
  implicit none
  include "visitfortransiminterface.inc"
  integer, intent(inout) :: value, sender
  integer ::IERR

#  if defined (MULTIPROCESSOR)

# ifdef VISIT_DEBUG
  BROADCASTINTCOUNT=BROADCASTINTCOUNT+1
  WRITE(VDB,*)"BEGIN VISITBROADCASTINTFUNCTION: PID=",VisitParRank,": COUNT&
       &=", BROADCASTINTCOUNT
# endif

  CALL MPI_BCAST(value,1,MPI_INTEGER,sender,MPI_COMM_WORLD,IERR)

# ifdef VISIT_DEBUG
  WRITE(VDB,*)"END VISITBROADCASTINTFUNCTION: PID=",VisitParRank,": COUNT&
       &=", BROADCASTINTCOUNT
# endif

  ! endif MPI
# endif


  visitbroadcastintfunction = 0


end function visitbroadcastintfunction

!c---------------------------------------------------------------------------
!c visitbroadcaststringfunction
!c---------------------------------------------------------------------------
integer function visitbroadcaststringfunction(str, lstr, sender)
# if defined (MULTIPROCESSOR)
  USE MPI
  USE MOD_VISIT
# endif
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(inout):: str
  integer, intent(in):: sender, lstr
  integer :: IERR


# if defined (MULTIPROCESSOR)

# ifdef VISIT_DEBUG
  BROADCASTSTRCOUNT= BROADCASTSTRCOUNT +1
  WRITE(VDB,*)"BEGIN VISITBROADCASTSTRINGFUNCTION: PID=",VisitParRank,": COUNT&
       &=", BROADCASTSTRCOUNT
#endif

  CALL MPI_BCAST(str,lstr,MPI_CHARACTER, sender ,MPI_COMM_WORLD,IERR)

# ifdef VISIT_DEBUG
  WRITE(VDB,*)"END VISITBROADCASTSTRINGFUNCTION: PID=",VisitParRank,": COUNT&
       &=", BROADCASTSTRCOUNT
# endif

# endif

  visitbroadcaststringfunction = 0

end function visitbroadcaststringfunction

!c---------------------------------------------------------------------------
!c visitslaveprocesscallback
!c---------------------------------------------------------------------------
subroutine visitslaveprocesscallback()
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  integer :: command

  command = VISIT_COMMAND_PROCESS
 
# ifdef VISIT_DEBUG
  WRITE(VDB,*) "Begin visitslaveprocesscallback: PID=",VisitParRank
  WRITE(VDB,*)"command=",command
# endif

  Call BroadcastSlaveCommand(command)

# ifdef VISIT_DEBUG
  WRITE(VDB,*) "End visitslaveprocesscallback: PID=",VisitParRank
# endif

end subroutine visitslaveprocesscallback


! HELPER FUNCTION FOR BROADCASTING COMMANDS FROM MASTER
subroutine BroadcastSlaveCommand(command)
# if defined (MULTIPROCESSOR)
  USE MPI
  USE MOD_VISIT
# endif
  implicit none
  integer, intent(in) :: command
  integer :: IERR

#  if defined (MULTIPROCESSOR)


# ifdef VISIT_DEBUG
  SLAVECALLBACKCOUNT=SLAVECALLBACKCOUNT+1
  WRITE(VDB,*) "Broadcast from master to slaves:",command,": PID=",VisitParRank,": COUNT&
       &=", SLAVECALLBACKCOUNT
# endif

  CALL MPI_BCAST(command,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)


# ifdef VISIT_DEBUG
  WRITE(VDB,*) "Broadcast from master to slaves complete: PID=",VisitParRank,": COUNT&
       &=", SLAVECALLBACKCOUNT
# endif

# endif

end subroutine BroadcastSlaveCommand



!c---------------------------------------------------------------------------
!c visitgetmetadata
!c---------------------------------------------------------------------------
integer function visitgetmetadata(handle)
  USE ALL_VARS, only : NPROCS, KBM1, KBM2, KB, PAR
  USE MOD_VISIT

  implicit none
  include "visitfortransiminterface.inc"
  integer, intent(in):: handle
  integer :: err, ind
  ! Mesh Variables
  integer :: tdim2, tdim3, sdim2,sdim3, mesh1,mesh2,mesh3,mesh4,mesh5&
       &, meshtype, meshlag
  character(LEN=20) :: UnitsName, LengthName, WidthName, DepthName,&
       & DomainName
  integer ::  LUnitsName, LLengthName, LWidthName, LDepthName, LDomainName 
  character(LEN=3) :: ch3

  character(LEN=20) :: matname
  character(LEN=20), allocatable :: matnamecopy(:)

  ! Scalar variables
  integer :: scalar

  ! Materials
  integer :: mat1, mat2
  ! Time stuff


  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! TIME AND CYCLE
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ifdef VISIT_DEBUG
  WRITE(VDB,*) "VisitSimCycle= ",Visit_Cycle,": PID=",VisitParRank
  WRITE(VDB,*) "VisitSimTime= ", Visit_Time_Ext,": PID=",VisitParRank
# endif

  err= visitmdsetcycletime(handle, Visit_Cycle, Visit_Time_Ext)  

# ifdef VISIT_DEBUG
  WRITE(VDB,*)"visitmdsetcycletime:result=",err
# endif

  if(VisitRunFlag.eq.0)then
     err = visitmdsetrunning(handle, VISIT_SIMMODE_RUNNING)
!     WRITE(VDB,*)"visitmdsetrunning(running):result=",err
  else
     err = visitmdsetrunning(handle, VISIT_SIMMODE_STOPPED)
!     WRITE(VDB,*)"visitmdsetrunning(stopped):result=",err
  end if


  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! MESH META DATA
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! General mesh variables
  meshtype = VISIT_MESHTYPE_UNSTRUCTURED
  tdim2 = 2
  tdim3 = 3

  sdim2 = 2
  sdim3 = 3

# if defined (SPHERICAL)


  UnitsName="Meters"
  LUnitsName=Len_Trim(UnitsName)

  LengthName="X"
  LLengthName=Len_Trim(LengthName)

  WidthName="Y"
  LWidthName=Len_Trim(WidthName)

  DepthName="Z"
  LDepthName=Len_Trim(DepthName)


# else
  
  UnitsName="Meters"
  LUnitsName=Len_Trim(UnitsName)
  
  LengthName="Zonal"
  LLengthName=Len_Trim(LengthName)
  
  WidthName="Meridional"
  LWidthName=Len_Trim(WidthName)
  
  DepthName="Depth"
  LDepthName=Len_Trim(DepthName)
  
# endif


  DomainName="Domain_"
  LDomainName=Len_Trim(DomainName)


# if defined (SPHERICAL)
  ! TWOD_Mesh
  mesh1 = visitmdmeshcreate(handle, "TWOD_Mesh", 9, meshtype, tdim2 &
       &,sdim3, NPROCS)
# else
  ! TWOD_Mesh
  mesh1 = visitmdmeshcreate(handle, "TWOD_Mesh", 9, meshtype, tdim2 &
       &,sdim2, NPROCS)
# endif

  if(mesh1.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh1, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh1, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName) ! Depth not used
     err = visitmdmeshsetblocktitle(handle, mesh1, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh1,DomainName , LDomainName)
  endif

  ! BATHYMETRY_Mesh
  mesh2 = visitmdmeshcreate(handle, "Bathymetry_Mesh", 15, meshtype, tdim2 &
       &,sdim3, NPROCS)

  if(mesh2.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh2, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh2, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh2, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh2,DomainName , LDomainName)
  endif

  ! SSH_Mesh
  mesh3 = visitmdmeshcreate(handle, "SSH_Mesh", 8, meshtype, tdim2 &
       &,sdim3, NPROCS)

  if(mesh3.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh3, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh3, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh3, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh3,DomainName , LDomainName)
  endif

  ! SIGMALAYER_Mesh
  mesh4 = visitmdmeshcreate(handle, "SigmaLayer_Mesh", 15, meshtype, tdim3 &
       &,sdim3, NPROCS)

  if(mesh4.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh4, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh4, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh4, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh4,DomainName , LDomainName)
  endif

  ! SIGMALEVEL_Mesh
  mesh5 = visitmdmeshcreate(handle, "SigmaLevel_Mesh", 15, meshtype, tdim3 &
       &,sdim3, NPROCS)

  if(mesh5.ne.VISIT_INVALID_HANDLE) then
     err = visitmdmeshsetunits(handle, mesh5, UnitsName, LUnitsName)
     err = visitmdmeshsetlabels(handle, mesh5, LengthName,LLengthName&
          &, WidthName, LWidthName, DepthName,LDepthName)
     err = visitmdmeshsetblocktitle(handle, mesh5, "Domains", 7)
     err = visitmdmeshsetblockpiecename(handle, mesh5,DomainName , LDomainName)
  endif


  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! ADD ONLINE LAGRANGIAN PARTICLE METADATA
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#  if defined (NETCDF_IO)
   if (LAG_ON) then

# if defined (SHPERICAL)
      if(MSR) WRITE(VDB,*)"VISIT CAN NOT USE LAGRANGIAN DATA IN SPERICAL CO&
           &ORDINATES"
      CALL PSTOP
# endif

      meshlag = visitmdmeshcreate(handle, "LAGRANGIAN_MESH", 15,&
           & VISIT_MESHTYPE_POINT,  0, 3, NPROCS)
      if(meshlag.ne.VISIT_INVALID_HANDLE) then
         err = visitmdmeshsetunits(handle, meshlag, UnitsName, LUnitsName)
         err = visitmdmeshsetlabels(handle, meshlag, LengthName,LLengthName&
              &, WidthName, LWidthName, DepthName,LDepthName)
         err = visitmdmeshsetblocktitle(handle, meshlag, "Domains", 7)
         err = visitmdmeshsetblockpiecename(handle, meshlag,DomainName , LDomainName)
      endif

      ! EASTWARD VELOCITY
      scalar = visitmdscalarcreate(handle, "LAG_U", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! NORTHWARD VELOCITY
      scalar = visitmdscalarcreate(handle, "LAG_V", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! UPWARD VELOCITY
      scalar = visitmdscalarcreate(handle, "LAG_W", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! Z LOCATION (METERS)
      scalar = visitmdscalarcreate(handle, "LAG_Z", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! SOME SCALAR
      scalar = visitmdscalarcreate(handle, "LAG_S", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)

      ! PATH LENGTH
      scalar = visitmdscalarcreate(handle, "LAG_D", 5, "LAGRANGIAN_MESH", 15, VISIT_VARCENTERING_NODE)


   end if
# endif

  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! ADD SCALAR METADATA
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  
  !  IF (VISIT_OPT .eq. 'basic') THEN
  scalar = visitmdscalarcreate(handle, "H", 1, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "EL", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
  
  scalar = visitmdscalarcreate(handle, "U", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "V", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "WW", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  

  scalar = visitmdscalarcreate(handle, "UA", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "VA", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)

  scalar = visitmdscalarcreate(handle, "T1", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "S1", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "RHO1", 4, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
 
  scalar = visitmdscalarcreate(handle, "VX", 2, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "VY", 2, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
 
  
#  if defined (MULTIPROCESSOR)
  scalar = visitmdscalarcreate(handle, "EL_PID", 6, "TWOD_Mesh", 9,&
       & VISIT_VARCENTERING_ZONE)


  if(PAR)then
     ! FVCOM NODE INDEX VARIABLES
     scalar = visitmdscalarcreate(handle, "Node_Index_TWOD", 15, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_Bathy", 16, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_SSH", 14, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_SigmaLevel", 21, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Node_Index_SigmaLayer", 21, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     
     ! FVCOM CELL INDEX VARIABLES
     scalar = visitmdscalarcreate(handle, "Cell_Index_TWOD", 15, "TWOD_Mesh", 9, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_Bathy", 16, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_SSH", 14, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_SigmaLevel", 21, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "Cell_Index_SigmaLayer", 21, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_ZONE)
  end if
# endif
  

# if defined (SPHERICAL)

  scalar = visitmdscalarcreate(handle, "XYZ_U", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_V", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_WW", 6, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)


  IF (VISIT_OPT .eq. 'advanced') THEN 
  scalar = visitmdscalarcreate(handle, "XYZ_UA", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_VA", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_WA", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  END IF

# endif

  
  IF (VISIT_OPT .eq. 'advanced') THEN 

     ! grid metrics
     scalar = visitmdscalarcreate(handle, "ART", 3, "TWOD_Mesh", 9, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ART1", 4, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "ART2", 4, "TWOD_Mesh", 9, VISIT_VARCENTERING_NODE)

     ! 2-d Arrays for the general vertical coordinates
     scalar = visitmdscalarcreate(handle, "Z", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "ZZ", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "DZ", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     ! Can not plot DZZ, it has KBM2 node centered values. There is
     ! no mesh for that!
     !scalar = visitmdscalarcreate(handle, "DZZ", 3, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)

     ! 2-d flow variable arrays at elements
     scalar = visitmdscalarcreate(handle, "UARK", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "VARK", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)

     scalar = visitmdscalarcreate(handle, "COR", 3, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !!$ DO NOT INCLUDE INTERPOLATED VALUES AT THE CELL CENTERS
     !  scalar = visitmdscalarcreate(handle, "H1", 2, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "D1", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "EL1", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "DT1", 3, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "ET1", 3, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     !  scalar = visitmdscalarcreate(handle, "ELRK1", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DTFA", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "CC_SPONGE", 9, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     ! 2-d flow variable arrays at nodes
     ! INCLUDED IN GENERAL VARIALBES
     scalar = visitmdscalarcreate(handle, "D", 1, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "DT", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "ET", 2, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "EGF", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "ELRK", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)

     ! SURFACE AND BOTTOM BOUNDARY CONDITIONS
     scalar = visitmdscalarcreate(handle, "CBC", 3, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "SWRAD", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "WUSURF2", 7, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WVSURF2", 7, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WUBOT", 5, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WVBOT", 5, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WUSURF", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WVSURF", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     ! Skipping WTSURF: What is it?
     ! GROUND WATER DISCHARGE
     scalar = visitmdscalarcreate(handle, "BFWDIS", 6, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "BFWDIS2", 7, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "BFWDIS3", 7, "Bathymetry_Mesh", 15, VISIT_VARCENTERING_NODE)
     ! RIVER DISCHARGE: add later!

     ! SURFACE MET FORCING
     scalar = visitmdscalarcreate(handle, "UUWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "VVWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# if defined (SPHERICAL)
     scalar = visitmdscalarcreate(handle, "XYZ_UWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "XYZ_VWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "XYZ_WWIND", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# endif
     scalar = visitmdscalarcreate(handle, "QPREC2", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "QPREC3", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "QEVAP2", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "QEVAP3", 6, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)

     ! 2 D Flow FLuxes
     scalar = visitmdscalarcreate(handle, "PSTX", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "PSTY", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVUA", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVVA", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADX2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADY2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DRX2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DRY2D", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVX", 4, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "ADVY", 4, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)

     ! INTERNAL MODE ARRAYS- ELEMENTS
     scalar = visitmdscalarcreate(handle, "W", 1, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_ZONE)

     scalar = visitmdscalarcreate(handle, "UF", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "VF", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "WT", 2, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_ZONE)

     ! INTERNAL MODE ARRAYS- NODES
     scalar = visitmdscalarcreate(handle, "Q2", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "L", 1, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Q2L", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "Q2", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)

# if defined (GOTM)
     scalar = visitmdscalarcreate(handle, "TKE", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "TKEF", 4, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "TEPS", 4, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "TEPSF", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
# endif

      !--------------------------------------------------------------
      !--------------------------------------------------------------
# if defined (ICE)

  scalar = visitmdscalarcreate(handle, "AICE", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "VICE", 4, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)
  scalar = visitmdscalarcreate(handle, "STRENGTH", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_NODE)

  scalar = visitmdscalarcreate(handle, "UICE2", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "VICE2", 5, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)

# if defined (SPHERICAL)
  scalar = visitmdscalarcreate(handle, "XYZ_UICE", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_VICE", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
  scalar = visitmdscalarcreate(handle, "XYZ_WICE", 8, "SSH_Mesh", 8, VISIT_VARCENTERING_ZONE)
# endif

# endif
     !--------------------------------------------------------------
     !--------------------------------------------------------------



     scalar = visitmdscalarcreate(handle, "KM", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "KH", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "KQ", 2, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "AAM", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)

     scalar = visitmdscalarcreate(handle, "TF1", 3, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "SF1", 3, "SigmaLayer_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "WTS", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)
     scalar = visitmdscalarcreate(handle, "WTTS", 3, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_NODE)

     ! Baroclinic Pressure Gradients
     scalar = visitmdscalarcreate(handle, "DRHOX", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)
     scalar = visitmdscalarcreate(handle, "DRHOY", 5, "SigmaLevel_Mesh", 15, VISIT_VARCENTERING_ZONE)

     ! Shape Coefficients: add later if needed

     ! Salinity and temperature bottom diffusion stuff

     ! Average arrays ?

     ! WHAT IS VISCOFH


     ! WHAT ABOUT BOUNDARY CONDITION VARIABLES?


  END IF

  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! VISIT MATERIALS
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  ! VISIT Materials are not working yet. Not sure Why? Causes Sig Sev...

  ! Add sigma layers material
  mat1 = visitmdmaterialcreate(handle, "Sigma_Layers", 12, "SigmaLayer_Mesh", 15)
  if(mat1.ne.VISIT_INVALID_HANDLE) then
     
     Do ind=1,KBM2
        
        write(ch3,'(i3)') ind
        
        if(ind<10 .and. ind>=1) then
           matname="Layer_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Layer_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Layer_"//trim(adjustl(ch3))
        else
           if(MSR) WRITE(VDB,*)"In VisitGetMetaData: Layer Number > 1000"
           if(MSR) WRITE(VDB,*)"Causing error in visit material"
        end if
        
        err = visitmdmaterialadd(handle, mat1, matname, 9)
        
     End Do
  endif


  mat2 = visitmdmaterialcreate(handle, "Sigma_Levels", 12, "SigmaLevel_Mesh", 15)
  if(mat2.ne.VISIT_INVALID_HANDLE) then
     
     Do ind=1,KBM1
        
        write(ch3,'(i3)') ind

        if(ind<10 .and. ind>=1) then
           matname="Level_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Level_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Level_"//trim(adjustl(ch3))
        else
           if(msr) WRITE(VDB,*)"In VisitGetMetaData: Level Number > 1000"
           if(msr) WRITE(VDB,*)"Causing error in visit material"
        end if

        err = visitmdmaterialadd(handle, mat2, matname, 9)


     End Do
  endif


  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! VISIT CONTROL COMMANDS
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  err = visitmdaddsimcommand(handle, "Halt/Step", 9, VISIT_CMDARG_NONE,1)

  err = visitmdaddsimcommand(handle, "Run", 3, VISIT_CMDARG_NONE,1)

  err = visitmdaddsimcommand(handle, "Dump Restart", 12, VISIT_CMDARG_NONE,1)

  err = visitmdaddsimcommand(handle, "Toggle Step", 11, VISIT_CMDARG_NONE,1)


!  err = visitmdaddsimcommand(handle, "i10", 3, VISIT_CMDARG_NONE,1)
!  err = visitmdaddsimcommand(handle, "i100", 4, VISIT_CMDARG_NONE,1)




  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! FVCOM EXPRESSIONS
  !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# if defined (SPHERICAL)

  err = visitmdexpressioncreate(handle, "Velocity",8,&
       & "{XYZ_U, XYZ_V, XYZ_WW}", 22, VISIT_VARTYPE_VECTOR) 

  err = visitmdexpressioncreate(handle, "Average_Velocity",16,&
       & "{XYZ_UA, XYZ_VA, XYZ_WA}", 24, VISIT_VARTYPE_VECTOR) 

  if(MSR) WRITE(VDB,*)"Expression create returned: ", err

# if defined (ICE)
  err = visitmdexpressioncreate(handle, "Ice_Velocity",12,&
       & "{XYZ_UICE, XYZ_VICE, XYZ_WICE}", 30, VISIT_VARTYPE_VECTOR) 
# endif

  err = visitmdexpressioncreate(handle, "Wind_Velocity",13,&
       & "{XYZ_UWIND, XYZ_VWIND, XYZ_WWIND}", 33, VISIT_VARTYPE_VECTOR) 

# else  
! if not spherical !!

  err = visitmdexpressioncreate(handle, "Velocity",8, "{U,V,WW}",&
       & 8, VISIT_VARTYPE_VECTOR) 

  err = visitmdexpressioncreate(handle, "Average_Velocity",16, "{UA,VA,0}",&
       & 9, VISIT_VARTYPE_VECTOR) 

# if defined (ICE)
  err = visitmdexpressioncreate(handle, "Ice_Velocity",12,&
       & "{UICE2,VICE2, 0}", 16, VISIT_VARTYPE_VECTOR) 
# endif

  err = visitmdexpressioncreate(handle, "Wind_Velocity",13,&
       & "{UUWIND,VVWIND,0}", 17, VISIT_VARTYPE_VECTOR) 

# endif


  visitgetmetadata = VISIT_OKAY


end function visitgetmetadata

!c---------------------------------------------------------------------------
!c visitgetmesh
!c---------------------------------------------------------------------------
integer function visitgetmesh(handle, domain, name, lname)
  USE MOD_VISIT
  USE ALL_VARS
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(in) :: name
  integer, intent(in) :: handle, domain, lname
  integer :: error
  TYPE(VISITMESHTYPE), POINTER :: vmp ! VISIT MESH POINTER

  TYPE(VISITLAG), POINTER :: lmp ! VISIT MESH POINTER

  
  error = VISIT_ERROR

  if (domain .neqv. VisitParRank) then
     if(MSR) WRITE(VDB,*) "Something screwy is going on: visit domain= ",domain, "&
          &; VisitParRank= ",VisitParRank 
     visitgetmesh=error
     return
  end if



  if (visitstrcmp(name, lname, "TWOD_Mesh", 9).eq.0) then
     
     vmp => VISIT_MESH(VISIT_TWODMESH)
     Call UPDATE_MESH(VMP,error)
     if (error .ne. 0) then ! CHECK FOR BAD RESULT
        visitgetmesh=NNode
        return
     end if

    
  elseif (visitstrcmp(name, lname, "Bathymetry_Mesh",15).eq.0) then
     
     vmp => VISIT_MESH(VISIT_BATHYMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=NNode
        return        
     end if
     
  elseif (visitstrcmp(name, lname, "SSH_Mesh",8).eq.0) then
     
     vmp => VISIT_MESH(VISIT_SSHMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=NNode
        return        
     end if
     
     
  elseif (visitstrcmp(name, lname, "SigmaLayer_Mesh",15).eq.0) then


     vmp => VISIT_MESH(VISIT_LAYERMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=NNode
        return        
     end if
     
     
  elseif (visitstrcmp(name, lname, "SigmaLevel_Mesh",15).eq.0) then
     
     vmp => VISIT_MESH(VISIT_LEVELMESH)
     Call UPDATE_MESH(VMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=NNode
        return        
     end if

#  if defined (NETCDF_IO)
  elseif (visitstrcmp(name, lname, "LAGRANGIAN_MESH", 15).eq.0) then
     
     lmp => VISIT_LAGDATA
     LMP%NDIMS=3
     Call UPDATE_LAG(LMP,error)
     if (error == -1) then ! CHECK FOR BAD RESULT
        visitgetmesh=NNode
        return        
     end if
     
     visitgetmesh = visitmeshpoint(handle, LMP%NDIMS, LMP%NODES, LMP%VX, LMP%VY, LMP%VZ)

     NULLIFY(LMP)
     return
     
# endif

  end if
  
   

  error = visitmeshunstructured3(handle,VMP%NDIMS, VMP%Nodes,&
       & VMP%Zones,VMP%GHOSTZONES, VMP%VX,&
       & VMP%VY, VMP%VZ, VMP%LCONN, VMP%NV&
       &,VMP%DATAOWNER)
  
  NULLIFY(VMP)

  visitgetmesh = error
end function visitgetmesh


!c---------------------------------------------------------------------------
!c visitgetscalar
!c---------------------------------------------------------------------------
integer function visitgetscalar(handle, domain, name, lname)
  USE ALL_VARS
  USE MOD_VISIT
  USE BCS
  USE MOD_PAR

      !--------------------------------------------------------------
      !--------------------------------------------------------------
# if defined (ICE)
      use ice_state
      use mod_ice
      use mod_ice2d
# endif
     !--------------------------------------------------------------
     !--------------------------------------------------------------


  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8),intent(in) :: name
  character(LEN=100) :: fnameVVV
  integer, intent(in) :: handle, domain, lname
  integer :: error
  integer, dimension(3) :: Vdims
  integer :: Vsize, ind


  integer, allocatable, dimension(:) :: V1di
  real*4, allocatable, dimension(:) :: V1df

  TYPE(VisitSphereVel), pointer :: VSV

  error = VISIT_ERROR


  ! BASIC Components of FVCOM 
  if (visitstrcmp(name, lname, "H",1).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, H(1:NTNode),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "EL",2).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, EL(1:NTNode),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "U",1).eq.0) then
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, U(1:MTElem,1:KBM1),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "V",1).eq.0) then
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V(1:MTElem,1:KBM1),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WW",2).eq.0) then
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WW(1:MTElem,1:KBM1),Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "T1",2).eq.0) then
     Vsize=KBM1*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, T1(1:NTNode,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "S1",2).eq.0) then
     Vsize=KBM1*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, S1(1:NTNode,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "RHO1",4).eq.0) then
     Vsize=KBM1*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, RHO1(1:NTNode,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VX",2).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VX(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VY",2).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VY(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# if defined(MULTIPROCESSOR)
  elseif (visitstrcmp(name, lname, "EL_PID",6).eq.0) then

     allocate(v1di(MTElem)); v1di=myid
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_TWOD",15).eq.0) then

     allocate(v1di(NTNode)); v1di=0
     v1di(1:NNode)=NGID(1:NNode)
     
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_Bathy",16).eq.0) then

     allocate(v1di(NTNode)); v1di=0
     v1di(1:NNode)=NGID(1:NNode)
     
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_SSH",14).eq.0) then

     allocate(v1di(NTNode)); v1di=0
     v1di(1:NNode)=NGID(1:NNode)
     
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_SigmaLevel",21).eq.0) then

     allocate(v1di(NTNode*KB)); v1di=0
     Do ind=1,KB
        v1di(((ind-1)*NTNode+1):((ind-1)*NTNode+NNode))=NGID(1:NNode)
     End Do
     
     Vsize=NTNode*KB
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Node_Index_SigmaLayer",21).eq.0) then

     allocate(v1di(NTNode*KBM1)); v1di=0
     Do ind=1,KBM1
        v1di(((ind-1)*NTNode+1):((ind-1)*NTNode+NNode))=NGID(1:NNode)
     End Do
     
     Vsize=NTNode*KBM1
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)


     ! Cell Index arrays!
  elseif (visitstrcmp(name, lname, "Cell_Index_TWOD",15).eq.0) then

     allocate(v1di(MTElem)); v1di=0
     v1di(1:MElem)=EGID(1:MElem)
     
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_Bathy",16).eq.0) then

     allocate(v1di(MTElem)); v1di=0
     v1di(1:MElem)=EGID(1:MElem)
     
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_SSH",14).eq.0) then

     allocate(v1di(MTElem)); v1di=0
     v1di(1:MElem)=EGID(1:MElem)
     
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_SigmaLevel",21).eq.0) then

     allocate(v1di(MTElem*KBM1)); v1di=0
     Do ind=1,KBM1
        v1di(((ind-1)*MTElem+1):((ind-1)*MTElem+MElem))=EGID(1:MElem)
     End Do
     
     Vsize=MTElem*KB
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

  elseif (visitstrcmp(name, lname, "Cell_Index_SigmaLayer",21).eq.0) then

     allocate(v1di(MTElem*KBM1)); v1di=0
     Do ind=1,KBM2
        v1di(((ind-1)*MTElem+1):((ind-1)*MTElem+MElem))=EGID(1:MElem)
     End Do
     
     Vsize=MTElem*KBM2
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, v1di,Vdims , 1 &
          &,VISIT_DATATYPE_INT, VISIT_OWNER_VISIT)
     deallocate(v1di)

# endif


# if defined (SPHERICAL)

! Velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_U",5).eq.0) then
     VSV =>VisitSphericalVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereVel(VSV)
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_V",5).eq.0) then
     VSV =>VisitSphericalVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereVel(VSV)
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WW",6).eq.0) then
     VSV =>VisitSphericalVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereVel(VSV)
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)


 !! Average velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_UA",6).eq.0) then
     VSV =>VisitSphericalAVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereAVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_VA",6).eq.0) then
     VSV =>VisitSphericalAVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereAVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WA",6).eq.0) then
     VSV =>VisitSphericalAVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereAVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

! WIND Velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_UWIND",9).eq.0) then
     VSV =>VisitSphericalWindVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereWindVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_VWIND",9).eq.0) then
     VSV =>VisitSphericalWindVel
     if (VSV%Updated_Time .NE. Visit_Time_INT) &
          & Call UpdateSphereWindVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WWIND",9).eq.0) then
     VSV =>VisitSphericalWindVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereWindVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

# if defined (ICE)
! ICE Velocity terms
  elseif (visitstrcmp(name, lname, "XYZ_UICE",8).eq.0) then
     VSV =>VisitSphericalIceVel
     if (VSV%Updated_Time .NE. Visit_Time_Int) &
          & Call UpdateSphereICEVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelX,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_VICE",9).eq.0) then
     VSV =>VisitSphericalIceVel
     if (VSV%Updated_Time .NE. Visit_Time_INT) &
          & Call UpdateSphereIceVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelY,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

  elseif (visitstrcmp(name, lname, "XYZ_WICE",9).eq.0) then
     VSV =>VisitSphericalIceVel
     if (VSV%Updated_Time .NE. Visit_Time_Ext) &
          & Call UpdateSphereIceVel(VSV)
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VSV%VelZ,Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_SIM)
     nullify(VSV)

# endif


# endif


#  if defined (NETCDF_IO)

     ! LAGRANGIAN PARTICLE TRACKING STUFF GOES HERE!
  elseif (visitstrcmp(name, lname, "LAG_U",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VU, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "LAG_V",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VV, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "LAG_W",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VW, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "LAG_Z",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VZ, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "LAG_S",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VS, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


  elseif (visitstrcmp(name, lname, "LAG_D",5).eq.0) then
     Vsize=Visit_LagData%Nodes
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Visit_LagData%VD, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# endif

     !=================================================================
     !=================================================================
     !=================================================================
     ! END BASIC VARIABLES   THE REST ARE VARIABLES ACCESSED BY RUNNING 
     ! FVCOM WITH THE VISIT_OPT FLAG SET TO ADVANCED. THESE VARIABLES 
     ! ARE FOR ADVANCDED USERS WHO KNOW THEIR USAGE AND MEANING IN FVCOM


     ! GRID METRICS
  elseif (visitstrcmp(name, lname, "ART",3).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ART(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ART1",4).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ART1(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ART2",4).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ART2(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     ! 2-d Arrays for the general vertical coordinates
  elseif (visitstrcmp(name, lname, "Z",1).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Z(1:NTNode,1:KB), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ZZ",2).eq.0) then
     Vsize=KBM1*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ZZ(1:NTNode,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DZ",2).eq.0) then
     Vsize=KBM1*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DZ(1:NTNode,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     !  elseif (visitstrcmp(name, lname, "DZZ",3).eq.0) then
     !     Vsize=KBM1*NTNode
     !     Vdims=(/Vsize,1,1/)
     !     error  = visitscalarsetdata(handle, DZZ(1:NTNode,1:KBM1), Vdims , 1 &
     !          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     ! 2-d flow variable arrays at elements

  elseif (visitstrcmp(name, lname, "UA",2).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UA(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VA",2).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VA(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "UARK",4).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UARK(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VARK",4).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UARK(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "COR",3).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, COR(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DTFA",4).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DTFA(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "CC_SPONGE",9).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, COR(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! 2-d flow variable arrays at nodes
  elseif (visitstrcmp(name, lname, "D",1).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, D(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DT",2).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DT(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ET",2).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ET(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "EGF",3).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, EGF(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ELRK",4).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ELRK(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

     ! SURFACE AND BOTTOM BOUNDARY CONDITIONS

  elseif (visitstrcmp(name, lname, "CBC",3).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, CBC(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "SWRAD",5).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, SWRAD(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WUSURF2",7).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WUSURF2(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WVSURF2",7).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WVSURF2(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WUBOT",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WUBOT(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WVBOT",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WVBOT(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WUSURF",6).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WUSURF(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WVSURF",6).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WVSURF(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! Ground water variables!!!
     ! BFWDIS variables are allocated BFWDIS(IBFW) 
         ! => to nodes mesh!
  elseif (visitstrcmp(name, lname, "BFWDIS",6).eq.0) then
     allocate(V1df(NTNode)); V1df=0
     V1df(NODE_BFW)=BFWDIS
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V1df, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     Deallocate(V1df)

  elseif (visitstrcmp(name, lname, "BFWDIS2",7).eq.0) then
     allocate(V1df(NTNode)); V1df=0
     V1df(NODE_BFW)=BFWDIS2
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V1df, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     Deallocate(V1df)

  elseif (visitstrcmp(name, lname, "BFWDIS3",7).eq.0) then
     allocate(V1df(NTNode)); V1df=0
     V1df(NODE_BFW)=BFWDIS3
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, V1df, Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
     Deallocate(V1df)
     ! END of Ground water variables!

     ! SURFACE MET FORCING
  elseif (visitstrcmp(name, lname, "UUWIND",6).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UUWIND(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VVWIND",6).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VVWIND(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QPREC2",6).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QPREC2(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QPREC3",6).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QPREC3(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QEVAP2",6).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QEVAP2(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "QEVAP3",6).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, QEVAP3(1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "PSTX",4).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, PSTX(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "PSTY",4).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, PSTY(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVUA",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVUA(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVVA",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVVA(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADX2D",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADX2D(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADY2D",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADY2D(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRX2D",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRX2D(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRY2D",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRY2D(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVX",4).eq.0) then
     Vsize=MTElem*KBM1
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVX(1:MTElem,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "ADVY",4).eq.0) then
     Vsize=MTElem*KBM1
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, ADVY(1:MTElem,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! INTERNAL MODE ARRAYS- ELEMENTS -- U, V, WW are in the basic
  elseif (visitstrcmp(name, lname, "W",1).eq.0) then
     Vsize=KBM2*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, W(1:MTElem,2:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "UF",2).eq.0) then
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VF(1:MTElem,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VF",2).eq.0) then
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VF(1:MTElem,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WT",2).eq.0) then
     Vsize=KBM2*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WT(1:MTElem,2:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


     ! INTERNAL MODE ARRAYS- NODES
  elseif (visitstrcmp(name, lname, "Q2",2).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Q2(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "L",1).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, L(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "Q2L",3).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, Q2L(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# if defined (GOTM)
  elseif (visitstrcmp(name, lname, "TKE",3).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TKE(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "TKEF",4).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TKEF(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "TEPS",4).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TEPS(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "TEPSF",5).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TEPSF(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)
# endif

  elseif (visitstrcmp(name, lname, "KM",2).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, KM(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "KH",2).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, KH(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "KQ",2).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, KQ(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "AAM",3).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, AAM(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "TF1",3).eq.0) then
     Vsize=KBM1*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, TF1(1:NTNode,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "SF1",3).eq.0) then
     Vsize=KBM1*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, SF1(1:NTNode,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WTS",3).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WTS(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WTTS",4).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WTTS(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "WTS",3).eq.0) then
     Vsize=KB*NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, WTS(1:NTNode,:), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRHOX",5).eq.0) then
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRHOX(1:MTElem,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "DRHOY",5).eq.0) then
     Vsize=KBM1*MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, DRHOY(1:MTElem,1:KBM1), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)


      !--------------------------------------------------------------
      !--------------------------------------------------------------
# if defined (ICE)
  elseif (visitstrcmp(name, lname, "AICE",4).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, AICE(1,1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VICE",4).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, AICE(1,1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "STRENGTH",8).eq.0) then
     Vsize=NTNode
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, STRENGTH(1,1:NTNode), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "UICE2",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, UICE2(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

  elseif (visitstrcmp(name, lname, "VICE2",5).eq.0) then
     Vsize=MTElem
     Vdims=(/Vsize,1,1/)
     error  = visitscalarsetdata(handle, VICE2(1:MTElem), Vdims , 1 &
          &,VISIT_DATATYPE_FLOAT, VISIT_OWNER_VISIT)

# endif
     !--------------------------------------------------------------
     !--------------------------------------------------------------

  endif


  visitgetscalar = error
end function visitgetscalar


!c---------------------------------------------------------------------------
!c visitgetcurve
!c---------------------------------------------------------------------------
integer function visitgetcurve(handle, name, lname)
  implicit none
  include "visitfortransiminterface.inc"
  character*8 name
  integer     handle, lname

  visitgetcurve = VISIT_ERROR
end function visitgetcurve

!c---------------------------------------------------------------------------
!c visitgetdomainlist
!c---------------------------------------------------------------------------
integer function visitgetdomainlist(handle)
  USE MOD_VISIT, only: VisitParRank
  USE LIMS, only : NPROCS
  implicit none
  include "visitfortransiminterface.inc"
  integer handle

  visitgetdomainlist = visitsetdomainlist(handle, NPROCS, VisitParRank,1)

  visitgetdomainlist = VISIT_OKAY
end function visitgetdomainlist



!c---------------------------------------------------------------------------
!c visitgetmaterial
!c---------------------------------------------------------------------------
integer function visitgetmaterial(handle, domain, name, lname)
  USE LIMS
  USE MOD_VISIT
  implicit none
  include "visitfortransiminterface.inc"
  character(LEN=8), intent(in) :: name
  integer, intent(in) :: handle, domain, lname
  integer:: rslt, ind, ind2
  character(LEN=3) :: ch3
  character(LEN=20) :: matname
  INTEGER, allocatable :: MatNo(:)
  integer, dimension(3):: dims

  if (visitstrcmp(name, lname, "Sigma_Layers",12).eq.0) then
     
     allocate(MatNo(KBM2))

     dims(1)=MTElem*KBM2
     dims(2)=1
     dims(3)=1

     rslt = visitMaterialSetDims(handle, dims, 1)
     Do ind =1, KBM2
        write(ch3,'(i3)') ind
        
        if(ind<10 .and. ind>=1) then
           matname="Layer_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Layer_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Layer_"//trim(adjustl(ch3))
        else
           if(MSR) WRITE(VDB,*)"In GetMaterial: Layer Number > 1000"
           if(MSR) WRITE(VDB,*)"Causing error in visit material"
        end if

        MatNo(ind) = visitMaterialAdd(handle, matname, 9)
        
     End Do
     
     Do ind=1, KBM2
        Do ind2=1,MTElem
           rslt = visitMaterialAddClean(handle,(ind-1)*MTElem+ind2,MatNo(ind))
        End Do
     End Do

     Deallocate(MatNo)
     rslt = VISIT_OKAY


  elseif (visitstrcmp(name, lname, "Sigma_Levels",12).eq.0) then
     
     allocate(MatNo(KBM1))

     dims(1)=MTElem*KBM1
     dims(2)=1
     dims(3)=1
     
     rslt = visitMaterialSetDims(handle, dims, 1)
     
     Do ind =1, KBM1

        write(ch3,'(i3)') ind
        
        if(ind<10 .and. ind>=1) then
           matname="Level_00"//trim(adjustl(ch3))
        else if(ind<100 .and. ind>=10) then
           matname="Level_0"//trim(adjustl(ch3))
        else if(ind<1000 .and. ind>=100) then
           matname="Level_"//trim(adjustl(ch3))
        else
           if(MSR) WRITE(VDB,*)"In GetMaterial: Level Number > 1000"
           if(MSR) WRITE(VDB,*)"Causing error in visit material"
        end if

        MatNo(ind) = visitMaterialAdd(handle, matname, 9)
     End Do
     

     Do ind=1, KBM1
        Do ind2=1,MTElem
           rslt = visitMaterialAddClean(handle,(ind-1)*MTElem+ind2,MatNo(ind))
        End Do
     End Do

     Deallocate(MatNo)
     rslt = VISIT_OKAY     


  else
     rslt = VISIT_ERROR
  end if
  


  visitgetmaterial = rslt
end function visitgetmaterial


#else
! if visit is not defined compile a dummy subroutine!

subroutine visit_dummy

implicit none

end subroutine visit_dummy

#endif

!end MODULE MOD_VISIT

























