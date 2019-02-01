#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/bash
            DEF_FLAGS     = -P  -traditional 
#           EXEC          = fvcom_serial
	    EXEC          = chesfvm_HFX
#==========================================================================
#  BEGIN USER DEFINITION SECTION
#==========================================================================
#        SELECT MODEL OPTIONS
#          SELECT FROM THE FOLLOWING OPTIONS BEFORE COMPILING CODE
#          SELECT/UNSELECT BY COMMENTING/UNCOMMENTING LINE (#)
#          CODE MUST BE CLEANED (with "make clean") AND
#          RECOMPILED IF NEW SET OF OPTIONS IS DESIRED
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#        PRECISION          DEFAULT PRECISION: SINGLE                     
#                           UNCOMMENT TO SELECT DOUBLE PRECISION
#--------------------------------------------------------------------------

#            FLAG_1 =  -DDOUBLE_PRECISION 

#--------------------------------------------------------------------------
#        SPHERICAL          SELECT SPHERICAL COORDINATES FOR INTEGRATION
#                           DEFAULT: CARTESIAN
#                           UNCOMMENT TO SELECT SPHERICAL COORDINATES
#--------------------------------------------------------------------------

#             FLAG_2 = -DSPHERICAL

#--------------------------------------------------------------------------
#        FLOODYING/DRYING   INCLUDE WET/DRY TREATMENT OF DOMAIN
#                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH
#                           INPUT FILE CONTROL.  (SEE exa_run.dat) FILE
#                           DEFAULT: NO FLOODYING/DRYING INCLUDED
#                           UNCOMMENT TO INCLUDE FLOODYING/DRYING
#--------------------------------------------------------------------------

             FLAG_3 = -DWET_DRY

#--------------------------------------------------------------------------
#        MULTI_PROCESSOR    INCLUDES PARALLELIZATION WITH MPI
#                           REQUIRES LINKING MPI LIBRARIES OR COMPILING
#                           WITH A PRELINKED SCRIPT (mpif90/mpf90/etc)
#                           DEFAULT: NO PARALLEL CAPABILITY
#                           UNCOMMENT TO INCLUDE MPI PARALLEL CAPABILITY
#--------------------------------------------------------------------------

              FLAG_4 = -DMULTIPROCESSOR


# For Spartina at UMCES
              PARLIB = -L/data/users/bclark/software/METIS_source -lmetis



#--------------------------------------------------------------------------
#        WATER_QUALITY      INCLUDE EPA WATER QUALITY MOD
#                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH
#                           VARIABLE WQM_ON IN INPUT FILE
#                           DEFAULT: NO WATER QUALITY MODEL
#                           UNCOMMENT TO INCLUDE WATER QUALITY MODEL
#--------------------------------------------------------------------------

#             FLAG_5 = -DWATER_QUALITY

#--------------------------------------------------------------------------
#        NETCDF OUTPUT      DUMP OUTPUT INTO NETCDF FILES (yes/no)
#                           REQUIRES SYSTEM DEPENDENT NETCDF LIBRARIES
#                           COMPILED WITH SAME F90 COMPILER
#                           SET PATH TO LIBRARIES WITH IOLIBS      
#                           SET PATH TO INCLUDE FILES (netcdf.mod) WITH IOINCS
#                           DEFAULT: NO NETCDF OUTPUT
#                           UNCOMMENT TO INCLUDE NETCDF OUTPUT CAPABILITY
#--------------------------------------------------------------------------
           FLAG_6       =  -DNETCDF_IO

# For Spartina at UMCES
           IOLIBS       =  -L/data/users/bclark/software/netcdf-3.6.1-src/intel/lib -lnetcdf
           IOINCS       =  -I/data/users/bclark/software/netcdf-3.6.1-src/intel/include

#             IOLIBS       =  -L/usr/local/install/netcdf/gcc_ifort/3.6.2/lib  -lnetcdf
#             IOINCS       =  -I/usr/local/install/netcdf/gcc_ifort/3.6.2/include

#--------------------------------------------------------------------------
#       NG_OI_ASSIMILATION  INCLUDE NUDGING BASED DATA ASSIMILATION AND 
#                           OPTIMAL INTERPOLATION(OI) BASED DATA 
#                           ASSIMILATION FOR CURRENT/TEMP/SALINITY/SST
#                           CAN BE ACTIVATED/DEACTIVATED AT RUN TIME WITH 
#                           INPUT FILE CONTROL.  (SEE exa_run.dat) FILE 
#                           DEFAULT: NO DATA ASSIMILATION INCLUDED 
#                           UNCOMMENT TO INCLUDE DATA ASSIMILATION 
#--------------------------------------------------------------------------

#             FLAG_7 = -DNG_OI_ASSIM

# For Build on em64t computer (Guppy) 
#             OILIB  = -L/usr/X11R6/lib64 -L/usr/local/lib64 -lmkl_lapack -lmkl_em64t -lguide -lpthread

# For Build on Cluster (Typhoeus and Hydra)
#             OILIB  = -L/usr/X11R6/lib64 -L/usr/local/lib/em64t -lmkl_lapack -lmkl_em64t -lguide -lpthread

#--------------------------------------------------------------------------
#       REMOVE SYSTEM CALLS IF RUNNING IN A DOS/WINDOWS ENVIRONMENT      
#                           OR IF COMPILER FAILS WHEN COMPILING SYSTEM
#                           CALLS, TURN ON THIS OPTION.  
#                           THIS WILL DEACTIVATE CERTAIN COMMANDS WHICH
#                           SETUP THE OUTPUT DIRECTORY STRUCTURE
#                           NOTE THAT THE USER WILL HAVE TO CREATE THE
#                           OUTPUT DIRECTORIES MANUALLY
#                           SEE THE README FILE INCLUDED WITH THE 
#                           DISTRIBUTION
#                           DEFAULT: INCLUDE SYSTEM CALLS
#                           UNCOMMENT TO SUPRESS SYSTEM CALLS
#--------------------------------------------------------------------------

#             FLAG_8 += -DDOS

#---------------------------------------------------------------------------
#       SOLID BOUNDARY     IF GCN, NO GHOST CELL
#                          IF GCY1, GHOST CELL IS SYMMETRIC RELATIVE TO BOUNDARY
#                          CELL EDGE
#                          IF GCY2, GHOST CELL IS SYMMETRIC RELATIVE TO MIDDLE
#                          POINT OF THE BOUNDARY CELL EDGE
#   !!!!!! ONLY ONE OF THE FLAGS BELOW CAN BE AND MUST BE CHOSEN
#---------------------------------------------------------------------------

             FLAG_9  = -DGCN
#             FLAG_9  = -DGCY1
#             FLAG_9 = -DGCY2

#--------------------------------------------------------------------------
#        TURBULENCE MODEL   USE GOTM TURBULENCE MODEL INSTEAD OF THE ORIGINAL
#                           FVCOM MELLOR-YAMADA 2.5 IMPLEMENTATION
#                           UNCOMMENT TO USE GOTM TURBULENCE MODEL
#--------------------------------------------------------------------------

#            FLAG_10 = -DGOTM -DDOUBLE_PRECISION
#            GOTMLIB       = -L../GOTM_source -lturbulence -lutil -lmeanflow
#            GOTMINCS      = -I../GOTM_source  

#--------------------------------------------------------------------------
#        pV3 Realtime Visualization Server (MEDM GROUP ONLY)                   
#--------------------------------------------------------------------------

#            FLAG_11 = -DPV3  
#            PV3LIB   = -L/hosts/salmon01/data00/medm/lib -lpV3 -lgpvm3 -lpvm3 -lfpvm3 -lpgftnrtl  -lpgc 

#--------------------------------------------------------------------------
#       EQUILIBRIUM TIDE
#--------------------------------------------------------------------------

#            FLAG_12 = -DEQUI_TIDE

#--------------------------------------------------------------------------
#       ATMOSPHERIC TIDE
#--------------------------------------------------------------------------

#            FLAG_13 = -DATMO_TIDE

#--------------------------------------------------------------------------
#       ARCTIC OCEAN INCLUDED (If you chose this flag, FLAG_2 should be 
#                              selected)
#--------------------------------------------------------------------------

#             FLAG_14 = -DNORTHPOLE

#--------------------------------------------------------------------------
#        Using A fully multidimensional positive definite advection
#        transport algorithm with small implicit diffusion. 
#        Based on Smolarkiewicz, P. K; Journal of Computational
#        Physics, 54, 325-362, 1984
#--------------------------------------------------------------------------

#            FLAG_15 = -DMPDATA

#--------------------------------------------------------------------------
#         Run Two-D Barotropic Mode Only
#--------------------------------------------------------------------------

#             FLAG_16 = -DTWO_D_MODEL

#--------------------------------------------------------------------------
#         Output 2-D Momentum Balance Checking
#--------------------------------------------------------------------------

#            FLAG_17 = -DBALANCE_2D  

#--------------------------------------------------------------------------
#           open boundary T/S time series nudging
#---------------------------------------------------------------------------

            FLAG_18 = -DTS_OBC

#--------------------------------------------------------------------------
#           OPEN BOUNDARY FORCING TYPE 
#           DEFAULT: OPEN BOUNDARY NODE WATER ELEVATION FORCING
#           UNCOMMENT TO SELECT BOTH OPEN BOUNDARY NODE WATER ELEVATION
#           FORCING AND OPEN BOUNDARY VOLUME TRANSPORT FORCING 
#---------------------------------------------------------------------------

#            FLAG_19 = -DMEAN_FLOW

#--------------------------------------------------------------------------
#           OUTPUT TIDAL INFORMATION AT NTIDENODE and NTIDECELL
#           FOR MEANFLOW CALCULATION.
#---------------------------------------------------------------------------

#             FLAG_20 = -DTIDE_OUTPUT

#--------------------------------------------------------------------------
#           dye release                                                                                                                            
#---------------------------------------------------------------------------

#            FLAG_21 = -DDYE_RELEASE  

#--------------------------------------------------------------------------
#        SUSPENDED SEDIMENT MODEL:   UNCOMMENT TO INCLUDE MODEL
#--------------------------------------------------------------------------

#            FLAG_22 =  -DSEDIMENT

#--------------------------------------------------------------------------
#       Kalman Filters     TWO STEPS INCLUDED FOR REDUCED KALMAN FILTERS DATA
#                          ASSIMILATION. FIRSTLY USERS SHOULD OPEN FLAG_22 IN
#                          ORDER TO LET PROGRAM PREPARE INPUT DATA. AND THEN,
#                          SECOND OPTION FLAG_12 NEED BY APPLIED FOR ASSIMILATION
#                          BE NOTED: ONLY OPEN ONE OPTION AT ONE TIME                  
#--------------------------------------------------------------------------

#             FLAG_23 = -DENKF_ASSIM
#             FLAG_23 = -DRRK_PRE
#             FLAG_23 = -DRRK_ASSIM


#             KFLIB   =  -L/opt/intel/mkl/8.0.1/lib/em64t -lmkl_lapack -lmkl_em64t  -lguide -llapack -lblas 
#             KFLIB   = $(KFLIB) -L/hosts/salmon01/data00/medm/lib -lnetcdf


#--------------------------------------------------------------------------
#         Run One-D Mode with Biological Model
#--------------------------------------------------------------------------

#             FLAG_24 = -DONE_D_MODEL

#--------------------------------------------------------------------------
#       GENERAL BIOLOGICAL MODEL:    UNCOMMENT TO INCLUDE MODEL
#-------------------------------------------------------------------------- 
#            FLAG_25 = -DBioGen
#            BIOLIB       = -L../BIO_source -lbiomodule
#            BIOINCS      = -I../BIO_source/

#--------------------------------------------------------------------------
#         Dynamic/Thermodynamic Ice                                                                                                                           
#--------------------------------------------------------------------------

# NOTE: Must use -DSPHERICAL  and -DHEAT_FLUX
#            FLAG_26 = -DICE  

#--------------------------------------------------------------------------
#         CALCULATE THE NET HEAT FLUX IN MODEL                                                                                                                           
#--------------------------------------------------------------------------

# NOTE: Must use -DSPHERICAL  and -DICE
            FLAG_27 = -DHEAT_FLUX  

#--------------------------------------------------------------------------
# Visit online simulation mode
#--------------------------------------------------------------------------  

#            FLAG_28 = -DVISIT
#            FLAG_29 = -DVISIT_DEBUG

# GROUP INSTALL VERSION on em64t (Guppy)
#      VISITLIBPATH = -L/hosts/salmon01/data00/medm/visit/1.5.5/linux-x86_64/libsim/V1/lib -L/hosts/salmon01/data00/medm/visit/1.5.5/linux-x86_64/lib  -L/usr/X11R6/lib64
#      VISITLIB     = -lm -ldl -lsimf -lsim -lpthread
#      VISITINC     = -I/hosts/salmon01/data00/medm/visit/1.5.5/linux-x86_64/libsim/V1/include

# GROUP INSTALL VERSION on i386 (Salmon and Minke)
#      VISITLIBPATH = -L/hosts/salmon01/data00/medm/visit/1.5.5/linux/libsim/V1/lib -L/hosts/salmon01/data00/medm/visit/1.5.5/linux/lib -L/usr/X11R6/lib
#      VISITLIB     = -lm -ldl -lsimf -lsim -lpthread
#      VISITINC     = -I/hosts/salmon01/data00/medm/visit/1.5.5/linux/libsim/V1/include

# GROUP INSTALL VERSION on Hydra and Typhoeus
#      VISITLIBPATH = -L/usr/local/src/visit/1.5.5/linux-x86_64/libsim/V1/lib -L/usr/local/src/visit/1.5.5/linux-x86_64/lib -L/usr/X11R6/lib64
#      VISITLIB     =  -lm -ldl -lsimf -lsim -lpthread
#      VISITINC     = -I/usr/local/src/visit/1.5.5/linux-x86_64/libsim/V1/include




# USE DEVELOPER INSTALL VISIT
#      VISITLIBPATH = -L/hosts/salmon01/data00/medm/visit/src_64/visit/archives -L/usr/X11R6/lib64 -L/hosts/salmon01/data00/medm/visit/src_64/visit/lib
#      VISITLIB     =  -lm -ldl -lsimf -lsim -lpthread
#      VISITINC     = -I/hosts/salmon01/data00/medm/visit/src_64/visit/include/visit

#--------------------------------------------------------------------------
#  Semi-Implicit time stepping method
#--------------------------------------------------------------------------

 #           FLAG_30 = -DSEMI_IMPLICIT

#	    PETSC_LIB = -L/data/users/bclark/software/petsc/lib -lpetsc
#	    PETSC_FC_INCLUDES= -I/data/users/bclark/software/petsc/include
#--------------------------------------------------------------------------
#       FOR UPWIND LEAST SQUARE SCHEME:
#           IF LIMITED_NO IS ACTIVATED, NO LIMITED USED
#           IF LIMITED_1  IS ACTIVATED, FIRST ORDER LIMITED SCHEME USED
#           IF LIMITED_2  IS ACTIVATED, SECOND ORDER LIMITED SCHEME USED
#   !!!!!! ONLY ONE OF THE FLAGS BELOW CAN BE AND MUST BE CHOSEN
#--------------------------------------------------------------------------

#            FLAG_31 = -DLIMITED_NO
            FLAG_31 = -DLIMITED_1
#            FLAG_31 = -DLIMITED_2

#--------------------------------------------------------------------------
#     CELL/EDGE VELOCITY BLOCK MODULE BY T.W. (NEW VERSION APRIL 2011)
#     THIS NEW MODULE IS REFERENCED TO MOD_WD.F
#     REMOVE # TO INCLUDE THE MODULE
#--------------------------------------------------------------------------

#            FLAG_32 = -DV_BLOCK


#--------------------------------------------------------------------------
#
# NOTE: To activate Kelp Module, by T.W., for Uranium extraction project, April 2013
# lastly edited in June 2013
#--------------------------------------------------------------------------
#            FLAG_33 = -DKELP

#--------------------------------------------------------------------------
#
# NOTE: To activate OUTPUT_FLUX Module, by Wen Long
#--------------------------------------------------------------------------
            FLAG_34 = -DOUTPUT_FLUX

#--------------------------------------------------------------------------
#        SELECT COMPILER/PLATFORM SPECIFIC DEFINITIONS
#          SELECT FROM THE FOLLOWING PLATFORMS OR USE "OTHER" TO DEFINE
#          THE FOLLOWING VARIABLES:
#          CPP:  PATH TO C PREPROCESSOR 
#           FC:  PATH TO FORTRAN COMPILER (OR MPI COMPILE SCRIPT)
#          OPT:  COMPILER OPTIONS
#       MPILIB:  PATH TO MPI LIBRARIES (IF NOT LINKED THROUGH COMPILE SCRIPT)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#        COMPAQ/ALPHA Definitions
#--------------------------------------------------------------------------
#        CPPFLAGS   = $(DEF_FLAGS)  -DCOMPAQ
#        CPP        = /bin/cpp
#        FC         = f90
#        DEBFLGS    =  # -check bounds -check overflow -g 
#        OPT        =  -fast -arch ev6 -fpe1
#--------------------------------------------------------------------------
#  CRAY Definitions         
#--------------------------------------------------------------------------
#         CPPFLAGS = $(DEF_FLAGS) -DCRAY 
#         CPP      = /opt/ctl/bin/cpp
#         FC       = f90
#         DEBFLGS  =  
#         OPT      =
#--------------------------------------------------------------------------
#   Linux/Portland Group Definitions 
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp
#         CPPFLAGS = $(DEF_FLAGS) 
#         FC       = pgf90
#         DEBFLGS  = -Mbounds -g -Mprof=func
#         OPT      = #-fast  -Mvect=assoc,cachesize:512000,sse  
#--------------------------------------------------------------------------
#  Intel Compiler Definitions
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp
#         CPPFLAGS = $(DEF_FLAGS) -DINTEL 
#         FC       = ifort
#         DEBFLGS  = #-check all
# Use 'OPT = -O0 -g' for fast compile to test the make 
# Use 'OPT = -xP' for fast run on em64t (Hydra, Typhoeus and Guppy)
# Use 'OPT = -xN' for fast run on em64t (Hydra, Typhoeus and Guppy)
#         OPT      = -O0 -g  
#         OPT      = -xP 

# Do not set static for use with visit!
#         VISOPT   = -Wl,--export-dynamic
#         LDFLAGS  = $(VISITLIBPATH)
#--------------------------------------------------------------------------
#  Intel/MPI Compiler Definitions (SMAST)      
#--------------------------------------------------------------------------
#          CPP      = /usr/bin/cpp
#          CPPFLAGS = $(DEF_FLAGS) -DINTEL   
#          FC       = mpif90 
#         DEBFLGS  = #-check all
# Use 'OPT = -O0 -g'  for fast compile to test the make 
# Use 'OPT = -xP'      for fast run on em64t (Hydra, Typhoeus and Guppy)
# Use 'OPT = -xN'      for fast run on em64t (Salmon and Minke)
#         OPT      = -O0 -g 
#          OPT      = -xP
#          OPT      = -xN
	 # CPP =   /opt/intel/bin/fpp
          CPP      = cpp
          CPPFLAGS = $(DEF_FLAGS) -DINTEL -traditional
          FC       = mpiifort -heap-arrays  
#           FC      = ifort
          DEBFLGS  = #-check all 
          OPT      = -DUSE_U_INT_FOR_XDR -DHAVE_RPC_RPC_H=1  #-O3 -xN -axN -tpp7
          #OPT      = -O3 #-xN -axN -tpp7
          OPT = -traceback -heap-arrays -i4 -O3 -vec-report0 -r8 -xHost -no-prec-div  -no-prec-sqrt -assume nominus0 -assume noprotect_parens -assume norealloc_lhs  -shared-intelLIB    

# Do not set static for use with visit!
#         VISOPT   = -Wl,--export-dynamic
#         LDFLAGS  = $(VISITLIBPATH)
#--------------------------------------------------------------------------
#  gfortran defs 
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp 
#         CPPFLAGS = $(DEF_FLAGS)  -DGFORTRAN
#         FC       = gfortran  -O3 
#         DEBFLGS  = 
#         OPT      = 
#         CLIB     = 
#--------------------------------------------------------------------------
#  absoft / mac os x defs
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp 
#         CPPFLAGS = $(DEF_FLAGS)  -DABSOFT  
#         FC       = f90  -O3  -lU77
#         DEBFLGS  = 
#         OPT      = 
#         CLIB     = 
#--------------------------------------------------------------------------
#  IBM/AIX Definitions
#--------------------------------------------------------------------------
#         CPPFLAGS = $(DEF_FLAGS) -DAIX 
#         CPP      = /usr/local/bin/cpp
#         FC       = mpxlf90 -qsuffix=f=f90
#         DEBFLGS  =  # -qcheck -C -g
#         OPT      = -O -qarch=pwr4 -qtune=pwr4  -bmaxdata:0x80000000 -qhot -qmaxmem=8096
#--------------------------------------------------------------------------
#  APPLE OS X/XLF Definitions (G5)
#--------------------------------------------------------------------------
#         CPPFLAGS = $(DEF_FLAGS) -DAIX
#         CPP      = /usr/bin/cpp
#         FC       =  /opt/ibmcmp/xlf/8.1/bin/xlf90 -qsuffix=f=f90
#         DEBFLGS  =  # -qcheck -C -g
#         OPT      = -O5 -qarch=g5 -qtune=g5  -qhot -qmaxmem=8096 -qunroll=yes -Wl,-stack_size,10000000
#--------------------------------------------------------------------------
#  OTHER                   
#--------------------------------------------------------------------------
#         CPPFLAGS = $(DEF_FLAGS)
#         CPP      = 
#         FC       = 
#         DEBFLGS  = 
#         OPT      = 
#==========================================================================
#  END USER DEFINITION SECTION
#==========================================================================

         FFLAGS = $(DEBFLGS) $(OPT)
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib

#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
           CPPARGS =  $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) \
   $(FLAG_3) $(FLAG_4) $(FLAG_5) $(FLAG_6)  \
   $(FLAG_7) $(FLAG_8) $(FLAG_9) $(FLAG_10) \
   $(FLAG_11) $(FLAG_12) $(FLAG_13) $(FLAG_14) \
   $(FLAG_15) $(FLAG_16) $(FLAG_17) $(FLAG_18) \
   $(FLAG_19) $(FLAG_20) $(FLAG_21) $(FLAG_22) \
   $(FLAG_23) $(FLAG_24) $(FLAG_25) $(FLAG_26) \
   $(FLAG_27) $(FLAG_28) $(FLAG_29) $(FLAG_30) \
   $(FLAG_31) $(FLAG_32) $(FLAG_33) $(FLAG_34)
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------

            LIBS  =	$(PV3LIB) $(CLIB)  $(PARLIB) $(IOLIBS)	\
			$(MPILIB) $(GOTMLIB) $(KFLIB) $(BIOLIB)	\
			$(OILIB) $(VISITLIB) $(PETSC_LIB)

            INCS  = 	$(IOINCS) $(GOTMINCS) $(BIOINCS) $(VISITINC) $(PETSC_FC_INCLUDES)

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPARGS) $(INCS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90 
#	\rm $*.f90
.F90.o:
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR) $< > $*.f  
	$(FC) -c $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR) $*.f  
#	\rm $*.f90


#--------------------------------------------------------------------------
#  FVCOM Source Code.
#--------------------------------------------------------------------------

# Note: Change the order of modules at your own risk!

MODS  = mod_prec.F	mod_types.F	mod_main.F	mod_utils.F	\
	mod_par.F	mod_inp.F	mod_wd.F	\
	mod_probe.F	mod_clock.F	mod_spherical.F	mod_obcs.F	\
	mod_wqm.F	mod_gotm.F	mod_equitide.F	mod_atmotide.F	\
	mod_heatflux.F	mod_northpole.F	mod_balance_2d.F		\
	mod_tsobc.F	mod_meanflow.F	mod_obcs2.F	mod_obcs3.F	\
	mod_dye.F	mod_newinp.F	mod_tridiag.F	mod_scal.F	\
	mod_sed.F	particle.F	linklist.F	mod_enkf.F	\
	mod_etkf.F	mod_rrk.F	mod_rrkassim.F			\
	mod_optimal_interpolation.F	mod_assim.F  	                \
        mod_ncdio.F mod_ncdave.F mod_lag.F mod_visit.F	mod_petsc.F     \
	mod_semi_implicit.F	mod_block.F	mod_kelp.F              \
        mod_output_flux.F

MAIN  = utilities.F	adjust2d3d.F	advave_edge_gcn.F adv_q.F	\
	advection_edge_gcn.F adv_s.F	adv_t.F		adv_uv_edge_gcn.F\
	alloc_vars.F	archive.F	arcrst.F	atg.F		\
	fct_s.F		fct_t.F 	baropg.F	bcmap.F		\
	bcond_bfw.F	bcond_gcn.F	bcond_ts.F	bcs_force.F	\
	brough.F	brough_ud.F	bracket.F	cell_area.F	\
	conv_over.F	data_run.F	dens.F		dens2.F		\
	depth_grad.F	domdec.F	extel_edge.F	extuv_edge.F	\
	genmap.F	getdim.F	initial_qql.F	hot_start_data.F\
	initial_ts.F	initial_uvel.F	iofiles.F	out_avge.F	\
	out_binary.F	out_sms_one.F	pdomdec.F	phy_baropg.F	\
	report.F	sectinf.F	set_bndry.F	set_sigma.F	\
	shape_coef_gcn.F sinter.F	startup.F	svan.F		\
	theta.F		us_fvcom.F	tge.F		vdif_q.F	\
	vdif_ts.F	vdif_uv.F	vertvl_edge.F	water_depth.F	\
	wreal.F		ghostuv.F	dens3.F	  	init_sed.F	\
	advave_edge_gcy.F	advection_edge_gcy.F	closefiles.F	\
	adv_uv_edge_gcy.F	bcond_gcy.F		shape_coef_gcy.F\
	depth_check.F	vdif_ts_gom.F	adjust_ts.F	shutdown_check.F\
	extelpf_edge.F	fct_q2.F	fct_q2l.F	viscofh.F	\
	print_vals.F	visitsim.F	rho_mean.F	adcor.F  	


MYSUBS = viz.F

BIOGEN	= mod_bio_3D.F	mod_onedtide.F

 SRCS = $(MODS) $(ICE) $(MAIN) $(MYSUBS) $(BIOGEN)

# object list is composed of free (*.F) and fixed (*.F90) source
 OBJS = $(SRCS:.F=.o) $(FIXED:.F90=.o)


#--------------------------------------------------------------------------
#  Linking Directives
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(VISOPT) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXEC)


#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		makedepf90  $(SRCS) >> makedepends


#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
	tar cvf fvcom.tar *.F *.F90  makefile exa_run.dat makedepends RELEASE_NOTES 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod

clobber:	clean
		/bin/rm -f *.f90 *.o fvcom

#--------------------------------------------------------------------------
#  Common rules for all Makefiles - do not edit.
#--------------------------------------------------------------------------

emptyrule::

#--------------------------------------------------------------------------
#  Empty rules for directories that do not have SUBDIRS - do not edit.
#--------------------------------------------------------------------------

install::
	@echo "install in $(CURRENT_DIR) done"

install.man::
	@echo "install.man in $(CURRENT_DIR) done"

Makefiles::

includes::
include ./makedepends
