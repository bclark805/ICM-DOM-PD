#--------------------------------------------------------------------------
#  KURT GLAESEMANN did minor changes to turn on -DMULTIPROCESSOR
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/bash
            DEF_FLAGS     = -P  -C -traditional 
            EXEC          =./wqm_dom4_NCDF
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
#        MULTI_PROCESSOR    INCLUDES PARALLELIZATION WITH MPI
#                           REQUIRES LINKING MPI LIBRARIES OR COMPILING
#                           WITH A PRELINKED SCRIPT (mpif90/mpf90/etc)
#                           DEFAULT: NO PARALLEL CAPABILITY
#                           UNCOMMENT TO INCLUDE MPI PARALLEL CAPABILITY
#--------------------------------------------------------------------------

              FLAG_1 = -DMULTIPROCESSOR
#              PARLIB = -L/hosts/salmon01/data00/medm/lib -lmetis 
	       PARLIB = -L/data/users/bclark/chesfvm/METIS_source -lmetis


		MKLLIB = -L/opt/intel/mkl/lib/intel64
		MKLINC = -I/opt/intel/mkl/include/intel64/ilp64 -I/opt/intel/mkl/include
#--------------------------------------------------------------------------
#        NETCDF OUTPUT      DUMP OUTPUT INTO NETCDF FILES (yes/no)
#                           REQUIRES SYSTEM DEPENDENT NETCDF LIBRARIES
#                           COMPILED WITH SAME F90 COMPILER
#                           SET PATH TO LIBRARIES WITH IOLIBS      
#                           SET PATH TO INCLUDE FILES (netcdf.mod) WITH IOINCS
#                           DEFAULT: NO NETCDF OUTPUT
#                           UNCOMMENT TO INCLUDE NETCDF OUTPUT CAPABILITY
#--------------------------------------------------------------------------
           FLAG_2       =  -DNETCDF_IO
           IOLIBS       =  -L/usr/local/netcdf-3.6.1-intel/lib/  -lnetcdf
           IOINCS       =  -I/usr/local/netcdf-3.6.1-intel/include/

#-------------------------------------------------------------------------
#        Optional calculatson for the sediment and water column
#       if Sediment DOM is off, than the calls to the sediment DOM module are ignored
#       if photodegration is off, than there is no photodegradation of colored carbon
#       This is for the CHESFVM-CDOM model V1.1 
#-------------------------------------------------------------------------


           FLAG_3        = -DSEDIMENT_DOM   
 
            FLAG_4        = -DPHOTODEGRADATION
            FLAG_5        = -DBUDGET
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
#  Intel Compiler Definitions
#--------------------------------------------------------------------------
#        CPP      = /usr/bin/cpp
#        CPPFLAGS = $(DEF_FLAGS) -DINTEL 
#        FC       = /share/apps/intel/Compiler/11.0/069/bin/intel64/ifort      
#        DEBFLGS  = #-check all 
#        OPT      = #-O3 -xN -axN -tpp7
#        CLIB     = #-static-libcxa 
#--------------------------------------------------------------------------
#  Intel Compiler Definitions
#--------------------------------------------------------------------------
          CPP      = /usr/bin/cpp
          CPPFLAGS = $(DEF_FLAGS) -DINTEL -DDOUBLE_PRECISION 
          FC       =  mpiifort 
        #  DEBFLGS = -g
        # DEBFLGS  = -r8 -i4 -g  -O0  -check all -g -debug -ftrapuv -debug-parameters -fp-stack-check -fpe0 -traceback -automatic -warn all
   #       OPT      = -i4 -traceback -heap-arrays  -O2 -vec-report0 -r8 -xHost -no-prec-div  -no-prec-sqrt -assume nominus0 -assume noprotect_parens -assume norealloc_lhs  -shared-intel # -warn unused -check uninit #  -p -g  -fno-inline-functions
          OPT      = -DUSE_U_INT_FOR_XDR -DHAVE_RPC_RPC_H=1  #-O3 -xN -axN -tpp7
          #OPT      = -O3 #-xN -axN -tpp7
          OPT = -traceback -heap-arrays  -i4 -O3 -vec-report0 -r8 -xHost -no-prec-div  -no-prec-sqrt -assume nominus0 -assume noprotect_parens -assume norealloc_lhs  -shared-intelLIB 

  #       CLIB     = -static-libcxa 
#--------------------------------------------------------------------------
#   Linux/Portland Group Definitions 
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp
#         CPPFLAGS = $(DEF_FLAGS) 
#         FC       = pgf90
#         DEBFLGS  = -Mbounds -g -Mprof=func
#         OPT      = #-fast  -Mvect=assoc,cachesize:512000,sse  
#--------------------------------------------------------------------------
#  gfortran defs 
#--------------------------------------------------------------------------
#         CPP      = /usr/bin/cpp 
#         CPPFLAGS = $(DEF_FLAGS)  -DGFORTRAN
#         FC       = gfortran  -O3 
#         DEBFLGS  = 
#         OPT      = 
#         CLIB     = 
#==========================================================================
#  END USER DEFINITION SECTION
#==========================================================================

         FFLAGS = $(DEBFLGS) $(OPT) 
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib

#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
           CPPARGS = $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) $(FLAG_3) $(FLAG_4)  $(FLAG_5)
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------

            LIBS  = $(CLIB)  $(PARLIB) $(IOLIBS) $(MPILIB)  # -L/opt/intel/mkl/lib/intel64/libmkl_blas95_ilp64.a -L/opt/intel/mkl/lib/intel64/libmkl_lapack95_ilp64.a $(MKLLIB) -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
 #    $(MKLLIB)/libmkl_blas95_ilp64.a  $(MKLLIB)/lib/intel64/libmkl_lapack95_ilp64.a -liomp5 -lpthread -lm -ldl 
#            LIBS  = -L../SPARSKIT2 -lskit $(PARLIB)
#             LIBS  = $(PARLIB) $(PETSC_LIB)
#           INCS  = $(IOINCS)$(PETSC_FC_INCLUDES)$(IOINCS)
           INCS  = $(IOINCS) $(MKLINC)

#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPARGS) $(INCS) $*.F > $*.f90
	$(FC)  -c $(FFLAGS) $(INCS) $*.f90
#	\rm $*.f90

#  FVCOM Source Code.
#--------------------------------------------------------------------------

MAIN  = mod_prec.F	mod_types.F	mod_utils.F	mod_control.F	\
	mod_lims.F	mod_bcs.F	mod_sizes.F	mod_hydrovars.F	\
	mod_buffers.F 	mod_fileinfo.F	      mod_inp.F  \
	mod_wqm.F	mod_zoop.F	mod_wqminit.F\
	utilities.F	mod_par.F	mod_sed_sav_exchange_vars.F	\
	bracket.F	mod_sed_sf_exchange_vars.F	mod_sf.F	\
	mod_sed_df_exchange_vars.F	mod_df.F	   mod_wc_dom.F \
	mod_sed_dom_exchange.F          mod_sed_dom.F   mod_owq.F	 \
	mod_sav.F	mod_ba.F	mod_sed.F	mod_algal.F	\
	mod_bcmap.F	mod_obcs.F	mod_ncd.F	ncdio.F		\
	mod_tge.F	cell_area.F	pdomdec.F			\
	domdec.F	genmap.F	bcs_force.F	wqm_inputs.F	\
	mod_kin.F			tvds.F		\
	adv_wqm.F	vdif_wqm.F	vertvl.F	viscofh.F	\
	bcond_wqm.F	fct_nut.F	mod_filenames.F	dens2.F mod_ncdio_new_BC.F		\
       mod_hydro.F wqm_main.F  							\



 SRCS = $(MAIN)  

 OBJS = $(SRCS:.F=.o)

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)

#--------------------------------------------------------------------------
#  Target to create dependecies.
#--------------------------------------------------------------------------

depend:
		makedepf90  $(SRCS) > makedepends


#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
	tar cvf fvcom.tar *.F  makefile exa_run.dat makedepends RELEASE_NOTES 

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod *.f90

clobber:	clean
		/bin/rm -f *.f90 *.o wqm

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
# DO NOT DELETE
