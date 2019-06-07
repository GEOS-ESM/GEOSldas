#
# recursive makefile for ESMA.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif


# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help :
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"


THIS := GEOSlandassim_GridComp
LIB   = lib$(THIS).a

#ALLDIRS = Shared GEOSmwRTM_GridComp GEOSlandobs_GridComp #GEOSenkf_GridComp
#SUBDIRS = $(wildcard $(ALLDIRS) )

TARGETS = esma_install install esma_clean clean esma_doc doc\
          help esma_distclean distclean


export ESMADIR BASEDIR ARCH SITE

#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

$(TARGETS):
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)

local_esma_clean local_clean:
	$(RM) -f *~ *.[aox] *.[Mm][Oo][Dd]

local_esma_distclean local_distclean:
	$(RM) -f *~ *.[aoxd] *.[Mm][Oo][Dd]

local_esma_doc local_doc:
	@ "not yet"

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := io_hdf5.F90 enkf_general.F90 adapt_types.F90 mwRTM_types.F90 clsm_ensupd_glob_param.F90 \
        mwRTM_routines.F90 clsm_ensupd_upd_routines.F90 clsm_ensdrv_drv_routines.F90 \
	clsm_ensupd_read_obs.F90 catch_bias_types.F90 clsm_bias_routines.F90 clsm_adapt_routines.F90 \
	clsm_ensupd_enkf_update.F90 clsm_ensdrv_out_routines.F90 GEOS_LandAssimGridComp.F90
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)
INC_LANDPERT=$(ESMAINC)/GEOSlandpert_GridComp
INC_ENS=$(ESMAINC)/GEOSens_GridComp
INC_DIRS = . $(INC_HDF5) $(INC_GEOS_LDASSHR) $(INC_NETCDF) $(INC_GMAO_SHARED) \
	     $(INC_GEOS_LANDSHR) $(INC_GEOS_SURF) $(INC_ESMF) $(INC_LANDPERT) $(INC_ENS)
MOD_DIRS =   $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(PP) $(D)LDAS_MPI

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	     $(AR) $(AR_FLAGS) $(LIB) $(OBJS)


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif


#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------

#    GEOSlandobs_GridComp_install:  Shared_install
#    GEOSenkf_GridComp_install: Shared_install
#    GEOSlandassim_GridComp_install:  Shared_install

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
