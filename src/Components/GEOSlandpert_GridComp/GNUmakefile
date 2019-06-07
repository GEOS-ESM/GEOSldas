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


THIS := $(shell basename `pwd`)
LIB   = lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)

esma_clean clean:
	$(RM) -f *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	$(RM) -f *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo "not yet"


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := nr_ran2_gasdev.F90 nr_jacobi.F90 nr_fft.F90 \
	random_fields.F90 land_pert.F90 force_and_cat_progn_pert_types.F90 LDAS_PertRoutines.F90 \
	GEOS_LandPertGridComp.F90

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

INC_DIRS = . $(ESMAINC)/GEOScatch_GridComp $(ESMAINC)/GEOSmetforce_GridComp \
	     $(INC_GEOS_LANDSHR) \
	     $(INC_GEOS_LDASSHR) $(INC_MAPL_BASE) $(INC_NETCDF) $(INC_ESMF)
MOD_DIRS =   $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

export USER_FDEFS += -DMKL_AVAILABLE

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


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
