#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := ../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # Generic stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
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

THIS = GEOS_LdasShared
LIB  = lib$(THIS).a
esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.mod *.x

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.mod *.x

esma_doc doc:
	@echo nothing to do here


#                  --------------------
#                  User Defined Targets
#                  --------------------


SRC  =  enkf_types.F90 catch_types.F90 LDAS_ensdrv_Globals.F90 LDAS_DateTimeMod.F90 LDAS_DriverTypes.F90 \
	LDAS_Convert.F90 LDAS_Exceptions.F90 LDAS_TileCoordType.F90 LDAS_PertTypes.F90 \
	LDAS_ensdrv_functions.F90 my_matrix_functions.F90 \
	LDAS_EASE_conv.F90 LDAS_TileCoordRoutines.F90 \
	LDAS_RepairForcing.F90  \
	LDAS_ensdrv_init_routines.F90 LDAS_ensdrv_mpi.F90 
SRC1 = nr_indexx.f my_lu_decomp.f
OBJS = $(SRC1:.f=.o) $(SRC:.F90=.o)
DEPS = $(OBJS:.o=.d) 
#OBJS = $(addsuffix .o,$(basename $(SRCS)))
#DEPS = $(addsuffix .d, $(basename $(SRCS)))

INC_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF)
MOD_DIRS =   $(INC_DIRS)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS  = $(foreach dir,$(INC_DIRS),$(M)$(dir)) 

vpath % $(INC_DIRS) 

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.



  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

