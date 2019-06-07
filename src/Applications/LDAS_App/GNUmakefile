#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 2014/09/24  pchakrab  New version for ESMF-ized LDAS


# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := ../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

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
	@echo "        FREAL = $(FREAL)"


THIS = GEOSldas
BINS = $(THIS).x
PRE_LDAS= preprocess_ldas.x
BIN2NC4 = tile_bin2nc4.x
MWRTMBIN2NC4 = mwrtm_bin2nc4.x
NETCDFDIR=$(BASEDIR)/$(ARCH)

esma_install install: $(BINS) $(PRE_LDAS) $(BIN2NC4) $(MWRTMBIN2NC4)
	$(MKDIR) $(ESMABIN)
	$(CP) -p $(BINS) $(ESMABIN)
	$(CP) -p ldas_setup $(ESMABIN)
	$(CP) -p LDASsa_*rc $(ESMAETC)
	$(CP) -p LDASsa_DEFAULT*nml $(ESMAETC)
	$(CP) -p lenkf.j.template $(ESMAETC)
	$(CP) -p $(PRE_LDAS) $(ESMABIN)
	$(CP) -p $(BIN2NC4) $(ESMABIN)
	$(CP) -p $(MWRTMBIN2NC4) $(ESMABIN)
	$(CP) -p process_hist.csh $(ESMABIN)
	$(CP) -p process_rst.csh $(ESMABIN)

esma_clean clean:
	$(RM) *~ *.[aox] *.mod *.x

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.mod *.x *.pyc

esma_doc doc:
	@echo No documentation yet


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = GEOSldas.F90
SRCS4 = preprocess_ldas.F90
SRCS5 = tile_bin2nc4.F90
SRCS6 = mwrtm_bin2nc4.F90

# Ensure that OBJS does not get deleted
.SECONDARY: $(OBJS)
OBJS = $(SRCS:.F90=.o)
DEPS = $(SRCS:.F90=.d)

LIB_SHARED = \
	$(ESMALIB)/libGEOS_Shared.a \
	$(ESMALIB)/libMAPL_Base.a \
        $(ESMALIB)/libMAPL_cfio_r4.a \
	$(ESMALIB)/libGMAO_gfio_r4.a \
	$(ESMALIB)/libGMAO_hermes.a \
	$(ESMALIB)/libGMAO_mpeu.a 

LIB_SURF = \
	$(wildcard \
	$(ESMALIB)/libGEOSsurface_GridComp.a \
        $(ESMALIB)/libGEOSlake_GridComp.a \
        $(ESMALIB)/libGEOSlandice_GridComp.a \
        $(ESMALIB)/libGEOSsaltwater_GridComp.a \
        $(ESMALIB)/libGEOSCICEThermo_GridComp.a \
        $(ESMALIB)/libGEOSland_GridComp.a \
        $(ESMALIB)/libGEOScatch_GridComp.a \
        $(ESMALIB)/libGEOScatchCN_GridComp.a \
        $(ESMALIB)/libGEOSvegdyn_GridComp.a \
	$(ESMALIB)/libGEOSmwRTM_GridComp.a \
	$(ESMALIB)/libGEOS_RadiationShared.a \
        $(ESMALIB)/libGEOS_SurfaceShared.a) \
	$(ESMALIB)/libGEOS_LandShared.a \

LIB_LANA = $(ESMALIB)/libGEOSlana_GridComp.a

LIB_LDAS = \
	$(ESMALIB)/libGEOSldas_GridComp.a \
	$(ESMALIB)/libGEOSmetforce_GridComp.a \
	$(ESMALIB)/libGEOSlandpert_GridComp.a \
	$(ESMALIB)/libGEOSens_GridComp.a \
	$(ESMALIB)/libGEOSlandassim_GridComp.a \
	$(ESMALIB)/libGEOS_LdasShared.a \
	$(ESMALIB)/libFVdycoreCubed_GridComp.a \
	$(ESMALIB)/libfvdycore.a \
	$(ESMALIB)/libGFDL_fms.a

LIB_DYN = \
	$(ESMALIB)/libFVdycoreCubed_GridComp.a \
	$(ESMALIB)/libfvdycore.a \
	$(ESMALIB)/libGFDL_fms.a 

LIB_MK = \
	$(ESMALIB)/libGEOS_LdasShared.a \
	$(ESMALIB)/libMAPL_Base.a \
	$(ESMALIB)/libMAPL_Base_stubs.a 
        
LIB_COMP = $(LIB_LDAS) $(LIB_SURF) $(LIB_SHARED)

FREAL = $(FREAL4)

INC_DIRS = . \
	$(INC_ESMF) $(INC_MAPL_BASE) $(NETCDFDIR)/include/netcdf \
	$(ESMAINC)/GEOSldas_GridComp \
	$(ESMAINC)/GEOS_LdasShared
	# $(ESMAINC)/GEOScatch_GridComp \

MOD_DIRS = $(INC_DIRS)

USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS = $(foreach dir,$(MOD_DIRS),$(M)$(dir))

export USER_FDEFS += -DDEBUG

%.x : $(OBJS) $(LIB_COMP)
	$(LD) -o $@ $(OBJS) \
	$(USER_FINCS) $(USER_FMODS) $(LIB_DYN) $(LIB_COMP) $(LIB_HDF5) $(LIB_DYN) $(LIB_ESMF) $(LIB_SDF) $(LIB_MPI) $(LIB_SCI) $(LIB_SYS)

$(PRE_LDAS) : $(SRCS4) $(LIB_MK)
	$(LD) -o $@ $(SRCS4) \
	$(USER_FINCS) $(USER_FMODS) $(LIB_MK) $(LIB_ESMF) $(LIB_SDF) $(LIB_SCI) $(LIB_MPI) $(LIB_SYS)


$(BIN2NC4) : $(SRCS5) $(LIB_MK)
	$(LD) -o $@ $(SRCS5) \
	$(USER_FINCS) $(USER_FMODS) $(LIB_MK) $(LIB_ESMF) $(LIB_SDF) $(LIB_SCI) $(LIB_MPI) $(LIB_SYS)

$(MWRTMBIN2NC4) : $(SRCS6) $(LIB_MK)
	$(LD) -o $@ $(SRCS6) \
	$(USER_FINCS) $(USER_FMODS) $(LIB_MK) $(LIB_ESMF) $(LIB_SDF) $(LIB_SCI) $(LIB_MPI) $(LIB_SYS)

#                  ---------------------
#                      Dependencies
#                  ---------------------


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
