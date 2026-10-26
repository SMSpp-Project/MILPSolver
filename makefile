##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MILPSolver                                                   #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by MILPSolver, i.e., core SMS++ and those of all the    #
#   individual *MILPSolver (currently Cplex, Gurobi, SCIP, HiGHS, DAQP,      #
#   qpOASES, Clarabel, MOSEK, ProxQP, PIQP, OSQP, qpmad, quadprog and, on    #
#   Linux, PIPS-IPM++). These are *not* copied into $(MILPSINC): adding      #
#   those -I directives to the compile commands will have to done by         #
#   whatever "main" makefile is using this.                                  #
#   Analogously, any external library and the corresponding -L< libdirs >    #
#   will have to be added to the final linking command by whatever "main"    #
#   makefile is using this.                                                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by MILPSolver.                                                      #
#                                                                            #
#   Input:  $(CC)           = compiler command                               #
#           $(SW)           = compiler options                               #
#           $(SMS++INC)     = the -I$( core SMS++ directory )                #
#           $(SMS++OBJ)     = the core SMS++ library                         #
#           $(libCPLEXINC)  = the -I$( Cplex library )                       #
#           $(libGUROBIINC) = the -I$( Gurobi library )                      #
#           $(libSCIPINC)   = the -I$( SCIP library )                        #
#           $(libHiGHSINC)  = the -I$( HiGHS library )                       #
#           $(CPLEX_PRESENT), $(GUROBI_PRESENT), $(SCIP_PRESENT),            #
#           $(HiGHS_PRESENT) = non-empty if the solver is installed, each    #
#                              *MILPSolver being built only in that case     #
#                              [see extlib/makefile-lib*]                    #
#           $(libPIPSINC)   = the -I$( PIPS-IPM++ library ), Linux only      #
#           $(libDAQPINC)   = the -I$( DAQP library )                        #
#           $(libqpOASESINC) = the -I$( qpOASES library )                    #
#           $(libClarabelINC) = the -I$( Clarabel library )                  #
#           $(libMOSEKINC)  = the -I$( MOSEK library )                       #
#           $(MOSEK_PRESENT) = non-empty if MOSEK is installed               #
#           $(libProxSuiteINC) = the -I$( ProxSuite library )                #
#           $(libPIQPINC)   = the -I$( PIQP library )                        #
#           $(libOSQPINC)   = the -I$( OSQP library )                        #
#           $(libqpmadINC)  = the -I$( qpmad library )                       #
#           $(quadprog_ROOT) = the directory where quadprog is installed     #
#           $(MILPSSDR)     = the directory where the source is              #
#                                                                            #
#   Output: $(MILPSOBJ)     = the final object(s) / library                  #
#           $(MILPSH)       = the .h files to include                        #
#           $(MILPSINC)     = the -I$( source directory )                    #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# tools: the parameter headers of the back-ends - - - - - - - - - - - - - - -
# the tools of the source directory, which may be given as a relative or as
# an absolute path, or be empty, i.e., the current one
TOOLSSDR := $(if $(MILPSSDR),$(MILPSSDR),.)/tools
MILPSIDR := $(MILPSSDR)/include

# PIPS-IPM++ source tree, needed only by PIPSMILPSolver and tools/pips_pars;
# $(PIPSIPMPP_ROOT) is set by extlib/makefile-libPIPS, override it from the
# command line if your layout differs

# non-empty if and only if the PIPS-IPM++ source tree is available; all the
# PIPS-IPM++ additions below are only active in that case, so that builds
# without a PIPS-IPM++ checkout are unaffected
PIPS_PRESENT := $(wildcard \
	$(PIPSIPMPP_ROOT)/PIPS-IPM/Core/Interface/PIPSIPMppInterface.hpp)

# non-empty if and only if DAQP is installed in $(DAQP_ROOT), which is set by
# extlib/makefile-libDAQP; DAQPMILPSolver is only built in that case
DAQP_PRESENT := $(wildcard $(subst $\",,$(DAQP_ROOT))/include/daqp/api.h)

# non-empty if and only if qpOASES is installed in $(qpOASES_ROOT), which is
# set by extlib/makefile-libqpOASES; QPOASESMILPSolver is only built in that
# case
QPOASES_PRESENT := $(wildcard \
	$(subst $\",,$(qpOASES_ROOT))/include/qpOASES.hpp)

# non-empty if and only if Clarabel is installed in $(Clarabel_ROOT), which is
# set by extlib/makefile-libClarabel; ClarabelMILPSolver is only built in that
# case
CLARABEL_PRESENT := $(wildcard \
	$(subst $\",,$(Clarabel_ROOT))/include/c/DefaultSolver.h)

# non-empty if and only if ProxSuite is installed in $(proxsuite_ROOT), which
# is set by extlib/makefile-libProxSuite; ProxQPMILPSolver is only built in
# that case
ProxSuite_INCDIR := $(subst $\",,$(proxsuite_ROOT))/include
ProxSuite_PRESENT := $(wildcard $(ProxSuite_INCDIR)/proxsuite/config.hpp)

# non-empty if and only if PIQP is installed in $(piqp_ROOT), which is set by
# extlib/makefile-libPIQP; PIQPMILPSolver is only built in that case
PIQP_PRESENT := $(wildcard $(subst $\",,$(piqp_ROOT))/include/piqp/piqp.hpp)

# non-empty if and only if OSQP 1.0 is installed in $(osqp_ROOT), which is set
# by extlib/makefile-libOSQP; OSQPMILPSolver is only built in that case
OSQP_PRESENT := $(wildcard \
	$(subst $\",,$(osqp_ROOT))/include/osqp/osqp_api_functions.h)

# non-empty if and only if qpmad is installed in $(qpmad_ROOT), which is set
# by extlib/makefile-libqpmad; QPMADMILPSolver is only built in that case
QPMAD_PRESENT := $(wildcard $(subst $\",,$(qpmad_ROOT))/include/qpmad/solver.h)

# non-empty if and only if the library of quadprog is in $(quadprog_ROOT)/lib,
# which is set by extlib/makefile-libquadprog; QuadProgMILPSolver is only
# built in that case
QUADPROG_PRESENT := $(wildcard \
	$(subst $\",,$(quadprog_ROOT))/lib/libquadprog.*)

# each *MILPSolver.h includes a pair of files listing the algorithmic
# parameters of its back-end, produced by the corresponding tool in
# $(TOOLSSDR); the name of the pair carries the version of the back-end, hence
# it is the preprocessor that is asked for the very names the *MILPSolver.h
# will look for, so that upgrading a back-end regenerates the pair rather than
# leaving the stale one behind (see $(TOOLSSDR)/pars_names.probe)
PARSNAMES := $(shell $(CC) -E -P -x c++ $(libCPLEXINC) $(libGUROBIINC) \
	$(libSCIPINC) $(libHiGHSINC) $(if $(MOSEK_PRESENT),$(libMOSEKINC)) \
	"$(TOOLSSDR)/pars_names.probe" \
	2>/dev/null | sed -n 's/^[ \t]*SMSpp_pars //p' | tr -d ' \t')

# a back-end whose version macros do not expand is not usable anyway: drop it
# rather than generating a bogusly named pair
PARSNAMES := $(foreach n,$(PARSNAMES),$(if $(findstring VERSION,$(n)),,$(n)))

CPXPARSD   := $(addprefix $(MILPSIDR)/,$(filter CPX%,$(PARSNAMES)))
GRBPARSD   := $(addprefix $(MILPSIDR)/,$(filter GRB%,$(PARSNAMES)))
SCIPPARSD  := $(addprefix $(MILPSIDR)/,$(filter SCIP%,$(PARSNAMES)))
HiGHSPARSD := $(addprefix $(MILPSIDR)/,$(filter HiGHS%,$(PARSNAMES)))
MOSEKPARSD := $(addprefix $(MILPSIDR)/,$(filter MOSEK%,$(PARSNAMES)))

# the parameters of PIPS-IPM++ are read out of its sources rather than out of
# a library, hence they are not versioned and pips_pars needs to be told where
# those sources are
PIPSPARSD := $(if $(PIPS_PRESENT),$(MILPSIDR)/PIPS_defs.h)
PIPSPARSS := $(PIPSIPMPP_ROOT)/PIPS-IPM/Core/Options/PIPSIPMppOptions.C

CPXPARSM   := $(CPXPARSD:_defs.h=_maps.h)
GRBPARSM   := $(GRBPARSD:_defs.h=_maps.h)
SCIPPARSM  := $(SCIPPARSD:_defs.h=_maps.h)
HiGHSPARSM := $(HiGHSPARSD:_defs.h=_maps.h)
MOSEKPARSM := $(MOSEKPARSD:_defs.h=_maps.h)
PIPSPARSM  := $(PIPSPARSD:_defs.h=_maps.h)

PARSH = $(CPXPARSD) $(CPXPARSM) $(GRBPARSD) $(GRBPARSM) \
	$(SCIPPARSD) $(SCIPPARSM) $(HiGHSPARSD) $(HiGHSPARSM) \
	$(PIPSPARSD) $(PIPSPARSM) $(MOSEKPARSD) $(MOSEKPARSM)

.PHONY: tools
tools: $(PARSH)

# each tool writes both files of its pair in one run; the *_maps.h is
# therefore made to depend on the *_defs.h, which carries the recipe, so that
# the two are never generated concurrently under a parallel make. A tool reads
# the parameters out of the back-end it is compiled against, hence it is
# rebuilt from scratch here: a leftover executable would happily write the
# headers of the version that was installed when it was last built

ifneq ($(CPXPARSD),)
$(CPXPARSD):
	@echo "[MILPSolver] generating $(@F) and $(@F:_defs.h=_maps.h)"
	@rm -f "$(TOOLSSDR)/cpx_pars" "$(TOOLSSDR)/cpx_pars.o"
	@$(MAKE) -C "$(TOOLSSDR)" cpx_pars
	@cd "$(TOOLSSDR)" && ./cpx_pars

$(CPXPARSM): $(CPXPARSD)
	@test -f $@ || ( cd "$(TOOLSSDR)" && ./cpx_pars )
endif

ifneq ($(GRBPARSD),)
$(GRBPARSD):
	@echo "[MILPSolver] generating $(@F) and $(@F:_defs.h=_maps.h)"
	@rm -f "$(TOOLSSDR)/grb_pars" "$(TOOLSSDR)/grb_pars.o"
	@$(MAKE) -C "$(TOOLSSDR)" grb_pars
	@cd "$(TOOLSSDR)" && ./grb_pars

$(GRBPARSM): $(GRBPARSD)
	@test -f $@ || ( cd "$(TOOLSSDR)" && ./grb_pars )
endif

ifneq ($(SCIPPARSD),)
$(SCIPPARSD):
	@echo "[MILPSolver] generating $(@F) and $(@F:_defs.h=_maps.h)"
	@rm -f "$(TOOLSSDR)/scip_pars" "$(TOOLSSDR)/scip_pars.o"
	@$(MAKE) -C "$(TOOLSSDR)" scip_pars
	@cd "$(TOOLSSDR)" && ./scip_pars

$(SCIPPARSM): $(SCIPPARSD)
	@test -f $@ || ( cd "$(TOOLSSDR)" && ./scip_pars )
endif

ifneq ($(HiGHSPARSD),)
$(HiGHSPARSD):
	@echo "[MILPSolver] generating $(@F) and $(@F:_defs.h=_maps.h)"
	@rm -f "$(TOOLSSDR)/highs_pars" "$(TOOLSSDR)/highs_pars.o"
	@$(MAKE) -C "$(TOOLSSDR)" highs_pars
	@cd "$(TOOLSSDR)" && ./highs_pars

$(HiGHSPARSM): $(HiGHSPARSD)
	@test -f $@ || ( cd "$(TOOLSSDR)" && ./highs_pars )
endif

ifneq ($(MOSEKPARSD),)
$(MOSEKPARSD):
	@echo "[MILPSolver] generating $(@F) and $(@F:_defs.h=_maps.h)"
	@rm -f "$(TOOLSSDR)/mosek_pars" "$(TOOLSSDR)/mosek_pars.o"
	@$(MAKE) -C "$(TOOLSSDR)" mosek_pars
	@cd "$(TOOLSSDR)" && ./mosek_pars

$(MOSEKPARSM): $(MOSEKPARSD)
	@test -f $@ || ( cd "$(TOOLSSDR)" && ./mosek_pars )
endif

ifneq ($(PIPSPARSD),)
$(PIPSPARSD): $(PIPSPARSS)
	@echo "[MILPSolver] generating $(@F) and $(@F:_defs.h=_maps.h)"
	@rm -f "$(TOOLSSDR)/pips_pars" "$(TOOLSSDR)/pips_pars.o"
	@$(MAKE) -C "$(TOOLSSDR)" pips_pars
	@cd "$(TOOLSSDR)" && ./pips_pars -s "$(PIPSIPMPP_ROOT)" ../include

$(PIPSPARSM): $(PIPSPARSD)
	@test -f $@ || \
	( cd "$(TOOLSSDR)" && ./pips_pars -s "$(PIPSIPMPP_ROOT)" ../include )
endif

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MILPSOBJ = $(MILPSSDR)/obj/MILPSolver.o

MILPSINC = -I$(MILPSSDR)/include

MILPSH = $(MILPSSDR)/include/MILPSolver.h

ifneq ($(CPLEX_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/CPXMILPSolver.o
    MILPSH += $(MILPSSDR)/include/CPXMILPSolver.h
endif

ifneq ($(GUROBI_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/GRBMILPSolver.o
    MILPSH += $(MILPSSDR)/include/GRBMILPSolver.h
endif

ifneq ($(SCIP_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/SCIPMILPSolver.o
    MILPSH += $(MILPSSDR)/include/SCIPMILPSolver.h
endif

ifneq ($(HiGHS_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/HiGHSMILPSolver.o
    MILPSH += $(MILPSSDR)/include/HiGHSMILPSolver.h
endif

ifneq ($(PIPS_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/PIPSMILPSolver.o
    MILPSH += $(MILPSSDR)/include/PIPSMILPSolver.h
endif

ifneq ($(DAQP_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/DAQPMILPSolver.o
    MILPSH += $(MILPSSDR)/include/DAQPMILPSolver.h
endif

ifneq ($(QPOASES_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/QPOASESMILPSolver.o
    MILPSH += $(MILPSSDR)/include/QPOASESMILPSolver.h
endif

ifneq ($(CLARABEL_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/ClarabelMILPSolver.o
    MILPSH += $(MILPSSDR)/include/ClarabelMILPSolver.h
endif

ifneq ($(MOSEK_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/MOSEKMILPSolver.o
    MILPSH += $(MILPSSDR)/include/MOSEKMILPSolver.h
endif

ifneq ($(ProxSuite_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/ProxQPMILPSolver.o
    MILPSH += $(MILPSSDR)/include/ProxQPMILPSolver.h
endif

ifneq ($(PIQP_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/PIQPMILPSolver.o
    MILPSH += $(MILPSSDR)/include/PIQPMILPSolver.h
endif

ifneq ($(OSQP_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/OSQPMILPSolver.o
    MILPSH += $(MILPSSDR)/include/OSQPMILPSolver.h
endif

ifneq ($(QPMAD_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/QPMADMILPSolver.o
    MILPSH += $(MILPSSDR)/include/QPMADMILPSolver.h
endif

ifneq ($(QUADPROG_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/QuadProgMILPSolver.o
    MILPSH += $(MILPSSDR)/include/QuadProgMILPSolver.h
endif

# the C++20 module smspp.milpsolver - - - - - - - - - - - - - - - - - - - - -
# [see makefile-mod of SMS++]

# the -D that select the headers of the back-ends in the module
MILPSMODD = $(addprefix -DMILPSolver_MODULE_,$(if $(CPLEX_PRESENT),CPLEX) \
	$(if $(GUROBI_PRESENT),GUROBI) $(if $(SCIP_PRESENT),SCIP) \
	$(if $(HiGHS_PRESENT),HiGHS) \
	$(if $(PIPS_PRESENT),PIPS) $(if $(DAQP_PRESENT),DAQP) \
	$(if $(QPOASES_PRESENT),qpOASES) $(if $(CLARABEL_PRESENT),Clarabel) \
	$(if $(MOSEK_PRESENT),MOSEK) $(if $(ProxSuite_PRESENT),ProxQP) \
	$(if $(PIQP_PRESENT),PIQP) $(if $(OSQP_PRESENT),OSQP) \
	$(if $(QPMAD_PRESENT),qpmad) $(if $(QUADPROG_PRESENT),quadprog))

$(eval $(call SMSppMODrule,MILPS,$(MILPSSDR),MILPSolver,smspp.milpsolver,\
	SMS++,$(MILPSMODD)))

# the interface includes the *MILPSolver.h, hence the parameter headers
$(MILPSSDR)/obj/MILPSolver_module.o: $(PARSH)

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	$(MAKE) -C "$(TOOLSSDR)" clean
	rm -f $(MILPSOBJ) $(MILPSSDR)/*~

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean
	rm -f $(MILPSSDR)/include/CPX*_defs.h $(MILPSSDR)/include/CPX*_maps.h
	rm -f $(MILPSSDR)/include/GRB*_defs.h $(MILPSSDR)/include/GRB*_maps.h
	rm -f $(MILPSSDR)/include/SCIP*_defs.h $(MILPSSDR)/include/SCIP*_maps.h
	rm -f $(MILPSSDR)/include/HiGHS*_defs.h \
	$(MILPSSDR)/include/HiGHS*_maps.h
	rm -f $(MILPSSDR)/include/PIPS*_defs.h $(MILPSSDR)/include/PIPS*_maps.h
	rm -f $(MILPSSDR)/include/MOSEK*_defs.h \
	$(MILPSSDR)/include/MOSEK*_maps.h
	rm -f $(TOOLSSDR)/.headers.stamp

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: tools clean distclean

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# each *MILPSolver.o also depends on the parameter headers of its back-end,
# which are generated by the rules above

$(MILPSSDR)/obj/CPXMILPSolver.o:   $(CPXPARSD) $(CPXPARSM)
$(MILPSSDR)/obj/GRBMILPSolver.o:   $(GRBPARSD) $(GRBPARSM)
$(MILPSSDR)/obj/SCIPMILPSolver.o:  $(SCIPPARSD) $(SCIPPARSM)
$(MILPSSDR)/obj/HiGHSMILPSolver.o: $(HiGHSPARSD) $(HiGHSPARSM)

ifneq ($(PIPS_PRESENT),)
$(MILPSSDR)/obj/PIPSMILPSolver.o:  $(PIPSPARSD) $(PIPSPARSM)
endif

ifneq ($(MOSEK_PRESENT),)
$(MILPSSDR)/obj/MOSEKMILPSolver.o: $(MOSEKPARSD) $(MOSEKPARSM)
endif

$(MILPSSDR)/obj/MILPSolver.o: $(MILPSSDR)/src/MILPSolver.cpp \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/MILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(SW)

$(MILPSSDR)/obj/CPXMILPSolver.o: $(MILPSSDR)/src/CPXMILPSolver.cpp \
	$(MILPSSDR)/include/CPXMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/CPXMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libCPLEXINC) $(SW)

$(MILPSSDR)/obj/SCIPMILPSolver.o: $(MILPSSDR)/src/SCIPMILPSolver.cpp \
	$(MILPSSDR)/include/SCIPMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/SCIPMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libSCIPINC) $(SW)

$(MILPSSDR)/obj/GRBMILPSolver.o: $(MILPSSDR)/src/GRBMILPSolver.cpp \
	$(MILPSSDR)/include/GRBMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/GRBMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libGUROBIINC) $(SW)

$(MILPSSDR)/obj/HiGHSMILPSolver.o: $(MILPSSDR)/src/HiGHSMILPSolver.cpp \
	$(MILPSSDR)/include/HiGHSMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/HiGHSMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libHiGHSINC) $(SW)

ifneq ($(PIPS_PRESENT),)
$(MILPSSDR)/obj/PIPSMILPSolver.o: $(MILPSSDR)/src/PIPSMILPSolver.cpp \
	$(MILPSSDR)/include/PIPSMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/PIPSMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libPIPSINC) $(SW)
endif

ifneq ($(DAQP_PRESENT),)
$(MILPSSDR)/obj/DAQPMILPSolver.o: $(MILPSSDR)/src/DAQPMILPSolver.cpp \
	$(MILPSSDR)/include/DAQPMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/DAQPMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libDAQPINC) $(SW)
endif

ifneq ($(QPOASES_PRESENT),)
$(MILPSSDR)/obj/QPOASESMILPSolver.o: $(MILPSSDR)/src/QPOASESMILPSolver.cpp \
	$(MILPSSDR)/include/QPOASESMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/QPOASESMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libqpOASESINC) $(SW)
endif

ifneq ($(CLARABEL_PRESENT),)
$(MILPSSDR)/obj/ClarabelMILPSolver.o: $(MILPSSDR)/src/ClarabelMILPSolver.cpp \
	$(MILPSSDR)/include/ClarabelMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/ClarabelMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libClarabelINC) $(SW)
endif

ifneq ($(MOSEK_PRESENT),)
$(MILPSSDR)/obj/MOSEKMILPSolver.o: $(MILPSSDR)/src/MOSEKMILPSolver.cpp \
	$(MILPSSDR)/include/MOSEKMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/MOSEKMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libMOSEKINC) $(SW)
endif

ifneq ($(ProxSuite_PRESENT),)
$(MILPSSDR)/obj/ProxQPMILPSolver.o: $(MILPSSDR)/src/ProxQPMILPSolver.cpp \
	$(MILPSSDR)/include/ProxQPMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/ProxQPMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libProxSuiteINC) $(SW)
endif

ifneq ($(PIQP_PRESENT),)
$(MILPSSDR)/obj/PIQPMILPSolver.o: $(MILPSSDR)/src/PIQPMILPSolver.cpp \
	$(MILPSSDR)/include/PIQPMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/PIQPMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libPIQPINC) $(SW)
endif

ifneq ($(OSQP_PRESENT),)
$(MILPSSDR)/obj/OSQPMILPSolver.o: $(MILPSSDR)/src/OSQPMILPSolver.cpp \
	$(MILPSSDR)/include/OSQPMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/OSQPMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libOSQPINC) $(SW)
endif

ifneq ($(QPMAD_PRESENT),)
$(MILPSSDR)/obj/QPMADMILPSolver.o: $(MILPSSDR)/src/QPMADMILPSolver.cpp \
	$(MILPSSDR)/include/QPMADMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/QPMADMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libqpmadINC) $(SW)
endif

ifneq ($(QUADPROG_PRESENT),)
$(MILPSSDR)/obj/QuadProgMILPSolver.o: $(MILPSSDR)/src/QuadProgMILPSolver.cpp \
	$(MILPSSDR)/include/QuadProgMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/QuadProgMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(SW)
endif

########################## End of makefile ###################################
