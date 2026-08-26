##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MILPSolver                                                   #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by MILPSolver, i.e., core SMS++ and those of all the    #
#   individual *MILPSolver (currently Cplex, Gurobi, SCIP, HiGHS and, on     #
#   Linux, PIPS-IPM++). These are *not* copied into $(MILPSINC): adding      #
#   those -I directives to the compile commands will have to done by         #
#   whatever "main" makefile is using this. Analogously, any external        #
#   library and the corresponding -L< libdirs > will have to be added to     #
#   the final linking command by whatever "main" makefile is using this.     #
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
#           $(libPIPSINC)   = the -I$( PIPS-IPM++ library ), Linux only      #
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
TOOLSSDR := ./$(MILPSSDR)/tools
MILPSIDR := $(MILPSSDR)/include

# PIPS-IPM++ source tree, needed only by PIPSMILPSolver and tools/pips_pars;
# $(PIPSIPMPP_ROOT) is set by extlib/makefile-libPIPS, override it from the
# command line if your layout differs

# non-empty if and only if the PIPS-IPM++ source tree is available; all the
# PIPS-IPM++ additions below are only active in that case, so that builds
# without a PIPS-IPM++ checkout are unaffected
PIPS_PRESENT := $(wildcard \
	$(PIPSIPMPP_ROOT)/PIPS-IPM/Core/Interface/PIPSIPMppInterface.hpp)

# each *MILPSolver.h includes a pair of files listing the algorithmic
# parameters of its back-end, produced by the corresponding tool in
# $(TOOLSSDR); the name of the pair carries the version of the back-end, hence
# it is the preprocessor that is asked for the very names the *MILPSolver.h
# will look for, so that upgrading a back-end regenerates the pair rather than
# leaving the stale one behind (see $(TOOLSSDR)/pars_names.probe)
PARSNAMES := $(shell $(CC) -E -P -x c++ $(libCPLEXINC) $(libGUROBIINC) \
	$(libSCIPINC) $(libHiGHSINC) "$(TOOLSSDR)/pars_names.probe" \
	2>/dev/null | sed -n 's/^[ \t]*SMSpp_pars //p' | tr -d ' \t')

# a back-end whose version macros do not expand is not usable anyway: drop it
# rather than generating a bogusly named pair
PARSNAMES := $(foreach n,$(PARSNAMES),$(if $(findstring VERSION,$(n)),,$(n)))

CPXPARSD   := $(addprefix $(MILPSIDR)/,$(filter CPX%,$(PARSNAMES)))
GRBPARSD   := $(addprefix $(MILPSIDR)/,$(filter GRB%,$(PARSNAMES)))
SCIPPARSD  := $(addprefix $(MILPSIDR)/,$(filter SCIP%,$(PARSNAMES)))
HiGHSPARSD := $(addprefix $(MILPSIDR)/,$(filter HiGHS%,$(PARSNAMES)))

# the parameters of PIPS-IPM++ are read out of its sources rather than out of
# a library, hence they are not versioned and pips_pars needs to be told where
# those sources are
PIPSPARSD := $(if $(PIPS_PRESENT),$(MILPSIDR)/PIPS_defs.h)
PIPSPARSS := $(PIPSIPMPP_ROOT)/PIPS-IPM/Core/Options/PIPSIPMppOptions.C

CPXPARSM   := $(CPXPARSD:_defs.h=_maps.h)
GRBPARSM   := $(GRBPARSD:_defs.h=_maps.h)
SCIPPARSM  := $(SCIPPARSD:_defs.h=_maps.h)
HiGHSPARSM := $(HiGHSPARSD:_defs.h=_maps.h)
PIPSPARSM  := $(PIPSPARSD:_defs.h=_maps.h)

PARSH = $(CPXPARSD) $(CPXPARSM) $(GRBPARSD) $(GRBPARSM) \
	$(SCIPPARSD) $(SCIPPARSM) $(HiGHSPARSD) $(HiGHSPARSM) \
	$(PIPSPARSD) $(PIPSPARSM)

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

MILPSOBJ = $(MILPSSDR)/obj/MILPSolver.o \
	$(MILPSSDR)/obj/CPXMILPSolver.o \
	$(MILPSSDR)/obj/GRBMILPSolver.o \
	$(MILPSSDR)/obj/SCIPMILPSolver.o \
	$(MILPSSDR)/obj/HiGHSMILPSolver.o

MILPSINC = -I$(MILPSSDR)/include

MILPSH = $(MILPSSDR)/include/MILPSolver.h \
	$(MILPSSDR)/include/CPXMILPSolver.h \
	$(MILPSSDR)/include/GRBMILPSolver.h \
	$(MILPSSDR)/include/SCIPMILPSolver.h \
	$(MILPSSDR)/include/HiGHSMILPSolver.h

ifneq ($(PIPS_PRESENT),)
    MILPSOBJ += $(MILPSSDR)/obj/PIPSMILPSolver.o
    MILPSH += $(MILPSSDR)/include/PIPSMILPSolver.h
endif

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	$(MAKE) -C "$(TOOLSSDR)" clean
	rm -f $(MILPSOBJ) $(MILPSSDR)/*~

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean
	rm -f $(MILPSSDR)/include/CPX*_defs.h $(MILPSSDR)/include/CPX*_maps.h
	rm -f $(MILPSSDR)/include/GRB*_defs.h $(MILPSSDR)/include/GRB*_maps.h
	rm -f $(MILPSSDR)/include/SCIP*_defs.h $(MILPSSDR)/include/SCIP*_maps.h
	rm -f $(MILPSSDR)/include/HiGHS*_defs.h $(MILPSSDR)/include/HiGHS*_maps.h
	rm -f $(MILPSSDR)/include/PIPS*_defs.h $(MILPSSDR)/include/PIPS*_maps.h
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

########################## End of makefile ###################################
