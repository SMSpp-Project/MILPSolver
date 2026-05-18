##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MILPSolver                                                   #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by MILPSolver, i.e., core SMS++ and those of all the    #
#   individual *MILPSolver (currently Cplex, Gurobi, SCIP, HiGHS). These are #
#   *not* copied into $(MILPSINC): adding those -I directives to the compile #
#   commands will have to done by whatever "main" makefile is using this.    #
#   Analogously, any external library and the corresponding -L< libdirs >    #
#   will have to be added to the final linking command by  whatever "main"   #
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

# tools: build + run headers generators - - - - - - - - - - - - - - - - - - -
TOOLSSDR := ./$(MILPSSDR)/tools
STAMP := $(TOOLSSDR)/.headers.stamp

.PHONY: tools
tools: $(STAMP)

$(STAMP):
	@echo "[MILPSolver] building and running *_pars in $(TOOLSSDR)"
	@set -e; \
	$(MAKE) -C "$(TOOLSSDR)"; \
	cd "$(TOOLSSDR)"; \
	for gen in $$(find . -maxdepth 1 -type f -perm -111 -name "*_pars"); do \
	  base=$${gen#./}; \
	  case "$$base" in \
	    cpx_pars)   p=CPX ;; \
	    grb_pars)   p=GRB ;; \
	    scip_pars)  p=SCIP ;; \
	    highs_pars) p=HiGHS ;; \
	    *)          p= ;; \
	  esac; \
	  if [ -n "$$p" ] && \
	     find ../include -maxdepth 1 -name "$${p}*_defs.h" -print -quit | grep -q . && \
	     find ../include -maxdepth 1 -name "$${p}*_maps.h" -print -quit | grep -q . ; then \
	    echo " -> $$base (skip: headers already present)"; \
	  else \
	    echo " -> $$base"; \
	    "./$$base"; \
	  fi; \
	done
	@echo "[MILPSolver] headers ready"
	@touch "$(STAMP)"

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MILPSOBJ = $(MILPSSDR)/obj/MILPSolver.o \
	$(MILPSSDR)/obj/CPXMILPSolver.o \
	$(MILPSSDR)/obj/GRBMILPSolver.o \
	$(MILPSSDR)/obj/SCIPMILPSolver.o \
	$(MILPSSDR)/obj/HiGHSMILPSolver.o \
	$(MILPSSDR)/obj/PIPSMILPSolver.o

MILPSINC = -I$(MILPSSDR)/include

MILPSH = $(MILPSSDR)/include/MILPSolver.h \
	$(MILPSSDR)/include/CPXMILPSolver.h \
	$(MILPSSDR)/include/GRBMILPSolver.h \
	$(MILPSSDR)/include/SCIPMILPSolver.h \
	$(MILPSSDR)/include/HiGHSMILPSolver.h \
	$(MILPSSDR)/include/PIPSMILPSolver.h

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
	rm -f $(STAMP)

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: tools clean distclean

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MILPSSDR)/obj/MILPSolver.o:      | $(STAMP)
$(MILPSSDR)/obj/CPXMILPSolver.o:   | $(STAMP)
$(MILPSSDR)/obj/SCIPMILPSolver.o:  | $(STAMP)
$(MILPSSDR)/obj/GRBMILPSolver.o:   | $(STAMP)
$(MILPSSDR)/obj/HiGHSMILPSolver.o: | $(STAMP)
$(MILPSSDR)/obj/PIPSMILPSolver.o:  | $(STAMP)

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

$(MILPSSDR)/obj/PIPSMILPSolver.o: $(MILPSSDR)/src/PIPSMILPSolver.cpp \
	$(MILPSSDR)/include/PIPSMILPSolver.h \
	$(MILPSSDR)/include/MILPSolver.h $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)/src/PIPSMILPSolver.cpp -o $@ \
	-I$(MILPSSDR)/include $(SMS++INC) $(libPIPSINC) $(SW)

########################## End of makefile ###################################
