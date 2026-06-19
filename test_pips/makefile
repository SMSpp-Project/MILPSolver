##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test_pips                                                     #
#                                                                            #
##############################################################################

# module name
NAME = test_pips

# basic directory
DIR = .

# common stuff, following tests/TwoStageStochasticBlock
include ../../tests/makefile_common

include ../../extlib/makefile-default-paths-linux

# define input macros for SMS++ complete makefile, then include it
SMS++SDR = ../../SMS++
include $(SMS++SDR)/lib/makefile-c

# MILPSolver
MILPSSDR = ..
include $(MILPSSDR)/makefile-s

# UCBlock, needed to deserialize EC/UC nc4 instances through the Block factory
UCBckDIR = ../../UCBlock
include $(UCBckDIR)/lib/makefile

# main module

# object files
MOBJ = $(MILPSOBJ) $(UCBckOBJ) $(SMS++OBJ)

# libraries
MLIB = $(SMS++LIB) $(MILPSLIB) $(UCBckLIB)

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test_pips.o
	$(CC) -o $(DIR)/$(NAME) $+ $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h

# include directives
MINC = $(SMS++INC) $(MILPSINC) $(UCBckINC)

# includes
MH = $(SMS++H) $(MILPSH) $(UCBckH)

# compile command
$(DIR)/test_pips.o: $(DIR)/test_pips.cpp $(MH)
	$(CC) -c $(DIR)/test_pips.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
