##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of the MILPSolver tests                                         #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# the executables built here
NAMES = test_cuts test_dual test_dynamic

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++20 -ferror-limit=1

# debug switches
SW_DEBUG = -g -glldb -fno-inline $(COMMON_SW)
# debug switches with address sanitizer and extra pedantic warning
SW_DEBUG_ASAN = -g3 -glldb -fno-inline $(COMMON_SW) -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG -DCLANG_1200_0_32_27_PATCH

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(addprefix $(DIR)/,$(NAMES))

# debug_asan target - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug_asan: SW = $(SW_DEBUG_ASAN)
debug_asan: $(addprefix $(DIR)/,$(NAMES))

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(addprefix $(DIR)/,$(NAMES))

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(addprefix $(DIR)/,$(NAMES))

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# core SMS++
SMS++SDR = ../../SMS++
include $(SMS++SDR)/lib/makefile-c

# MILPSolver (without core SMS++)
MILPSSDR = ..
include $(MILPSSDR)/makefile-s

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(MILPSOBJ) $(SMS++OBJ)

# libraries
MLIB =  $(MILPSLIB) $(SMS++LIB)

# every test is built from its own .o linked against the same modules
$(DIR)/test_cuts: $(MOBJ) $(DIR)/test_cuts.o
	$(CC) -o $@ $(DIR)/test_cuts.o $(MOBJ) $(MLIB) $(SW)

$(DIR)/test_dual: $(MOBJ) $(DIR)/test_dual.o
	$(CC) -o $@ $(DIR)/test_dual.o $(MOBJ) $(MLIB) $(SW)

$(DIR)/test_dynamic: $(MOBJ) $(DIR)/test_dynamic.o
	$(CC) -o $@ $(DIR)/test_dynamic.o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(SMS++INC) $(MILPSINC)

# includes
MH =    $(SMS++H) $(MILPSH)

# compile command

$(DIR)/%.o: $(DIR)/%.cpp $(MH)
	$(CC) -c $< -o $@ $(MINC) $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug debug_asan release clean distclean

############################ End of makefile #################################
