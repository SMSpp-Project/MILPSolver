##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test_cuts (for MILPSolver)                                   #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = cuts_test

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++17 -DCLANG_1200_0_32_27_PATCH -ferror-limit=1

# debug switches
SW_DEBUG = -g3 -glldb -fno-inline -Wno-enum-compare $(COMMON_SW)
# debug switches with address sanitizer and extra pedantic warning
SW_DEBUG_ASAN = -g3 -glldb -fno-inline -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter -Wno-enum-compare $(COMMON_SW)
# production switches with address sanitizer
SW_RELEASE_ASAN = -O3 -DNDEBUG -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -Wno-enum-compare $(COMMON_SW)
# production switches
SW_RELEASE = -O3 -DNDEBUG -Wno-enum-compare $(COMMON_SW)

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(DIR)/$(NAME)

# debug with address sanitizer target- - - - - - - - - - - - - - - - - - - -

debug_asan: SW = $(SW_DEBUG_ASAN)
debug_asan: $(DIR)/$(NAME)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(DIR)/$(NAME)

# release with address sanitizer target- - - - - - - - - - - - - - - - - - -

release_asan: SW = $(SW_RELEASE_ASAN)
release_asan: $(DIR)/$(NAME)

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# define input macros for SMS++ complete makefile, then include it
SMS++SDR = ../../SMS++
include $(SMS++SDR)/lib/makefile-c

# define input macros for MILPSolver makefile (minus SMS++), then include it
MILPSSDR = ..
include $(MILPSSDR)/makefile-s

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(MILPSOBJ) $(SMS++OBJ)

# libraries
MLIB =  $(SMS++LIB) $(MILPSLIB)

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test_cuts.o
	$(CC) -o $(DIR)/$(NAME) $^ $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(SMS++INC) $(MILPSINC)

# includes
MH =    $(SMS++H) $(MILPSH)

# compile command

$(DIR)/test_cuts.o: $(DIR)/test_cuts.cpp $(MH)
	$(CC) -c $(DIR)/test_cuts.cpp -o $@ $(MINC) $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug debug_asan release release_asan clean distclean

############################ End of makefile #################################
