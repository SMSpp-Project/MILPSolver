##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test_dynamic (for MILPSolver)                                #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = DynamicLP_test

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++17 -ferror-limit=1 -Wno-enum-compare

# debug switches
SW_DEBUG = -g3 $(COMMON_SW)
# debug switches with address sanitizer and extra pedantic warning
SW_DEBUG_ASAN = -g3 -glldb -fno-inline -fsanitize=undefined -fsanitize=address -fno-omit-frame-pointer -Wpedantic -Wextra -Wno-unused-parameter $(COMMON_SW)
# production switches with address sanitizer
SW_RELEASE_ASAN = -O3 -DNDEBUG -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined $(COMMON_SW)
# production switches
SW_RELEASE = -O3 -DNDEBUG $(COMMON_SW)

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(NAME)

# debug with address sanitizer target- - - - - - - - - - - - - - - - - - - -

debug_asan: SW = $(SW_DEBUG_ASAN)
debug_asan: $(NAME)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(NAME)

# release with address sanitizer target- - - - - - - - - - - - - - - - - - -

release_asan: SW = $(SW_RELEASE_ASAN)
release_asan: $(NAME)

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(NAME)

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

# MILPSolver
MILPSSDR = ..
include $(MILPSSDR)/makefile-s

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(SMS++OBJ) $(MILPSOBJ)

# libraries
MLIB =  $(SMS++LIB) $(MILPSLIB)

$(NAME): $(DIR)/test.o $(MOBJ)
	$(CC) -o $(DIR)/$(NAME) $^ $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(SMS++INC) $(MILPSINC)

# includes
MH =    $(SMS++H) $(MILPSH)

# compile command

$(DIR)/test.o: $(DIR)/test.cpp $(MH)
	$(CC) -c $(DIR)/test.cpp -o $@ $(MINC) $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug debug_asan release release_asan clean distclean

############################ End of makefile #################################
