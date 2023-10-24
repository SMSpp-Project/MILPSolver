##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of cpx_pars / scip_pars / grb_pars / highs_pars                 #
#                                                                            #
#                              Antonio Frangioni                             #
#                              Enrico Calandrini                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# basic directory
DIR = .

# define the directory of "core" SMS++ as required by the following includes
SMS++SDR = $(DIR)/../../SMS++

# OS-specific settings
# defines some general values, like UNAME_S, which identify the type of OS
# and architecture for which the project is being compiled, so that other
# makefiles can rely on these to automatically adapt without a need for the
# user to manually edit them
include $(SMS++SDR)/lib/makefile-OS

# default paths
# include a makefile that defines all the *_ROOT paths needed by every other
# makefile (possibly after automatic OS adaptation, see above) relatively to
# the places where all external libraries are to be found; use the previously
# defined UNAME_S to load the one of the current architecture
ifeq ($(UNAME_S),Linux)
    include $(SMS++SDR)/../extlib/makefile-default-paths-linux
endif
ifeq ($(UNAME_S),Darwin)
    include $(SMS++SDR)/../extlib/makefile-default-paths-macos
endif
ifeq ($(UNAME_S),Windows)
    include $(SMS++SDR)/../extlib/makefile-default-paths-win
endif

# non-default paths
# optionally (note the "-" in front) include a makefile that defines all the
# *_ROOT paths needed by every other makefile (possibly after automatic OS
# adaptation, see above) to cater for external libraries being located in
# non-standard locations.
-include $(SMS++SDR)/../extlib/makefile-paths

# debug switches
#SW = -g -glldb -fno-inline -std=c++17 -ferror-limit=1
# production switches
SW = -O3 -std=c++17 -DNDEBUG

# compiler
CC = clang++

# include the makefiles of the individual *_pars- - - - - - - - - - - - - - -
# if you are not interested in some of them, just comment them away here

# cpx_pars
include $(DIR)/makefile_cpx

# grb_pars
include $(DIR)/makefile_grb

# highs_pars
include $(DIR)/makefile_highs

# scip_pars
include $(DIR)/makefile_scip

# collate all (defined) executables
NAME = $(CPX_PARS) $(GRB_PARS) $(HIGHS_PARS) $(SCIP_PARS)

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(NAME)

############################ End of makefile #################################
