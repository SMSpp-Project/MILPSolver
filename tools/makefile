##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of cpx_pars / scip_pars / grb_pars / highs_pars                 #
#                                                                            #
#                              Antonio Frangioni                             #
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
    ifeq ($(UNAME_P),x86_64)
        include $(SMS++SDR)/../extlib/makefile-default-paths-macos-x86_64
    else
        include $(SMS++SDR)/../extlib/makefile-default-paths-macos-arm64
    endif
endif

# non-default paths
# optionally (note the "-" in front) include a makefile that defines all the
# *_ROOT paths needed by every other makefile (possibly after automatic OS
# adaptation, see above) to cater for external libraries being located in
# non-standard locations.
-include $(SMS++SDR)/../extlib/makefile-paths

# common flags
COMMON_SW = -std=c++17

# debug switches
SW_DEBUG = -g -glldb -fno-inline
# production switches
SW_RELEASE = -O3 -DNDEBUG

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(COMMON_SW) $(SW_DEBUG)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(COMMON_SW) $(SW_RELEASE)

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

#default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#debug: SW = $(SW_DEBUG)
debug: $(NAME)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#release: SW = $(SW_RELEASE)
release: $(NAME)

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(NAME)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug release clean distclean

############################ End of makefile #################################
