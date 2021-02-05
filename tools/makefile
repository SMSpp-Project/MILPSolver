##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of cpx_pars                                                     #
#                                                                            #
#                                VERSION 1.00                                #
#                               25 - 10 - 2020                               #
#                                                                            #
#                              Implementation by:                            #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
NAME = cpx_pars

# basic directory
DIR = ./

# debug switches
#SW = -g -glldb -fno-inline -std=c++17 -ferror-limit=1 $(MCFCS)
# production switches
SW = -O3 -std=c++17 -DNDEBUG $(MCFCS)

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# define input macros for libCPLEX complete makefile, then include it
libCPLEXDIR = $(MILPSSDR)../../extlib/
include $(libCPLEXDIR)makefile-libCPLEX

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  

# libreries
MLIB =  $(libCPLEXLIB)

$(NAME): $(MOBJ) $(DIR)$(NAME).o
	$(CC) -o $(NAME) $(DIR)$(NAME).o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

# include directives
MINC =  $(libCPLEXINC)

# includes
MH =    

# compile command

$(DIR)$(NAME).o: $(DIR)$(NAME).cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
