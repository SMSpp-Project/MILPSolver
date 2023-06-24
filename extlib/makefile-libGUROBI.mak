##############################################################################
############################ makefile-libGUROBI ##############################
##############################################################################
#                                                                            #
#   makefile of libGUROBI                                                    #
#                                                                            #
#   Input:  none, this is a pre-built library                                #
#                                                                            #
#   Output: accordingly, there is no *H and *OBJ in output, since there is   #
#           no need to check for changes in the .h and rebuild the .o / .a   #
#           $(libGUROBILIB) = external libreries + -L<libdirs> for libGUROBI #
#           $(libGUROBIINC) = the -I$(include directories) for libGUROBI     #
#                                                                            #
#                             Antonio Frangioni                              #
#                         Dipartimento di Informatica                        #
#                            Universita' di Pisa                             #
#                                                                            #
#                             Enrico Calandrini                              #
#                         Dipartimento di Matematica                         #
#                            Universita' di Pisa                             #
#                                                                            #
##############################################################################

# internal macros - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# basic Gurobi directory
libGUROBIBSCDIR = $(GUROBI_HOME)/

# lib Gurobi directory
libGUROBIINCDIR = $(libGUROBIBSCDIR)lib/

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

libGUROBILIB = -L$(libGUROBIINCDIR) $(libGUROBIINCDIR)libgurobi_g++5.2.a $(libGUROBIINCDIR)libgurobi_c++.a $(libGUROBIINCDIR)libgurobi100.so
libGUROBIINC = -I$(libGUROBIBSCDIR)include/

############################# End of makefile ################################