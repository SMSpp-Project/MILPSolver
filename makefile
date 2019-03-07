##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MILPSolver                                                   #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by MILPSolver, i.e., core SMS++ and Cplex. These are    #
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
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(libCPLEXINC) = the -I$( Cplex library )                        #
#           $(MILPSSDR)  = the directory where the source is                 #
#                                                                            #
#   Output: $(MILPSOBJ) = the final object(s) / library                      #
#           $(MILPSH)   = the .h files to include                            #
#           $(MILPSINC) = the -I$( source directory )                        #
#                                                                            #
#                                VERSION 1.00                                #
#                               07 - 03 - 2019                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Operations Research Group                         #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################


# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MILPSOBJ = $(MILPSSDR)obj/MILPSolver.o

MILPSINC = -I$(MILPSSDR)include/

MILPSH   = $(MILPSSDR)include/MILPSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MILPSOBJ) $(MILPSSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MILPSOBJ): $(MILPSSDR)src/MILPSolver.cpp $(MILPSH) $(SMS++OBJ)
	$(CC) -c $(MILPSSDR)src/MILPSolver.cpp -o $@ \
	$(MILPSINC) $(SMS++INC) $(libCPLEXINC) $(SW)

########################## End of makefile ###################################
