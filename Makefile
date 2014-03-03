# define the C compiler to use
CC = gcc

# define any compile-time flags
CFLAGS = -Wall -g

# define any directories containing header files other than /usr/include
#
INCLUDES = -I./include

# define the C source files
APIS = ./src/rngs.c ./src/words.c ./src/rvgs.c ./src/math_api.c ./src/dataMatrix.c
MAIN = ./src/GS2A.c 

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
API_OBJS = $(APIS:.c=.o)
MAIN_OBJS = $(MAIN:.c=.o)

# define the executable file 
MAIN_APP = ./bin/GS2A

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

all:    $(MAIN_APP) 

$(MAIN_APP): $(API_OBJS) $(MAIN_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN_APP) $(API_OBJS) $(MAIN_OBJS) -lm 

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) $(API_OBJS) $(MAIN_OBJS) 

depend: $(SRCS)
	makedepend $(INCLUDES) $^
