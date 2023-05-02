# the compiler: gcc for C program, define as g++ for C++
CC = gcc

#compiler flags:
# 	-g		adds debugging information to the executable file with gdb (GNU DeBugging)
# 	-Wall	turns on most, but not all, compiler warnings
CFLAGS = -g -Wall -Wextra

# defines any directories containing header files other than /usr/include
INCLUDES = -I/usr/local/include -I../include

# define library paths in addition to /usr/lib
LFLAGS = -L/usr/local/include -L../lib

# define any libraries to link into executable
#		-lgsl			links GSL executables
#		-lgslcblas	links GSL Cblas executabes for Cblas operations
#		-lm 			links system math
LIBS = -lgsl -lgslcblas -lm

# define the source path
SRCDIR = ../src/
# define the obj path
OBJDIR = ../obj/
# define the binary path
BINDIR = ../bin/

VPATH = $(SRCDIR) $(OBJDIR)

# define the C source and object files
FILES = main array_utils matrix_utils matrix_calc_utils graph_alloc_utils graph_utils extern_var opt_utils random_matrix_utils model model_graphs model_mc model_min
SRCS  = $(addsuffix .c, $(FILES))
OBJS = $(addsuffix .o, $(FILES))

# define the executable file
MAIN = main.exe

.PHONY: depend clean

all: $(MAIN)
	@echo SUCCESS: $(MAIN) has been compiled
	
$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(BINDIR)$(MAIN) $(addprefix $(OBJDIR), $(OBJS)) $(LFLAGS) $(LIBS)
	
# this is a suffix replacement rule for building ,o's from .c's 
# it uses automatic variables $<: the name of the prerequisite of
# the rule (a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $^ -o $(OBJDIR)$@
	
clean:
	$(RM) *.o $(OBJDIR)*.o *~ $(BINDIR)$(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^
	
# DO NOT DELETE THIS LINE -- make depend needs it