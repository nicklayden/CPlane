#compiler
cxx = g++

# executable name
program = plane

# source files
srcs = main.cpp plotter.cpp

# directory to sources and compiled obj files
srcdir = src
objdir = obj

objs = $(srcs:.cpp=.o)
bins = $(src:.cpp=)

# global cpp flags
cppflags = -std=c++11 -O3 -Wno-deprecated-declarations

# specialized cpp flags - per compile basis
cppflags += $(cppflags-$@)
# flags needed only to compile executable at final step
cppflags-plane += -lsfml-graphics -lsfml-system -lsfml-window 

# linker flags:
linkflags = -stdlib=libc++


# giving folder prefix to the executable files.
exeobjs = $(addprefix $(objdir)/, $(objs))

# recipe to make
all: $(program)

$(program): $(exeobjs)
	$(cxx) -o $(program) $(exeobjs) $(cppflags) 

$(objdir)/%.o: $(srcdir)/%.cpp | objdirmk
	$(cxx) $(cppflags) -c $< -o $@

# check to see if obj/ is already a directory. if not, make it.
objdirmk:
	@mkdir -p $(objdir)

.PHONY: clean
# removes the executable and the compiled obj .o files
clean:
	@echo Removing executable and cleaning object directories.
	$(RM) $(program) $(EXEOBJS) 
	$(RM) -r $(objdir)

