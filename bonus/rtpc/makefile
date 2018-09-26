OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CXX = `root-config --cxx`
CFLAGS = `root-config --cflags` \
	-O3 -W -Wall -Wextra -Wno-long-long \
	-fno-common \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
# CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

RTPC_garf: RTPC.C
	$(CXX) $(CFLAGS) -c RTPC.C
	$(CXX) $(CFLAGS) -o RTPC_garf RTPC.o $(LDFLAGS)
	rm RTPC.o

