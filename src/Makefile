-include Makefile.inc
#CXX=clang++
VPATH=./include
CXX=g++
OPTFLAGS=-g
VERBOSITY=NOVERBOSE
EIGEN_DIR=/usr/local/include/eigen3
PYBIND11_DIR=../extern/pybind11/include/
STDFLAGS=-std=c++20
CXXFLAGS+=$(STDFLAGS) $(OPTFLAGS) -Wall 

CPPFLAGS=-DDIMENSION=$(DIM) -I./include -I$(EIGEN_DIR)

DEPEND=make.dep

CC=$(CXX)
#
# The general setting is as follows:
# mains are identified bt main_XX.cpp
# all other files are XX.cpp
#

# get all files *.cpp
SRCS=$(wildcard *.cpp)
# get the corresponding object file
OBJS_ALL = $(SRCS:.cpp=.o)
# take out python interface
HEADERS=$(wildcard *.hpp)
#
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)
OBJS_NOEXEC=$(filter-out main%.o,$(OBJS))
OBJS=$(OBJS_ALL)

#========================== ORA LA DEFINIZIONE DEGLI OBIETTIVI
.phony= 2d 3d objs all alldim clean distclean doc

.DEFAULT_GOAL = all

all: $(EXEC)

objs: $(DEPEND) $(OBJS_NOEXEC)

clean:
	$(RM) -f $(EXEC) $(OBJS) 

distclean:
	$(MAKE) clean
	$(RM) -f ./doc $(DEPEND)
	$(RM) *.out *.bak *~
	find 2D -type f \! -name README.txt | xargs /bin/rm -f
	find 3D -type f \! -name README.txt | xargs /bin/rm -f
doc:
	doxygen $(DOXYFILE)

$(EXEC): $(OBJS)


$(DEPEND): $(SRCS)
	$(RM) $(DEPEND)
	for f in $(SRCS); do \
	$(CXX) $(STDFLAGS) $(CPPFLAGS) -MM $$f >> $(DEPEND); \
	done

-include $(DEPEND)




