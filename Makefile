# Author: Yipeng Sun, Alex Fernez

#################
# Configuration #
#################

BINPATH 	:=	bin
VPATH		:=	include:src
CPP_FILES	:=	$(wildcard src/*.cpp)
EXE_FILES	:=	$(patsubst src/%.cpp,$(BINPATH)/%.exe,$(CPP_FILES))

# Compiler settings
COMPILER		:=	$(shell root-config --cxx)
CXXFLAGS		:=	$(shell root-config --cflags)
LINKFLAGS		:=	$(shell root-config --libs)
ADDCXXFLAGS		:=	-O2 -march=native -mtune=native -Iinclude
ADDLINKFLAGS	:=	-lRooFitCore -lRooFit -lRooStats -lHistFactory 

OS := $(shell uname)
ifeq ($(OS),Darwin)
  $(info OS is $(OS) (macOS), adding -lc++fs to LINKFLAGS)
  LINKFLAGS := $(shell root-config --libs) -lc++fs
endif



########
# Misc #
########
.PHONY: exe clean

exe: $(EXE_FILES)

clean:
	@rm -rf ./$(BINPATH)/*.exe
	@rm -rf ./gen/*
	@remove_backups.sh



###############
# Compile C++ #
###############

$(BINPATH)/%.exe: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) $(ADDLINKFLAGS)
