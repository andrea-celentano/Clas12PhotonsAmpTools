
INC_DIR :=  -I. -I$(AMPTOOLS) -I../include

CXX := gcc
CXX_FLAGS := $(shell root-config --cflags)

SRCDIRS := src 
TARGET_LIBS := $(addsuffix .a, $(addprefix lib, $(SRCDIRS)))

SRCDIRS_GPU := $(SRCDIRS)
TARGET_LIBS_GPU :=  $(addsuffix _GPU.a, $(addprefix lib, $(SRCDIRS_GPU)))

#To build GPU-accelerated code type: make GPU=1
ifdef GPU

NVCC := nvcc
CUDA_FLAGS := -arch=compute_13 -code=compute_13
INC_DIR += -I. -I$(CUDA_INSTALL_PATH)/include

CXX_FLAGS += -DGPU_ACCELERATION
DEFAULT := libClas12PhotonsAmpTools_GPU.a

else

DEFAULT := libClas12PhotonsAmpTools.a

endif

export

.PHONY: default clean

default: lib $(DEFAULT)

lib:
	mkdir lib

libClas12PhotonsAmpTools.a: $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd lib; ar -x $(lib) ) )
	@cd lib && ar -rv $@ *.o
	@cd lib && rm -f *.o

libClas12PhotonsAmpTools_GPU.a: $(TARGET_LIBS_GPU)
	$(foreach lib_GPU, $(TARGET_LIBS_GPU), $(shell cd lib; ar -x $(lib_GPU) ) )
	@cd lib && ar -rv $@ *.o
	@cd lib && rm -f *.o

lib%_GPU.a: 
	@$(MAKE) -C $(subst lib,, $(subst _GPU.a,, $@ )) LIB=$@
	@cp $(subst lib,, $(subst _GPU.a,, $@))/$@ lib/

lib%.a: 
	@$(MAKE) -C $(subst lib,, $(subst .a,, $@ )) LIB=$@
	@cp $(subst lib,, $(subst .a,, $@))/$@ lib/

clean: $(addprefix clean_, $(SRCDIRS))
	-rm -f lib/*.a

clean_%:
	@-cd $(subst clean_,, $@) && $(MAKE) clean
