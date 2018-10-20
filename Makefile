#GCLIB := ../gclib

GCLIB := $(if $(GCLIB),$(GCLIB),../gclib)

INCDIRS := -I${GCLIB}

BASEFLAGS  = -Wall -Wextra ${INCDIRS} -std=c++11 -D_REENTRANT -fno-exceptions -fno-rtti

GCCV8 := $(shell expr `g++ -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCV8)" "1"
 BASEFLAGS += -Wno-class-memaccess
endif

CXX   := $(if $(CXX),$(CXX),g++)
LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

CXXFLAGS := $(if $(CXXFLAGS),$(BASEFLAGS) $(CXXFLAGS),$(BASEFLAGS))

#

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  LIBS := 
  ifneq (,$(findstring static,$(MAKECMDGOALS)))
    LDFLAGS += -static-libstdc++ -static-libgcc
  endif
  CXXFLAGS := -O3 -DNDEBUG $(CXXFLAGS)
else # debug build
  CXXFLAGS += -g -O0 -DDEBUG -D_DEBUG -DGDEBUG
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     #make memcheck : use the statically linked address sanitizer in gcc 4.9.x
     GCCVER49 := $(shell expr `g++ -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address
     GCCVER5 := $(shell expr `g++ -dumpversion | cut -f1 -d.` \>= 5)
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     #CFLAGS += $(BASEFLAGS)
     CXXFLAGS += -fno-common -fstack-protector
     #LIBS := -Wl,-Bstatic -lasan -lubsan -Wl,-Bdynamic -ldl $(LIBS)
     LIBS := -lasan -lubsan -ldl $(LIBS)
  endif
endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := g++

OBJS = ${GCLIB}/GFastaIndex.o ${GCLIB}/GFaSeqGet.o ${GCLIB}/gff.o \
 ${GCLIB}/gdna.o ${GCLIB}/codons.o ${GCLIB}/GBase.o \
 ${GCLIB}/GStr.o ${GCLIB}/GArgs.o

.PHONY : all
all:    gffcompare trmap
debug:  gffcompare trmap
release: gffcompare trmap
static: gffcompare trmap
memcheck: gffcompare trmap
memdebug: gffcompare trmap

${GCLIB}/gff.o  : ${GCLIB}/gff.h
./gtf_tracking.o : ./gtf_tracking.h
./gffcompare.o : ./gtf_tracking.h

gffcompare: ${OBJS} ./t_classify.o ./gtf_tracking.o ./gffcompare.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

trmap: ${OBJS} ./t_classify.o ./GIntervalTree.o ./trmap.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

.PHONY : clean
clean:: 
	@${RM} core core.* gffcompare gffcompare.exe trmap trmap.exe ${OBJS} *.o* 
