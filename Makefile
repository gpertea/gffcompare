# Useful directories

THISCODEDIR := .
GCLDIR := ../gclib
# Directory where libz.a can be found
# (please build that first before making this package) 
# ZDIR := ../zlib
# Directories to search for header files
#SEARCHDIRS := -I${ZDIR} -I${THISCODEDIR} -I${GCLDIR}
SEARCHDIRS := -I${THISCODEDIR} -I${GCLDIR}

SYSTYPE :=     $(shell uname)

# C compiler

MACHTYPE :=     $(shell uname -m)
ifeq ($(MACHTYPE), i686)
    MARCH = -march=i686 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
else
    MARCH = 
endif    

# CVS checked in
CC      := g++
BASEFLAGS  = -Wall -Wextra ${SEARCHDIRS} $(MARCH) \
 -fno-exceptions -fno-rtti -D_REENTRANT -D_DARWIN_C_SOURCE

#

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  CFLAGS = -O3 -DNDEBUG $(BASEFLAGS)
  LDFLAGS = 
  LIBS = 
  ifneq (,$(findstring static,$(MAKECMDGOALS)))
    LDFLAGS += -static-libstdc++ -static-libgcc
  endif
else # debug build
  ifneq (,$(filter %memcheck %memdebug, $(MAKECMDGOALS)))
     #make memcheck : use the statically linked address sanitizer in gcc 4.9.x
     GCCVER49 := $(shell expr `g++ -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CFLAGS := -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address
     GCCVER5 := $(shell expr `g++ -dumpversion | cut -f1 -d.` \>= 5)
     ifeq "$(GCCVER5)" "1"
       CFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CFLAGS += $(BASEFLAGS)
     CFLAGS := -g -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector $(CFLAGS)
     LDFLAGS := -g
     #LIBS := -Wl,-Bstatic -lasan -lubsan -Wl,-Bdynamic -ldl $(LIBS)
     LIBS := -lasan -lubsan -ldl $(LIBS)
  else
    #ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    #   BASEFLAGS += -DGMEMTRACE
    #   GMEMTRACE=1
    #endif
    #--- just plain debug build ---
     CFLAGS = -g -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
     LDFLAGS = -g
     LIBS = 
  endif
endif

%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cc
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.C
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cxx
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := g++

OBJS = ${GCLDIR}/GFastaIndex.o ${GCLDIR}/GFaSeqGet.o ${GCLDIR}/gff.o \
 ./gtf_tracking.o ${GCLDIR}/gdna.o ${GCLDIR}/codons.o ${GCLDIR}/GBase.o \
 ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o

.PHONY : all
all:    gffcompare
debug:  gffcompare
release: gffcompare
static: gffcompare
memcheck: gffcompare
memdebug: gffcompare

${GCLDIR}/gff.o  : ${GCLDIR}/gff.h
./gtf_tracking.o : ./gtf_tracking.h
./gffcompare.o : ./gtf_tracking.h
gffcompare: ${OBJS} ./gffcompare.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

.PHONY : clean
clean:: 
	@${RM} core core.* gffcompare gffcompare.exe ${OBJS} *.o* 


