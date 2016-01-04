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

ifeq ($(findstring release,$(MAKECMDGOALS)),)
  CFLAGS = -g -DDEBUG $(BASEFLAGS)
  LDFLAGS = -g
  #LDFLAGS = -g -Wl,--eh-frame-hdr -L/opt/geo/lib
  LIBS = 
  # use these instead, for HEAP Profiling with gproftools 
  #CFLAGS = -g -DDEBUG -DHEAPROFILE -I/opt/geo/include $(BASEFLAGS)
  #LDFLAGS = -g -Wl,--eh-frame-hdr -L/opt/geo/lib
  #LIBS = -ltcmalloc
else
  CFLAGS = -O2 -DNDEBUG $(BASEFLAGS)
  LDFLAGS = 
  LIBS = 
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
${GCLDIR}/gff.o  : ${GCLDIR}/gff.h
./gtf_tracking.o : ./gtf_tracking.h
./gffcompare.o : ./gtf_tracking.h
gffcompare: ${OBJS} ./gffcompare.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}

.PHONY : clean
clean:: 
	@${RM} core core.* gffcompare gffcompare.exe ${OBJS} *.o* 


