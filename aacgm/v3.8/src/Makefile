# Makefile to create the AACGM dlm .so file.
# This is for MAC but should be easily modified for linux too.
# 
# Set the ARCH variable to either "i386" or "x86_64"
# Set the IDL_DIR appropriately
#
# DLG
# greendl1@ornl.gov
#

ARCH = x86_64 
CFLAGS = -arch ${ARCH} -fPIC -dynamic -fno-common -D_REENTRANT #-no-cpp-precomp 
IDL_DIR="/usr/local/exelis/idl82"
INC_DIR = -I ${IDL_DIR}/external/include
CC = gcc

OBJS = $(patsubst %.c,%.o,$(wildcard *.c))


LIBS = -L${IDL_DIR}/bin/bin.darwin.${ARCH} -lidl 

all: ${OBJS}
	${CC} -shared ${CFLAGS} *.o ${LIBS} ${INC_DIR} -o idl_aacgm.so
	cp idl_aacgm.so ../

%.o: %.c
	${CC} -c ${CFLAGS} $< -o $@ ${INC_DIR}

clean:
	rm *.o *.so



