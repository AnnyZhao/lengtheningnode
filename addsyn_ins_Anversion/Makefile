#	Makefile for keeping additive synthesis instrument up to date
#	Written by: James Beauchamp, 1993
#
#  06/05/12 jwb Modified for OSX 10.6 installation

CFLAGS = -O -arch i386 -c

all: test.wav

clean:
	rm -f *.o m4c.* *.list as.c test.wav S1* S2*

test.wav: m4c.as test.sc
	rm -f test.wav .list
	$(M4CDIR)/util/runm4c m4c.as test.wav test.sc test.list

m4c.as: as_orch.o as.o
	$(M4CDIR)/util/mx as_orch.o as.o
	mv m4c m4c.as

as.o: as.c anread.c header.c header.h macro.h
	cc $(CFLAGS) as.c

as.c: as.t addsynfunc.c
	$(M4CDIR)/util/insdes as.t
