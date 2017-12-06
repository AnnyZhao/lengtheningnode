OBJSE = anread.o sndhdr.o wavhdr.o byteorder.o header.o fft2.o plotseg.o g_raph.a
INCL = header.h byteorder.h sndhdr.h wavhdr.h macro.h g_raph.h

extend: extend.c $(OBJSE) $(INCL)
	cc -m32 -o extend extend.c $(OBJSE)
plotsegp.o: plotsegp.c $(INCL)
	cc -m32 -c plotsegp.c
plotseg.o: plotseg.c $(INCL)
	cc -m32 -c plotseg.c

