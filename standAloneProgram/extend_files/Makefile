OBJSE = anread.o sndhdr.o wavhdr.o byteorder.o header.o  plotseg.o g_raph.a
INCL = header.h byteorder.h sndhdr.h wavhdr.h macro.h

extend: extend.c $(OBJSE) $(INCL)
		cc -m32 -o extend extend.c $(OBJSE)
clean:
		rm extend fft2.o anread.o header.o plotseg.o
