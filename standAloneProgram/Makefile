OBJSE = anread.o sndhdr.o wavhdr.o byteorder.o header.o fft2.o plotseg.o g_raph.a
INCL = header.h byteorder.h sndhdr.h wavhdr.h macro.h

addsynextend: addsynextend.c $(OBJSE) $(INCL)
		cc -m32 -o addsynextend addsynextend.c $(OBJSE)
clean:
		rm addsynextend fft2.o anread.o header.o plotseg.o

anread.o: anread.c $(INCL)
		cc -m32 -c anread.c -o anread.o

header.o: header.c $(INCL)
		cc -m32 -c header.c -o header.o

plotseg.o: plotseg.c $(INCL)
		cc -m32 -c plotseg.c -o plotseg.o

fft2.o: fft2.c $(INCL)
		cc -m32 -c fft2.c -o fft2.o
