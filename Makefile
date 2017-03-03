OBJSE = anread.o sndhdr.o wavhdr.o byteorder.o header.o
OBJST = sndhdr.o wavhdr.o byteorder.o plotseg.o fft2.o g_raph.a
INCL = header.h byteorder.h sndhdr.h wavhdr.h macro.h

extend: addsynextend.c $(OBJSE) $(INCL)
		cc -m32 -o extend addsynextend.c $(OBJSE)
addsyn: addsyn.c $(OBJSE) $(INCL)
		cc -m32 -o addsyn addsyn.c $(OBJSE)
test: test2.c $(OBJST) $(INCL)
		cc -m32 -o test test2.c $(OBJST)
fft2.o: fft2.c
		cc -m32 -c fft2.c
anread.o: anread.c
		cc -m32 -c anread.c
header.o: header.c
		cc -m32 -c header.c
plotseg.o: plotseg.c g_raph.h macro.h
		cc -m32 plotseg.c
