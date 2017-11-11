OBJSE = anread.o sndhdr.o wavhdr.o byteorder.o header.o  plotseg.o g_raph.a
OBJST = sndhdr.o wavhdr.o byteorder.o plotseg.o fft2.o g_raph.a
INCL = header.h byteorder.h sndhdr.h wavhdr.h macro.h

all: extend addsyn test fft2.o anread.o header.o plotseg.o

extend: addsynextend.c $(OBJSE) $(INCL)
		cc -g -m32 -o extend addsynextend.c $(OBJSE)
addsyn: addsyn.c $(OBJSE) $(INCL)
		cc -g -m32 -o addsyn addsyn.c $(OBJSE)
clean:
		rm extend addsyn
