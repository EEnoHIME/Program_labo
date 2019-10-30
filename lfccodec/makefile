CC = gcc
ifeq ($(CC), icc)
  #CFLAG = -help
  CFLAG = -Wall -O3 -unroll -parallel -par-threshold100 -g -qopt-report -mcmodel=medium -shared-intel
else
  CFLAG = -O4 -funroll-loops -Wall
  ifeq ($(shell uname -m), i686)
    CFLAG += -mpentiumpro
  endif
endif

ENCODER = encmrp
ENCOBJ = encmrp.o common.o rc.o
DECODER = decmrp
DECOBJ = decmrp.o common.o rc.o


all : $(ENCODER) $(DECODER)

$(ENCODER) : $(ENCOBJ)
	$(CC) $(CFLAG) -o $@ $(ENCOBJ) -lm

$(DECODER) : $(DECOBJ)
	$(CC) $(CFLAG) -o $@ $(DECOBJ) -lm

.c.o :
	$(CC) $(CFLAG) -c $<

encmrp.o : encmrp.c mrp.h
decmrp.o : decmrp.c mrp.h
common.o : common.c mrp.h
rc.o : rc.c mrp.h

clean :
	rm -f $(ENCODER) $(DECODER) $(ENCOBJ) $(DECOBJ) *~
