#Sharon needs to use the following:
CC=c++
CCFLAGS=-Wall -pedantic -O3 

LDFLAGS=-lm

NAMES=readin

OBJ=Margin.o bayes.o nrutil.o # defining the objects
ROBJ=Readin.o fitness.o bayes.o nrutil.o # defining the objects
POBJ=ProbCont.o nrutil.o
MOBJ=Metrop.o nrutil.o bayes.o
GOBJ=GaussTest1d.o nrutil.o bayes.o
LOBJ=LogLike1d.o nrutil.o bayes.o
GGOBJ=GaussTest2d.o nrutil.o bayes.o
MMOBJ=Metrop2d.o nrutil.o bayes.o

all: $(NAMES)

margin: Margin.o $(OBJ)
	$(CC) $(CCFLAGS) $(OBJ) $(LDFLAGS) -o margin

readin: Readin.o $(ROBJ)
	$(CC) $(CCFLAGS) $(ROBJ) $(LDFLAGS) -o readin

contours: ProbCont.o $(POBJ)
	$(CC) $(CCFLAGS) $(POBJ) $(LDFLAGS) -o contours

metrop: Metrop.o $(MOBJ)
	$(CC) $(CCFLAGS) $(MOBJ) $(LDFLAGS) -o metrop

gauss1d: GaussTest1d.o $(GOBJ)
	$(CC) $(CCFLAGS) $(GOBJ) $(LDFLAGS) -o gauss1d

loglike1d: LogLike1d.o $(LOBJ)
	$(CC) $(CCFLAGS) $(LOBJ) $(LDFLAGS) -o loglike1d

gauss2d: GaussTest2d.o $(GGOBJ)
	$(CC) $(CCFLAGS) $(GGOBJ) $(LDFLAGS) -o gauss2d

metrop2d: Metrop2d.o $(MMOBJ)
	$(CC) $(CCFLAGS) $(MMOBJ) $(LDFLAGS) -o metrop2d

Margin.o: \
	Margin.cpp \
	nrutil.h \
	Makefile
	$(CC) $(CCFLAGS) -c Margin.cpp

Readin.o: \
	Readin.cpp \
	nrutil.h \
	Struct.h \
	fitness.h \
	bayes.h \
	Makefile
	$(CC) $(CCFLAGS) -c Readin.cpp

ProbCont.o: \
	ProbCont.cpp \
	nrutil.h \
	Makefile
	$(CC) $(CCFLAGS) -c ProbCont.cpp

Metrop.o: \
	Metrop.cpp \
	nrutil.h \
	Makefile
	$(CC) $(CCFLAGS) -c Metrop.cpp

GaussTest1d.o: \
	GaussTest1d.cpp \
	nrutil.h \
	Makefile
	$(CC) $(CCFLAGS) -c GaussTest1d.cpp

LogLike1d.o: \
	LogLike1d.cpp \
	nrutil.h \
	Makefile
	$(CC) $(CCFLAGS) -c LogLike1d.cpp


GaussTest2d.o: \
	GaussTest2d.cpp \
	nrutil.h \
	Makefile
	$(CC) $(CCFLAGS) -c GaussTest2d.cpp

Metrop2d.o: \
	Metrop2d.cpp \
	nrutil.h \
	Makefile
	$(CC) $(CCFLAGS) -c Metrop2d.cpp

nrutil.o: \
	nrutil.h \
	nrutil.c
	$(CC) $(CCFLAGS) -c nrutil.c

bayes.o: \
	bayes.cpp \
	bayes.h
	$(CC) $(CCFLAGS) -c bayes.cpp

fitness.o: \
	fitness.cpp \
	fitness.h \
	Struct.h 
	$(CC) $(CCFLAGS) -c fitness.cpp

clean:
	rm -f core *~ $(OBJ) $(APPOBJ)

veryclean:
	rm -f core *~ $(OBJ) $(APPOBJ) $(NAMES)
