## need module load mkl gcc48

#INCL = -I/usr/include/gcc/darwin/3.3-fast \
	 -I/System/Library/Frameworks/vecLib.framework/Headers/

INCL = -I$(MKL_HOME)/include

LIB = -L$(MKL_HOME)/lib/intel64

CFLAGS =    -O3  -funroll-loops -fstrict-aliasing -mtune=native -m64 -ffast-math #-Wno-long-double #-DNDEBUG 

LIBS =  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lpthread -lm  -lstdc++

CC       = g++

MPE_C =       mpe.cpp protein.cpp BD.cpp pdb.cpp mcoeff.cpp gradcoeff.cpp \
	      xform.cpp shcoeff.cpp tricoeff.cpp torqcoeff.cpp rotcoeff.cpp \
		transcoeff.cpp

FIT_C =       mcoeff.cpp gradcoeff.cpp nrc_svd.cpp polarfit.cpp \
	      shcoeff.cpp tricoeff.cpp fit.cpp transcoeff.cpp

PARAMFIT_C =   mcoeff.cpp gradcoeff.cpp nrc_svd.cpp paramfit.cpp \
	      shcoeff.cpp tricoeff.cpp 	

TESTFIT_C =   expansion.cpp molecule.cpp testfit.cpp lautil.cpp \
		rotationbase.cpp rotation.cpp drottheta.cpp drotphi.cpp


MPE_OBJS      =       $(MPE_C:.cpp=.o)
FIT_OBJS      =       $(FIT_C:.cpp=.o)
PARAMFIT_OBJS      =       $(PARAMFIT_C:.cpp=.o)
TESTFIT_OBJS      =       $(TESTFIT_C:.cpp=.o)

mpe:  $(MPE_OBJS) Makefile
	rm -f ../bin/mpe
	$(CC)  $(CFLAGS) $(MPE_OBJS) $(LIB) $(LIBS) -o mpe
	mv mpe ../bin/

fit:  $(FIT_OBJS) Makefile
	rm -f ../bin/fit
	$(CC)  $(CFLAGS) $(FIT_OBJS) $(LIB) $(LIBS) -o fit
	mv fit ../bin

paramfit:  $(PARAMFIT_OBJS) Makefile
	rm -f ../bin/paramfit
	$(CC)  $(CFLAGS) $(PARAMFIT_OBJS) $(LIB) $(LIBS) -o paramfit -framework veclib
	mv paramfit ../bin

testfit:  $(TESTFIT_OBJS) Makefile
	rm -f ../bin/testfit
	$(CC)  $(CFLAGS) $(TESTFIT_OBJS) $(LIB) $(LIBS) -o testfit -framework veclib
	mv testfit ../bin/

.cpp.o:
	$(CC) $(CFLAGS) $(INCL) -c $<

clean:
	rm *.o

depend:
	makedepend $(INCL) $(MPE_C) $(FIT_C) $(PARAMFIT_C) $(TESTFIT_C)


# DO NOT DELETE THIS LINE -- make depend depends on it.

#mpe.o: /usr/include/sys/time.h /usr/include/sys/appleapiopts.h
#mpe.o: /usr/include/sys/types.h /usr/include/sys/cdefs.h
#mpe.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#mpe.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#mpe.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#mpe.o: /usr/include/time.h /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#mpe.o: /usr/include/ppc/limits.h
mpe.o: mpe.h util.h mcoeff.h gradcoeff.h tricoeff.h
mpe.o: torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h protein.h pdb.h   
protein.o: protein.h pdb.h util.h 
#protein.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h /usr/include/sys/time.h
#protein.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#protein.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#protein.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#protein.o: /usr/include/ppc/endian.h /usr/include/time.h
#protein.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#protein.o: /usr/include/ppc/limits.h
protein.o: torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h  mpe.h mcoeff.h gradcoeff.h tricoeff.h
BD.o: BD.h util.h 
#/usr/include/sys/time.h 
#BD.o: /usr/include/sys/types.h /usr/include/sys/cdefs.h /usr/include/sys/appleapiopts.h
#BD.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#BD.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#BD.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#BD.o: /usr/include/time.h /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#BD.o: /usr/include/ppc/limits.h 
BD.o: pdb.h protein.h mpe.h mcoeff.h gradcoeff.h
BD.o: tricoeff.h torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h
pdb.o: pdb.h util.h 
#/usr/include/sys/time.h /usr/include/sys/appleapiopts.h
#pdb.o: /usr/include/sys/types.h /usr/include/sys/cdefs.h
#pdb.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#pdb.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#pdb.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#pdb.o: /usr/include/time.h /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#pdb.o: /usr/include/ppc/limits.h 
pdb.o: protein.h mpe.h mcoeff.h gradcoeff.h
pdb.o: tricoeff.h torqcoeff.h xform.h shcoeff.h rotcoeff.h transcoeff.h
mcoeff.o: mcoeff.h util.h 
#/usr/include/sys/time.h
#mcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#mcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#mcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#mcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#mcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#mcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#mcoeff.o: /usr/include/ppc/limits.h 
mcoeff.o: shcoeff.h
gradcoeff.o: gradcoeff.h tricoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#gradcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#gradcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#gradcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#gradcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#gradcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#gradcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#gradcoeff.o: /usr/include/ppc/limits.h
xform.o: xform.h mcoeff.h util.h 
#/usr/include/sys/time.h
#xform.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#xform.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#xform.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#xform.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#xform.o: /usr/include/ppc/endian.h /usr/include/time.h
#xform.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#xform.o: /usr/include/ppc/limits.h 
xform.o: gradcoeff.h tricoeff.h torqcoeff.h
xform.o: shcoeff.h rotcoeff.h transcoeff.h
shcoeff.o: shcoeff.h util.h 
#/usr/include/sys/time.h
#shcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#shcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#shcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#shcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#shcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#shcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#shcoeff.o: /usr/include/ppc/limits.h mcoeff.h
tricoeff.o: tricoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#tricoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#tricoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#tricoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#tricoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#tricoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#tricoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#tricoeff.o: /usr/include/ppc/limits.h
torqcoeff.o: torqcoeff.h tricoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#torqcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#torqcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#torqcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#torqcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#torqcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#torqcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#torqcoeff.o: /usr/include/ppc/limits.h 
torqcoeff.o: shcoeff.h
rotcoeff.o: rotcoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#rotcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#rotcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#rotcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#rotcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#rotcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#rotcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#rotcoeff.o: /usr/include/ppc/limits.h 
rotcoeff.o: tricoeff.h shcoeff.h
transcoeff.o: shcoeff.h util.h 
#/usr/include/sys/time.h
#transcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#transcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#transcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#transcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#transcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#transcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#transcoeff.o: /usr/include/ppc/limits.h 
transcoeff.o: mcoeff.h transcoeff.h
mcoeff.o: mcoeff.h util.h 
#/usr/include/sys/time.h
#mcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#mcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#mcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#mcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#mcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#mcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#mcoeff.o: /usr/include/ppc/limits.h 
mcoeff.o: shcoeff.h
gradcoeff.o: gradcoeff.h tricoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#gradcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#gradcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#gradcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#gradcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#gradcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#gradcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#gradcoeff.o: /usr/include/ppc/limits.h
#nrc_svd.o: /usr/include/unistd.h /usr/include/sys/cdefs.h
#nrc_svd.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
#nrc_svd.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#nrc_svd.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#nrc_svd.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#nrc_svd.o: /usr/include/sys/unistd.h /usr/include/signal.h
#nrc_svd.o: /usr/include/sys/signal.h /usr/include/machine/signal.h
#nrc_svd.o: /usr/include/ppc/signal.h nrc_svd.h
#polarfit.o: mcoeff.h util.h /usr/include/sys/time.h
#polarfit.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#polarfit.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#polarfit.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#polarfit.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#polarfit.o: /usr/include/ppc/endian.h /usr/include/time.h
#polarfit.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#polarfit.o: /usr/include/ppc/limits.h transcoeff.h shcoeff.h nrc_svd.h fit.h
shcoeff.o: shcoeff.h util.h 
#/usr/include/sys/time.h
#shcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#shcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#shcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#shcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#shcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#shcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#shcoeff.o: /usr/include/ppc/limits.h 
shcoeff.o: mcoeff.h
tricoeff.o: tricoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#tricoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#tricoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#tricoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#tricoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#tricoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#tricoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#tricoeff.o: /usr/include/ppc/limits.h
#fit.o: nrc_svd.h mcoeff.h util.h /usr/include/sys/time.h
#fit.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#fit.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#fit.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#fit.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#fit.o: /usr/include/ppc/endian.h /usr/include/time.h
#fit.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#fit.o: /usr/include/ppc/limits.h shcoeff.h fit.h
transcoeff.o: shcoeff.h util.h 
#/usr/include/sys/time.h
#transcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#transcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#transcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#transcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#transcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#transcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#transcoeff.o: /usr/include/ppc/limits.h 
transcoeff.o: mcoeff.h transcoeff.h
mcoeff.o: mcoeff.h util.h 
#/usr/include/sys/time.h
#mcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#mcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#mcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#mcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#mcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#mcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#mcoeff.o: /usr/include/ppc/limits.h 
mcoeff.o: shcoeff.h
gradcoeff.o: gradcoeff.h tricoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#gradcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#gradcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#gradcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#gradcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#gradcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#gradcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#gradcoeff.o: /usr/include/ppc/limits.h
#nrc_svd.o: /usr/include/unistd.h /usr/include/sys/cdefs.h
#nrc_svd.o: /usr/include/sys/types.h /usr/include/sys/appleapiopts.h
#nrc_svd.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#nrc_svd.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#nrc_svd.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#nrc_svd.o: /usr/include/sys/unistd.h /usr/include/signal.h
#nrc_svd.o: /usr/include/sys/signal.h /usr/include/machine/signal.h
#nrc_svd.o: /usr/include/ppc/signal.h nrc_svd.h
#paramfit.o: nrc_svd.h mcoeff.h util.h /usr/include/sys/time.h
#paramfit.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#paramfit.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#paramfit.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#paramfit.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#paramfit.o: /usr/include/ppc/endian.h /usr/include/time.h
#paramfit.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#paramfit.o: /usr/include/ppc/limits.h shcoeff.h
#paramfit.o: /System/Library/Frameworks/vecLib.framework/Headers/clapack.h
shcoeff.o: shcoeff.h util.h 
#/usr/include/sys/time.h
#shcoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#shcoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#shcoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#shcoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#shcoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#shcoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#shcoeff.o: /usr/include/ppc/limits.h 
shcoeff.o: mcoeff.h
tricoeff.o: tricoeff.h mcoeff.h util.h 
#/usr/include/sys/time.h
#tricoeff.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#tricoeff.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#tricoeff.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#tricoeff.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#tricoeff.o: /usr/include/ppc/endian.h /usr/include/time.h
#tricoeff.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#tricoeff.o: /usr/include/ppc/limits.h
#expansion.o: expansion.h util.h /usr/include/sys/time.h
#expansion.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#expansion.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#expansion.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#expansion.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#expansion.o: /usr/include/ppc/endian.h /usr/include/time.h
#expansion.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#expansion.o: /usr/include/ppc/limits.h
#molecule.o: molecule.h expansion.h util.h /usr/include/sys/time.h
#molecule.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#molecule.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#molecule.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#molecule.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#molecule.o: /usr/include/ppc/endian.h /usr/include/time.h
#molecule.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#molecule.o: /usr/include/ppc/limits.h
#testfit.o: nrc_svd.h molecule.h expansion.h util.h /usr/include/sys/time.h
#testfit.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#testfit.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#testfit.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#testfit.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#testfit.o: /usr/include/ppc/endian.h /usr/include/time.h
#testfit.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#testfit.o: /usr/include/ppc/limits.h rotation.h rotationbase.h drottheta.h
#testfit.o: drotphi.h
#lautil.o: /System/Library/Frameworks/vecLib.framework/Headers/clapack.h
#rotationbase.o: rotationbase.h expansion.h util.h /usr/include/sys/time.h
#rotationbase.o: /usr/include/sys/appleapiopts.h /usr/include/sys/types.h
#rotationbase.o: /usr/include/sys/cdefs.h /usr/include/machine/types.h
#rotationbase.o: /usr/include/ppc/types.h /usr/include/machine/ansi.h
#rotationbase.o: /usr/include/ppc/ansi.h /usr/include/machine/endian.h
#rotationbase.o: /usr/include/ppc/endian.h /usr/include/time.h
#rotationbase.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#rotationbase.o: /usr/include/ppc/limits.h
#rotation.o: rotation.h rotationbase.h expansion.h util.h
#rotation.o: /usr/include/sys/time.h /usr/include/sys/appleapiopts.h
#rotation.o: /usr/include/sys/types.h /usr/include/sys/cdefs.h
#rotation.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#rotation.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#rotation.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#rotation.o: /usr/include/time.h
#rotation.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#rotation.o: /usr/include/ppc/limits.h drottheta.h drotphi.h
#drottheta.o: drottheta.h rotationbase.h expansion.h util.h
#drottheta.o: /usr/include/sys/time.h /usr/include/sys/appleapiopts.h
#drottheta.o: /usr/include/sys/types.h /usr/include/sys/cdefs.h
#drottheta.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#drottheta.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#drottheta.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#drottheta.o: /usr/include/time.h
#drottheta.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#drottheta.o: /usr/include/ppc/limits.h rotation.h drotphi.h
#drotphi.o: drotphi.h rotationbase.h expansion.h util.h
#drotphi.o: /usr/include/sys/time.h /usr/include/sys/appleapiopts.h
#drotphi.o: /usr/include/sys/types.h /usr/include/sys/cdefs.h
#drotphi.o: /usr/include/machine/types.h /usr/include/ppc/types.h
#drotphi.o: /usr/include/machine/ansi.h /usr/include/ppc/ansi.h
#drotphi.o: /usr/include/machine/endian.h /usr/include/ppc/endian.h
#drotphi.o: /usr/include/time.h
#drotphi.o: /usr/include/gcc/darwin/3.3-fast/machine/limits.h
#drotphi.o: /usr/include/ppc/limits.h rotation.h drottheta.h
