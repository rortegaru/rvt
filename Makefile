FCMP=gfortran -g   -ffixed-line-length-none
GCMP=gcc -g 
#FCMP=f77 
MCHDEP= mchdep.o
MCHCMD= mnmarg.o mgtarg.o
RM = rm -f

##### DO NOT CHANGE #####

all: tdcal rvcal dorvmat fscal tscal xmultiplerv

tdcal:	tdcal.o prop.o tdsim.o filt.o  lgstr.o \
		boore.o qinv.o recipes.o $(MCHCMD) $(MCHARG)
	$(FCMP) tdcal.o prop.o tdsim.o filt.o  \
		boore.o qinv.o recipes.o lgstr.o $(MCHCMD) $(MCHARG) -o tdcal

rvcal:	rvcal.o prop.o rvsim.o filt.o  lgstr.o \
		boore.o qinv.o recipes.o $(MCHCMD) $(MCHARG)
	$(FCMP) rvcal.o prop.o rvsim.o filt.o  \
		boore.o qinv.o recipes.o lgstr.o $(MCHCMD) $(MCHARG) -o rvcal

dorvmat:	dorvmat.o rvmat.o prop.o rvsim.o filt.o  lgstr.o \
		boore.o qinv.o recipes.o $(MCHCMD) $(MCHARG)
	$(FCMP) dorvmat.o rvmat.o prop.o rvsim.o filt.o  \
		boore.o qinv.o recipes.o lgstr.o $(MCHCMD) $(MCHARG) -o dorvmat

fscal:	fscal.o prop.o fssim.o filt.o  lgstr.o \
		boore.o qinv.o recipes.o $(MCHCMD) $(MCHARG)
	$(FCMP) fscal.o prop.o fssim.o filt.o  \
		boore.o qinv.o recipes.o lgstr.o $(MCHCMD) $(MCHARG) -o fscal

tscal:	tscal.o prop.o tssim.o filt.o  lgstr.o boore.o \
		qinv.o recipes.o sacsubf.o $(MCHCMD) $(MCHARG)
	$(FCMP) tscal.o prop.o tssim.o filt.o  boore.o \
		qinv.o recipes.o lgstr.o sacsubf.o $(MCHCMD) $(MCHARG) -o tscal

xmultiplerv: xmultiplerv.o multiplerv.o prop.o rvsim.o filt.o lgstr.o ran0.o \
		boore.o qinv.o recipes.o $(MCHCMD) $(MCHRG) 
	$(FCMP) xmultiplerv.o multiplerv.o prop.o rvsim.o filt.o ran0.o \
		boore.o qinv.o recipes.o lgstr.o $(MCHCMD) $(MCHARG) -o xmultiplerv

clean:
	${RM} *.o


.c.o:
	$(GCMP) -c $<

.f.o:
	$(FCMP) -c $<
