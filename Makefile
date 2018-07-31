EXEC   = 2LPTicLCDM

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o  \
         read_table.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile



#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing


#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type
OPT   +=  -DHUBBLE_TABLE # switch this on if you want tabulated H diagram
#OPT   +=  -DHUBBLE_USER  # switch this on if you want to give your own H/H0 value
OPT   +=  -DDMMASS_TABLE # switch this on if you want tabulated DM mass
#OPT   +=  -DOMEGA_USER   # switch this on if you want to give your own omega_matter value at the starting redshift
#OPT   +=  -DACCURATE_DA  # allow user to specify dD/da value for displacement velocity

OPT   +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles start from a glass (as opposed to grid)

#OPT += -DONLY_ZA # swith this on if you want ZA initial conditions (2LPT otherwise)



OPTIONS =  $(OPT)


FFTW_INCL = -I/home/jjzhang/Library/fftw/include
FFTW_LIBS = -L/home/jjzhang/Library/fftw/lib

GSL_LIBS=   -L/home/jjzhang/Library/gsl/lib
GSL_INCL =  -I/home/jjzhang/Library/gsl/include


CC       = mpicc
MPICHLIB = -L/home/jjzhang/Library/openmpi/lib -lmpi

OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)


FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas

CFLAGS =   $(OPTIONS)  $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)



