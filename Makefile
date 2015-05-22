#Unix makefile for fortran-file	

# Parameters
# name of the target program here
MAKEFILE = Makefile
EXE = monolayer

FC = mpif90 #${F90}
#FC = gfortran
#FC = /shared/software/openmpi-1.6.1/bin/mpif77
#FC = /shared/software/openmpi-1.6.1/bin/mpif90

# This flags are used in the compilation stage (name should be CFLAGS)
# To debug:
FFLAGS= -cpp -g -p -fbacktrace -fcheck=all -Wall  
## -g produce debugging information in the operating system's native format.
## -pg generate extra code to write profile information suitable for the analysis program gprof
# To run
# -fpp: Calls first to the C Preprocessor
#FFLAGS= -cpp -O2 -fno-toplevel-reorder
#FFLAGS= -cpp -O3 -fno-toplevel-reorder

SRC = module_globales.f90 \
      module_Csys.f90 \
      module_FreeEnergy.f90 \
      module_pore.f90 \
      set_bulk_properties.f90 \
      units_adaptation.f90 monolayer.f90 \
      read_input.f90 \
      savetodisk.f90 \
      checknum.f90 \
      save_data.f90 \
      open_files.f90 \
      close_files.f90 \
      imprimir_cadenas.f90 \
      printstate.f90\
      check_run.f90 \
      creador.f90 pxs.f90 kai.f90 \
      cadenas72mr.f90 rota36.f90\
      mrrrr.f90 rands.f90 \
      call_kinsol.f90 fkfun.f90 factorcurv.f90 \
      fkpset.f90 fkpsol.f90 \
      allocating.f90 \
      set_initial_guess.f90 \
      set_pore_distrib.f90 \
      calc_conductance.f90 \
      calc_mean_values.f90 \
      calc_energy.f90 \
      Fmix.f90 Fmixs.f90 Fconf.f90 \
      Fmixpos.f90 Fmixneg.f90 F_vdW.f90\
      Fospi.f90 Fpol_sup.f90 \
      FmixHplus.f90 FmixOHmin.f90 \
      Fchem_eq.f90 \
      Fchem_eq_wall.f90 \
      pong_energy.f90

OBJS = $(SRC:.f90=.o)
##  OBJS = module_globales.o \
##        module_Csys.o \
##        module_FreeEnergy.o \
##        module_pore.o \
##        module_kinsolparam.o \
##        set_bulk_properties.o \
##        units_adaptation.o monolayer.o \
##        read_input.o \
##        savetodisk.o \
##        open_files.o \
##        close_files.o \
##        imprimir_cadenas.o \
##        save_data.o \
##        printstate.o\
##        check_run.o \
##        creador.o pxs.o kai.o \
##        cadenas72mr.o rota36.o\
##        mrrrr.o rands.o \
##        call_kinsol.o fkfun.o factorcurv.o \
##        fkpset.o fkpsol.o \
##        allocating.o \
##        set_initial_guess.o \
##        set_pore_distrib.o \
##        calc_conductance.o \
##        calc_mean_values.o \
##        calc_energy.f90 \
##        Fmixs.o Fconf.o \
##        Fmixpos.o Fmixneg.o \
##        FmixHplus.o FmixOHmin.o \
##        Fchem_eq.o \
##        pong_energy.o

objs_mpfun90 = mpfun90.o mpmod90.o mpmodm90.o mpmodx90.o

.SUFFIXES:            # this deletes the default suffixes 
.SUFFIXES: .f90 .o    # this defines the extensions I want
 
# some definitions

SHELL = /bin/bash

# This flags are used in the linker stage (LOAD LIBRARIES = LDFLAGS) 
# ***************
# Next FLAGS cames from: sundials-config -m kinsol -t s -l f
# FLAGS in Desktop ARG 
#LDFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.9.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.9.1/../../../x86_64-linux-gnu 
#LDFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.9.1/../../../../lib  -L/usr/lib/gcc/x86_64-linux-gnu/4.9.1/../../..
#LDFLAGS+= -lgfortran -lm -lgcc_s -lquadmath
 LDFLAGS= -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu 
 LDFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu
 LDFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lgfortran -lm -lgcc_s -lquadmath
 LDFLAGS+= -L/usr/local/share/lib 
 LDFLAGS+= -L/usr/local/lib
 LDFLAGS+= -L/lib/../lib -L/usr/lib/../lib 
 LDFLAGS+= -lsundials_fkinsol -lsundials_fnvecserial -lsundials_kinsol -lsundials_nvecserial -lm
# LDFLAGS+= -L/shared/software/sundials-2.5.0-openmpi/lib -lm 
# LDFLAGS+= -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib 

# LFLAGS=  -L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm 
# LFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.9.1/../../../x86_64-linux-gnu 
# LFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.9.1/../../../../lib 
# LFLAGS+= -L/lib/x86_64-linux-gnu 
# LFLAGS+= -L/lib/../lib 
# LFLAGS+= -L/usr/lib/x86_64-linux-gnu 
# LFLAGS+= -L/usr/lib/../lib
# LFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.9.1/../../.. -lgfortran -lm -lgcc_s -lquadmath
 LFLAGS=  -L/shared/software/sundials-2.5.0-openmpi/lib 
 LFLAGS+= -lsundials_fkinsol -lsundials_fnvecserial -lsundials_kinsol -lsundials_nvecserial -lm 
 LFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu 
 LFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib 
 LFLAGS+= -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu 
 LFLAGS+= -L/usr/lib/../lib 
 LFLAGS+= -L/usr/local/share/lib 
 LFLAGS+= -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath


# Actions
#all:   $(EXE) 
 
.f90.o:
	${FC} $(cflags) ${FFLAGS} -c $< -o $@ $(LDFLAGS)
#	${FC} $(cflags) -c ${FFLAGS} -c $< -o $@ $(LDFLAGS)
#	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(LDFLAGS)


$(EXE): $(OBJS) Makefile control_run.h 
#	$(FC) $(OBJS) -o $(EXE) $(LFLAGS) #$(LDFLAGS)
	$(FC) $(FFLAGS) $(OBJS) $(objs_mpfun90) -o $(EXE) $(LDFLAGS) $(LFLAGS) 
#	$(FF) $(FFLAGS) $(LFLAGS) $(LDFLAGS) $(OBJS) -o $(EXE) 

$(OBJS): control_run.h
monolayer: control_run.h Makefile control_run.h module_globales.o module_Csys.o control_run.h
module_globales.o: module_globales.f90 control_run.h
module_Csys.o: module_Csys.f90 control_run.h
fkfun.o: fkfun.f90 control_run.h
fkpsol.o: fkpsol.f90 control_run.h
fkpset.o: fkpset.f90 control_run.h
factorcurv.o: factorcurv.f90 control_run.h
set_bulk_properties.o : set_bulk_properties.f90 module_Csys.o control_run.h
kai.o : kai.f90 module_Csys.o control_run.h
cadenas72mr.o: cadenas72mr.f90 module_Csys.o control_run.h
creador.o: creador.f90 module_Csys.o control_run.h
mrrrr.o: mrrrr.f90 module_Csys.o control_run.h
pxs.o: pxs.f90 module_Csys.o control_run.h
rands.o: rands.f90 module_Csys.o control_run.h
read_input.o: read_input.f90 module_Csys.o control_run.h
rota36.o: rota36.f90 module_Csys.o control_run.h
units_adaptation.o: units_adaptation.f90 control_run.h
monolayer.o: module_globales.o module_Csys.o control_run.h Makefile control_run.h

#install: all
#	cp $(TARGET) ~/bin
 print:
	@echo "$(SRC)"
	@echo "$(OBJS)"

 clean:	
	@rm -rf csys.mod freeenergy.mod globales.mod pore.mod   $(SRC:.f90=.o) #$(SRC:.f90=.d) $(TARGET) *~

 realclean: clean
	@rm -f .depend

 depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 
 package:
	tar -czvf pore_progMP.tar.gz Makefile* *.f* monolayer fort.8 control_run.h NOTAS mp*
	@echo "***** Program packaged in pore_prog.tar.gz  *****"

ifeq (.depend, $(wildcard .depend))
include .depend
endif

