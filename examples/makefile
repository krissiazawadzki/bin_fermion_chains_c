CCOMPILER = gcc
CXXCOMPILER = g++

DIR_UTIL_ROUTINES = ../../dft-nrg/src/util_etc/
DIR_SRC_LIBS = ../src

INC_BLAS = "/usr/local/opt/openblas/include/"
LAPACKFLAGS += -llapack -lblas -llapacke


CFLAGS += -Wall -O2   -I$(INC_BLAS) -I$(DIR_UTIL_ROUTINES) -I$(DIR_SRC_LIBS)
CXXFLAGS += -Wall -O2   -I$(INC_BLAS) -I$(DIR_UTIL_ROUTINES) -I$(DIR_SRC_LIBS) -std=c++11



MAINSOURCE = testing_new_module


OBJECTS =	$(DIR_UTIL_ROUTINES)/matrix_ops/matrix_ops.o \
			$(DIR_SRC_LIBS)/spin_configurations.o \
			$(DIR_SRC_LIBS)/qs_hilbert.o	\
			$(DIR_SRC_LIBS)/bin_sectors.o	


DEPDIR = .deps
df = $(DEPDIR)/$(*F)



# LSOC
testing_new_module: $(DEPDIR) $(OBJECTS) testing_new_module.o makefile
	$(CXX) $(CXXFLAGS) $(OBJECTS) testing_new_module.o -o testing_new_module.x $(LAPACKFLAGS) 
	
testing_new_module-debug:  CFLAGS += -g 	
testing_new_module-debug:  CXXFLAGS += -g 	
testing_new_module-debug:  testing_new_module



clean:
	rm -f $(OBJECTS) *.o
distclean: clean
	rm -rf $(DEPDIR)

%.o: %.c makefile
	$(COMPILE.c) -MMD -o $@ $<
	@cp $*.d $(df).P; \
		sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		    -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
		rm -f $*.d
%.o: %.cpp makefile
	$(COMPILE.cpp) -MMD -o $@ $<
	@cp $*.d $(df).P; \
		sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		    -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
		rm -f $*.d
$(DEPDIR):
	mkdir -p $(DEPDIR)
-include $(OBJECTS:%.o=$(DEPDIR)/%.P)
