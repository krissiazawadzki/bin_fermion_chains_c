CCOMPILER = gcc
CXXCOMPILER = g++

DIR_UTIL_ROUTINES = ../../cpp_util_modules/util_etc/
DIR_SRC_LIBS = ../src

INC_BLAS = "/usr/local/opt/openblas/include/"
LAPACKFLAGS += -DWITH_LAPACK -llapack -lblas -lm -lc  -lstdc++  -framework Accelerate

LDFLAGS = -L/usr/lib64 -L/usr/lib/ -L/home/krissia/lib/ -lblas -llapack -lm -lc  -lstdc++  -framework Accelerate

CFLAGS += -Wall -O2   -I$(INC_BLAS) -I$(DIR_UTIL_ROUTINES) -I$(DIR_SRC_LIBS) #-std=c++11
CXXFLAGS += -Wall -O2   -I$(INC_BLAS) -I$(DIR_UTIL_ROUTINES) -I$(DIR_SRC_LIBS) #-std=c++11



MAINSOURCE = testing_new_module


OBJECTS =	$(DIR_UTIL_ROUTINES)/matrix_ops/matrix_ops.o \
			$(DIR_UTIL_ROUTINES)/file_ops/file_ops.o \
			$(DIR_UTIL_ROUTINES)/file_ops/mkdir.o \
			$(DIR_UTIL_ROUTINES)/string_ops/string_ops.o \
			$(DIR_SRC_LIBS)/spin_configurations.o \
			$(DIR_SRC_LIBS)/qs_hilbert.o	\
			$(DIR_SRC_LIBS)/bin_sectors.o	\
			$(DIR_SRC_LIBS)/QSmatrix.o	\
			$(DIR_SRC_LIBS)/first_iter.o \
			$(DIR_SRC_LIBS)/prim.o \
			$(DIR_SRC_LIBS)/io_bin.o \
			$(DIR_SRC_LIBS)/translation_eigenstates_of_S.o \
			$(DIR_SRC_LIBS)/int_array_util.o


DEPDIR = .deps
df = $(DEPDIR)/$(*F)



# test new module
testing_new_module: $(DEPDIR) $(OBJECTS) testing_new_module.o makefile
	$(CXX) $(CXXFLAGS) $(OBJECTS) testing_new_module.o -o testing_new_module.x $(LAPACKFLAGS) 
	
testing_new_module-debug:  CFLAGS += -g 	
testing_new_module-debug:  CXXFLAGS += -g 	
testing_new_module-debug:  testing_new_module


reduced_hilbert_QS: $(DEPDIR) $(OBJECTS) reduced_hilbert_QS.o makefile
	$(CXX) $(CXXFLAGS) $(OBJECTS) reduced_hilbert_QS.o -o reduced_hilbert_QS.x $(LAPACKFLAGS) $(LDFLAGS)
	
reduced_hilbert_QS-debug:  CFLAGS += -g 	
reduced_hilbert_QS-debug:  CXXFLAGS += -g 	
reduced_hilbert_QS-debug:  reduced_hilbert_QS


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
