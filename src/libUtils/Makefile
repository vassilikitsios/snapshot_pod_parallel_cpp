# -------------------------
BASE_DIR=../..
include $(BASE_DIR)/Makefile.in

LIBS = -lm

# --------Suffixes---------
.SUFFIXES: .cpp 

.cpp.o:
	$(CPP) $(CPPFLAGS) $(LIBS) -c $*.cpp

.cpp :
	$(CPP) $(CPPFLAGS) $(LIBS) -o $< $@

# -------Dependencies------
ARCHIVE_NAME = libUtils.a
HEADER_NAME = Utils.h

SRC = binary_read_or_write.cpp file_opening_closing.cpp memory_allocations.cpp misc.cpp param.cpp string.cpp time_stats.cpp tridiagonal_solver.cpp matrix.cpp cubic_spline.cpp finite_difference.cpp interpolation.cpp interpolation_1D.cpp

OBJS = $(addsuffix .o, $(basename $(SRC)))

HEADER_FILES = $(addsuffix .h, $(basename $(SRC)))

default: $(ARCHIVE_NAME)

$(ARCHIVE_NAME): $(OBJS)
	mkdir -p $(LIB_DIR) $(INCLUDE_DIR)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@
	cat $(HEADER_FILES) > $(INCLUDE_DIR)/$(HEADER_NAME)
	mv -v $(ARCHIVE_NAME) $(LIB_DIR)

# This automatically updates all the dependencies shown below
depend:
	makedepend -Y $(SRC)

clean:
	rm -fv *.o *~ ._* $(LIB_DIR)/$(ARCHIVE_NAME) $(INCLUDE_DIR)/$(HEADER_NAME)

# DO NOT DELETE

