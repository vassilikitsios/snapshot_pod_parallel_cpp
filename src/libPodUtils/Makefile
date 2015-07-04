# Vassili Kitsios 27/02/2007
# -------------------------
BASE_DIR = ../..
include $(BASE_DIR)/Makefile.in

LIBS = -lm -L$(LIB_DIR)/ $(UTILS_LIBS) $(CLAPACK_LIBS) $(VTK_LIBS)
INCL = -I$(INCLUDE_DIR)/ $(CLAPACK_INC) $(VTK_INC)

# --------Suffixes---------
.SUFFIXES: .cpp 

.cpp.o:
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -c $*.cpp

.cpp :
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -o $< $@

# -------Dependencies------
ARCHIVE_NAME = libPodUtils.a
HEADER_NAME = PodUtils.h

SRC = Snapshot_Data_t.cpp Snapshot_Derivative_Data_t.cpp Inner_Product_Calculator_t.cpp Spatial_Mode_Calculator_t.cpp Field_Reconstructor_t.cpp real_sym_eig_solver.cpp POD_t.cpp

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
