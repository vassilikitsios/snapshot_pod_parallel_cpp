# Vassili Kitsios 1/12/2006
# -------------------------
BASE_DIR = ../..
include $(BASE_DIR)/Makefile.in

LIBS = -lm -L$(LIB_DIR)/ $(POD_LIBS) $(UTILS_LIBS) $(VTK_LIBS) $(CLAPACK_LIBS)
INCL = -I$(INCLUDE_DIR)/ $(UTILS_INC) $(VTK_INC) $(CLAPACK_INC)

# --------Suffixes---------
.SUFFIXES: .cpp 

.cpp.o:
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -c $*.cpp

.cpp :
	$(CPP) $(CPPFLAGS) $(LIBS) $(INCL) -o $< $@

# -------Dependencies------
BIN_NAME = snapshot_pod
INPUT_DECK_NAME = $(BIN_NAME).in
RUN_NAME = $(BIN_NAME).run

SRC = input_data.cpp main.cpp

OBJS = $(addsuffix .o, $(basename $(SRC))) 

# -------Make Commands------
default: $(BIN_NAME)
	(mkdir -p $(BASE_DIR)/bin ; mv $(BIN_NAME) ./bin ; cp ./bin/$(BIN_NAME)* $(BASE_DIR)/bin)

$(BIN_NAME): $(OBJS)
	$(CPP) $(CPPFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -fv ./bin/$(BIN_NAME) $(BASE_DIR)/bin/$(BIN_NAME)* ~/bin/$(BIN_NAME)* *.o *~ ._*  ./bin/*~

# This automatically updates all the dependencies shown below
depend:
	makedepend -Y $(SRC)
	
