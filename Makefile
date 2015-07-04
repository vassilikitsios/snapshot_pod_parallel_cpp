#-------------------------------------------
# To add in new system options:
#               1) copy the ./makefiles/Makefile.generic.in file to ./makefiles/Makefile.YOUR_SYSTEM.in
#               2) make required changes to ./makefiles/Makefile.YOUR_SYSTEM.in
#               3) in the present Makefile modify the default print statement to includei "YOUR_SYSTEM"
#               4) in the present Makefile add "YOUR_SYSTEM" to the line currently containing "generic" and "cygwin"
#-------------------------------------------

default:
	@echo "-------------------------------------------"
	@echo "The system options are:"
	@echo "  make MacBook"
	@echo "  make edda"
	@echo "  make tango"
	@echo "  make dell"
	@echo "  make all (uses current architecture)"
	@echo "-------------------------------------------"

all:
	(cd ./src ; make ; cd ../drivers ; make)

MacBook edda tango alienor dell:
	(rm -f Makefile.in ; ln -s makefiles/Makefile.$@.in Makefile.in ; make clean ; make all)

clean:
	(rm -fv *~ ._* ; cd ./src ; make clean ; cd ../drivers ; make clean ; cd ../makefiles ; rm -fv *~)

