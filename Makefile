# shaowei July 10/2018
.PHONY : all

# "-" at the beginning of each line asks make to ignore error message and not to stop
# "@" asks make not to print the commands that are being executed.
all :
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""
	@-cd ./exec-subtr/; make --keep-going
	@-cd ./exec-subtrx/; make --keep-going
	@-cd ./exec-fit/; make --keep-going
	#@-cd ./exec-TMVA/; make --keep-going
	@-echo ""
	@-echo "****////   NOTE: Run make repeatedly until 'Nothing to be done for' every directory.   ////****"
	@-echo ""

.PHONY : clean
clean :
	rm -f ./obj/*/*.o
	rm -f ./bin/*
