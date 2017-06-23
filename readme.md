Multifluid, magnetospheric plasma simulations -- info for getting started.

1. Edit your .bashrc to contain the following lines:
	# Lifts restrictions on use of RAM
	ulimit -l unlimited
	# Adds compiler file locations to $PATH
	PATH=$PATH:/opt/intel/composer_xe_2013/bin:/usr/local/ncarg/bin:/opt/dx/bin
	# Identifies needed libraries to compiler
	export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013.2.146/compiler/lib/intel64
	# Optional: makes runscripts a little more convenient to use from anywhere
	alias run='$HOME/multifluid/scripts/run_multi.run'

2. Check that repository file path matches what the scripts expect. This path of this file should be: $HOME/multifluid/readme.md

3. Compile the executable by navigating to the execution directory ($HOME/multifluid) and typing the command: make
	This executes the makefile script to create an executable, by default named: master.x

4. Verify the settings of the input file. At the top, the 'start' variable should be set to .t. for an initial run.
	If running the executable directly, the input file should have the path $HOME/multifluid/input
	If you are queueing multiple runs, the first input file should have the path $HOME/multifluid/input_$runname_001, where $runname is the name of your executable (master by default).

5. Set the error message email to be your own. This is set by a variable at the top of both run_multi.run and archive.run, in the scripts directory: $HOME/multifluid/scripts

6. Start your run with the following syntax:
	For solo executable:
		cd ~/multifluid
		./master.x >> output.log ; echo "master.x finished running!" | mail -s "Multifluid run complete" your@email.com &
	For queued runs (note the use of the optional alias in .bashrc listed above):
		run $runname $start $stop &
	where $runname is the pre-extension executable name, $start and $stop are the indices of the first and last integers in the queued series.
	Example: run master 1 5 &
	makes 5 iterations of master.x.
	Input files for each iteration must be present in $execdir ($HOME/multifluid by default), named input_$runname_XXX, where XXX can range from 001 to 999.
	If the input file for a given iteration is not present, the input file for the previous iteration is copied over. Non-initializing runs must have the 'start' variable set to .f. in their input files.

7. Inspect the output data with the following command:
	idt gmeta
	OR
	idt $runname_XXX
	Example: idt master_001


Additional notes:
@ run_multi.run moves all data files into the data directory, $HOME/multifluid/data, once the executable has finished running.
	The output data file, named gmeta while the code is running, is renamed to: $runname_XXX

@ run_multi.run also depends on another script, archive.run. The latter compressed all relevant files from a single run into 1 tar archive in the completed directory: $HOME/multifluid/data/complete
	The gmeta file for each run is saved both inside and alongside the archive in the completed directory.
	All other files are moved to the trash directory: $HOME/multifluid/trash

@ Restarted and queued runs require 2 fluid data files to be present in addition to the correctly-named input file.
	Queued runs copy the files fluid13 and fluid14 from the previous iteration into fluid11 and fluid12 respectively.
	A restarted run, or the first of a set of queued runs, expects fluid11 and fluid12 to be present in the exection directory.
