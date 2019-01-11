This repository contains software written mostly in Fortran for simulating magnetospheric plasmas. A multifluid model is used that treats electrons, proton solar wind, hydrogen, and oxygen as resistive, charged fluids. Requires the Intel Fortran compiler ifort and OpenMP for parallelization.

There are scripts serving 5 purposes:
	1. To compile the code (see step 4 below)
	2. To run the simulations (step 7 below)
	3. To print images of physical quantities (step 3 below)
	4. To archive binary output data and logs (step 9 below)
	5. To parse fluid files into ASCII from their default binary format

Info for getting started:

	1. Edit your .bashrc to contain the following lines:
		# Lifts restrictions on use of RAM
		ulimit -l unlimited
		# Adds compiler file locations to $PATH
		PATH=$PATH:/opt/intel/composer_xe_2013/bin:/usr/local/ncarg/bin:/opt/dx/bin
		# Identifies needed libraries to compiler
		export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013.2.146/compiler/lib/intel64
		# Optional: makes runscripts a little more convenient to use from anywhere
		alias run='$HOME/g3m3/scripts/run_multi.run'

	2. Check that repository file path matches what the scripts expect. This path of this file should be: $HOME/g3m3/readme.md

	3. Install miniconda3 with the modules below and ensure that your $PATH points to the miniconda3 version of python3. This is the default during installation. Use the following commands to get python3 ready:
		@ cd ~
		@ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
		@ chmod +x Miniconda3-latest-Linux-x86_64.sh
		@ ./Miniconda3-latest-Linux-x86_64.sh
		@ Accept default prompts, except for recommended install location: ~/usr/local/miniconda3
		@ conda install matplotlib numpy scipy jupyter pandas
		@ Find the following line in the file equivalent to ~/usr/local/miniconda3/pkgs/matplotlib-2.2.2-py36h0e671d2_1/lib/python3.6/site-packages/matplotlib/streamplot.py (with approximate line numbers):
			390 def _update_trajectory(self, xm, ym):
			391         """Update current trajectory position in mask.
			392 
			393         If the new position has already been filled, raise `InvalidIndexError`.
			394         """
			395         if self._current_xy != (xm, ym):
			396             if self[ym, xm] == 0:
			397                 self._traj.append((ym, xm))
			398                 self._mask[ym, xm] = 1
			399                 self._current_xy = (xm, ym)
			400             else:
			401                 pass # <----- Add this line
			402                 #raise InvalidIndexError # <----- Comment out this line

	4. Compile the executable by navigating to the execution directory ($HOME/g3m3) and typing the command: make ; make gifs
		This executes the makefile script to create an executable, by default named: master.x
		make gifs compiles an executable for parallelizing gif stacking.

	5. Verify the settings of the input file. At the top, the 'start' variable should be set to .t. for an initial run.
		If running the executable directly, the input file should have the path $HOME/g3m3/input
		If you are queueing multiple runs, the first input file should have the path $HOME/g3m3/input_$runname_001, where $runname is the name of your executable (master by default).

	6. Set the error message email to be your own. This is set by a variable at the top of both run_multi.run and archive.run, in the scripts directory: $HOME/g3m3/scripts

	7. Start your run with the following syntax:
		For solo executable:
			cd ~/g3m3
			./master.x >> output.log 2>&1 ; echo "master.x finished running!" | mail -s "Multifluid run complete" your@email.com &
		For queued runs (note the use of the optional alias in .bashrc listed above):
			run $runname $start $stop &
		where $runname is the pre-extension executable name, $start and $stop are the indices of the first and last integers in the queued series.
		Example: run master 1 5 &
		makes 5 iterations of master.x.
		Input files for each iteration must be present in $execdir ($HOME/g3m3 by default), named input_$runname_XXX, where XXX can range from 001 to 999.
		If the input file for a given iteration is not present, the input file for the previous iteration is copied over. Non-initializing runs must have the 'start' variable set to .f. in their input files.

	8. Inspect the output data images in the figures/images directory, e.g. using evince.
	
	9. Archive your data with the following command:
		From the scripts directory:
			./archive.run $runname $runnum
		where $runname is as above, and runnum matches one of the start or stop numbers.
		Example: ./archive.run master 1
	This will compress the 4 fluid files and the associated input file, output log, etc. into a single .tar.gz archive in the g3m3/data/complete directory. The fluid files and associated files will be moved to the g3m3/trash directory. This directory should be cleared periodically. Files are not directly deleted because if there are missing files, the tar command may not work correctly and files would then be lost.


Additional notes:
@ run_multi.run moves all data files into the data directory, $HOME/g3m3/data, once the executable has finished running.
	The output data file, named gmeta while the code is running, is renamed to: $runname_XXX

@ run_multi.run also depends on another script, archive.run. The latter compressed all relevant files from a single run into 1 tar archive in the completed directory: $HOME/g3m3/data/complete
	The gmeta file for each run is saved both inside and alongside the archive in the completed directory.
	All other files are moved to the trash directory: $HOME/g3m3/trash

@ Restarted and queued runs require 2 fluid data files to be present in addition to the correctly-named input file.
	Queued runs copy the files fluid13 and fluid14 from the previous iteration into fluid11 and fluid12 respectively.
	A restarted run, or the first of a set of queued runs, expects fluid11 and fluid12 to be present in the exection directory.
	Subsequent runs must have 'input_out' copied to 'input' for ut values to work properly.
