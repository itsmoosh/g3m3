#!/bin/bash
# Combines data files output by code/3d_write_graphing_data into a compressed tar archive.
# Required inputs:
#	path/to/data_dir
#	run_name
#	# in seq (I0.3 format)

# Save the current directory to go back to after we're done
pwdvar="$(pwd)"
# Change to the data dir so we don't have annoying folder in our archive
cd $1
# Archive the data
tar -czf ${2}_t${3}_data.tar.gz ${2}_*_t${3}.dat
# Get rid of the (now compressed) data files
rm ${2}_*_t${3}.dat
# Go back to our original directory
cd $pwdvar
