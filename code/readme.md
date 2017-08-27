The Fortran code in this repository is split into the following categories:

3d:	Primary program sequences, including the main time-stepping loop.
a:	Anisotropy subroutines:	associated with pressure anisotropy in the plasma, resulting from helical motion along magnetic field lines.
c:	Spacecraft subroutines:	reading, writing, and calculations associated with spacecraft inputs and outputs
g:	Graphics subroutines:	mid-simulation graphics outputs, which are stored in the gmeta file.
s:	Standard (general) subroutines.
cut:	Stand-alone code for parsing binary fluid files into ASCII output data. Unlike the other categories, this category is not in the code/ directory.

Each category of subroutines was contained in a single large file until 31-07-2017. The commit history for each large file is now associated with the file listed first alphabetically, except for 3d because 3d_main,f90 was retained.
