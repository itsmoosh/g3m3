The Fortran code in this repository is split into the following categories:

3d:		Primary program sequence, including the main time-stepping loop.
cft:	Spacecraft subroutines:	reading, writing, and calculations associated
		with spacecraft inputs and outputs
fnc:	Functions: In-line functions.
gfx:	Graphics subroutines:	mid-simulation graphics outputs, which are
		stored in the gmeta file.
mod:	Modules containing global constants.
sub:	Standard (general) subroutines.
cut:	Stand-alone code for parsing binary fluid files into ASCII output data.
		Unlike the other categories, this category is not in the code/
		directory.
