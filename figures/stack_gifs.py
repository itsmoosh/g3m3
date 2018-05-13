"""
Compiles gifs of all .pngs from 001 to the number given.
Required arguments:
	run_name
	nplots (I0.3 format)

Example from terminal: python3 ./figures/stack_gifs.py 
"""

import gfx_functions as gfx
import sys

qty_list = ('bandalf', 'flows', 'qpres', 'hpres', 'opres', 'elec', 'temps', 'dens')
qty_list = ('temps',)

gfx.upd_gifs(qty_list, sys.argv[1], sys.argv[2], n_grids=5, n_skipped=int(sys.argv[2]))
