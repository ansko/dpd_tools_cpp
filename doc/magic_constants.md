# Magic constants

Some comments may be missed in source code (or placed in separate files).


### From the experiment or literature

CEC = 92.6 mequiv/100g = 92.6 * 1e-8 equiv/g
Cation exchange capacity, mequiv per 100 g. The value 92.6 is taken from 
literature (Anoukou et al.). The value is used below:

dpd_rho  = 3
DPD density. Is taken from literature.


### From our MD calculations

lx = 90
ly = 80
lz = 65
Size (approximate but close to the actual) of the L system from the MD
calculations. It is in a good agreement with reality (density of polymer differs).

exchange_surface_density = 0.15 beacuse:
cell_square = 400 (per one cell)
cell_mass = 720 (per one cell)
number of cations to be exchanged in the system equal to 162 cells:
92.6 * 1e-8 * 1e3 * (720 * 162) = 108  (1e3: g<->kg)
exchanges surface densuty:
exchange_surface_density = 108 excahnges / lx / ly = 108 / 80 / 90 = 0.15

mmt_real_thickness = 8
Approximate value. Real thickness is less, but mmt is surrounded by some 
emptiness that divides it from modifier or polymer.

md_soft_atoms = 41940
Number of soft (not mmt; modifier + polymer) atoms in the L system. Since 
modifier's and polymer's atoms are treated as having equal radii, here it is.
