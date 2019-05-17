## The way masses are being computed

### Mass values in forcefields (Clayff and CVFF):

Al  26.98154   (ao)
Mg  24.305     (mgo)
Si  28.0855    (st)
O   15.9994    (ob, obos, oh, ohs; o') 
H   1.00797    (ho; h, hn)
N   14.006700  (n4, n2, n)
C   12.011150  (c2, c3, c')


### Some masses in real units

MMT mass in mix.data during MD:
720 atoms that are
18 cells x 40 atoms per cell
18 x (ao x 4 + st x 8 + o* x 24 + h x 4) + 12 (mgo - ao)

18 * (4*26.98154 + 8*28.0855 + 24*15.9994 + 4*1.00797) + 12 * (24.305 - 26.98154)
= 12939.179039999999

lx = 31.08
ly = 26.94

dens = 12939.179039999999 / 31.08 / 26.94 = 15.453544986284404 amu / AA


### Some masses in LJ units (periodic)

lx = ly = xy_size * real_bead_radius


### Into real

mmt mass = dens * lx * ly (amu, dens == 15.45)
bead mass = mmt mass / mmt_atoms

= 66 for lj_bead_radius == 1


### Monomer mass

6*c* + n4 + o' + 13*h*
= 6*12.01115 + 14.0067 + 15.999400 + 13*1.00797 = 115
(1.74 * 66)


### Modifier beads

2 * (2*c* + o* + 5*h*) + c* + 3*h* + n* + 16*c* + 33*h*
= 21*c* + 2*o* + n* + 46*h*
= 21*12.01115 + 2*15.9994 + 14.0067 + 46*1.00797 = 344 (per 3 beads)
1 bead = 115
(1.74 * 66)
