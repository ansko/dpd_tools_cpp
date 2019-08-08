# DPD tools


https://github.com/ansko/dpd_tools_cpp

Some tools for structure creation for DPD calculations, written in C++ for the
the faster structure composition.


## source code

### isolated_parallel.cpp

Create isolated mmt platelet(s) surrounded by modifier in polymer. Parallelism 
is done as split of the parallelepiped cell into sub-parallelepipeds along axes,
their count is nx*ny*nz. Though the program works even when nx, ny and nz are 
notably anisotropic, it is not a good case beacause may lead to a strong 
orientation of polymer chains.

On my notebook complete threads count:
4 is ok for 41k atoms (74 to 14 s)
8 is ok for 100k atoms (403 to 85 s)


### periodic.cpp

Create periodic structure made of mmt, modifier and polymer, similar to those
studied by means of MD. Mainly is done for the parameterization adjusting.


## Some useful links


https://github.com/lammps/lammps/blob/master/tools/msi2lmp/frc_files/clayff.frc
