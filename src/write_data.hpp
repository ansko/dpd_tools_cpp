#ifndef WRITE_DATA_HPP
#define WRITE_DATA_HPP


#include <fstream>
#include <set>
#include <string>

#include "structure.hpp"


void write_data(std::string out_fname, Structure &structure)
{
    std::ofstream ofs(out_fname);
    ofs << "LAMMPS data file by Anton\n\n";
    size_t atoms_count = structure.atoms().size();
    size_t bonds_count = structure.bonds().size();
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "Writing datafile " << out_fname << std::endl;
    std::cout << "Atoms count = " << atoms_count << std::endl;
    std::cout << "Bonds count = " << bonds_count << std::endl;
    #endif

    if (atoms_count > 0)
      {
        std::set<size_t> atom_types;
        for (auto & entry : structure.atoms())
          {
            atom_types.insert(entry.second.type);
          }
        ofs << atoms_count << " atoms\n";
        ofs << atom_types.size() << " atom types\n";
        #ifdef DEBUG
        std::cout << "Atom types count = " << atom_types.size() << std::endl;
        #endif
      }
    if (bonds_count > 0)
      {
        std::set<size_t> bond_types;
        for (auto & entry : structure.bonds())
          {
            bond_types.insert(entry.second.type);
          }
        ofs << bonds_count << " bonds\n";
        ofs << bond_types.size() << " bond types\n";
        #ifdef DEBUG
        std::cout << "Bond types count = " << bond_types.size() << std::endl;
        #endif
      }
    ofs << "\n";
    ofs << structure.xlo << " " << structure.xhi << " xlo xhi\n";
    ofs << structure.ylo << " " << structure.yhi << " ylo yhi\n";
    ofs << structure.zlo << " " << structure.zhi << " zlo zhi\n";
    ofs << "\n";
    if (atoms_count > 0)
      {
        ofs << "Atoms\n\n";
        for (auto & entry : structure.atoms())
          {
            ofs << entry.first << " " << 1
                << " " << entry.second.type
                << " " << entry.second.q
                << " " << entry.second.x
                << " " << entry.second.y
                << " " << entry.second.z
                << " " << entry.second.nx
                << " " << entry.second.ny
                << " " << entry.second.nz << "\n";
          }
      }
    if (bonds_count > 0)
      {
        ofs << "\nBonds\n\n";
        for (auto & entry : structure.bonds())
          {
            ofs << entry.first
                << " " << entry.second.type << " " << entry.second.atom_one
                << " " << entry.second.atom_two << "\n";
          }
      }
}

#endif  // WRITE_DATA_HPP
