#ifndef STRUCTURE_ADD_MMT_CIRCULAR_HPP
#define STRUCTURE_ADD_MMT_CIRCULAR_HPP


#include "structures.hpp"


// Add an "infinite" (because of periodic boundary conditions) MMT platelet
// into existing Structure (==*this)
bool Structure::add_mmt_circular(OptionsParser &o, float x, float y, float z,
size_t charged_count)
{
    // Print (or not) passed parameters of MMT addition
    #ifdef DETAILED_OUTPUT
        std::cout << "**********\n"
                  << "\nStructure.hpp add_mmt_circular input parameters:\n"
                  << "\nplatelet_radius = " << o.platelet_radius
                  << "\nbead_radius_clay = " << o.lj_bead_radius_clay
                  << "\natom_type = " << o.mmt_atom_type
                  << "\nmmt_edge_bond_type = " << o.mmt_edge_bond_type
                  << "\nmmt_diagonal_bond_type = " << o.mmt_diagonal_bond_type
                  << "\nx = " << x
                  << "\ny = " << y
                  << "\nz = " << z
                  << "\ncharged_count = " << charged_count
                  << "\nbead_charge = " << o.bead_charge
                  << "\n**********\n";
    #endif

    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    size_t new_idx = 1;

    // Add atoms and bonds between top atom and neighboring bottom atom
    std::map<int, std::map<int, std::map<std::string, size_t> > > atom_ids;
    for (int idxx = -(int)o.platelet_radius + 1; idxx < (int)o.platelet_radius;
         ++idxx)
      {
        atom_ids[idxx] = std::map<int, std::map<std::string, size_t> >();
        float dx = idxx * 2 * o.lj_bead_radius_clay;
        for (int idxy = -(int)o.platelet_radius + 1;
             idxy < (int)o.platelet_radius; ++idxy)
          {
            if (pow(o.platelet_radius, 2) < pow(idxx, 2) + pow(idxy, 2))
              {
                continue;
              }
            atom_ids[idxx][idxy] = std::map<std::string, size_t>();
            float dy = idxy * 2 * o.lj_bead_radius_clay;
            new_atoms.push_back(Atom(x + dx, y + dy, z - o.lj_bead_radius_clay,
                0, o.mmt_atom_type, 0, 0, 0, "filler"));
            atom_ids[idxx][idxy]["bottom"] = new_idx;
            new_atoms.push_back(Atom(x + dx, y + dy, z + o.lj_bead_radius_clay,
                0, o.mmt_atom_type, 0, 0, 0, "filler"));
            atom_ids[idxx][idxy]["top"] = new_idx + 1;
            new_bonds.push_back(Bond(o.mmt_edge_bond_type, new_idx,
                new_idx + 1));
            new_idx += 2;
          }
      }

    // Add other bonds
    for (auto itx = atom_ids.begin(); itx != atom_ids.end(); ++itx)
      {
        for (auto ity = atom_ids[itx->first].begin();
             ity != atom_ids[itx->first].end(); ++ity)
          {
            // Edge bonds
            if (atom_ids.find(itx->first - 1) != atom_ids.end()
                && atom_ids[itx->first - 1].find(ity->first) !=
                   atom_ids[itx->first - 1].end())
              {
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[itx->first - 1][ity->first]["top"],
                    atom_ids[itx->first][ity->first]["top"]));
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[itx->first - 1][ity->first]["bottom"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids[itx->first].find(ity->first - 1) !=
                atom_ids[itx->first].end())
              {
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[itx->first][ity->first - 1]["top"],
                    atom_ids[itx->first][ity->first]["top"]));
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[itx->first][ity->first - 1]["bottom"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            // Diagonal bonds
            if (atom_ids.find(itx->first - 1) != atom_ids.end()
                && atom_ids[itx->first - 1].find(ity->first - 1) !=
                   atom_ids[itx->first - 1].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[itx->first - 1][ity->first - 1]["top"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids.find(itx->first - 1) != atom_ids.end()
                && atom_ids[itx->first - 1].find(ity->first + 1) !=
                   atom_ids[itx->first - 1].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[itx->first - 1][ity->first + 1]["top"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids.find(itx->first + 1) != atom_ids.end()
                && atom_ids[itx->first + 1].find(ity->first - 1) !=
                   atom_ids[itx->first + 1].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[itx->first + 1][ity->first - 1]["top"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids.find(itx->first + 1) != atom_ids.end()
                && atom_ids[itx->first + 1].find(ity->first + 1) !=
                   atom_ids[itx->first + 1].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[itx->first + 1][ity->first + 1]["top"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
          }
      }

    // Assign (or not) charges
    if (charged_count != 0)
      {
        if (new_atoms.size() < charged_count)
          {
            std::cerr << "Cannot charge " << charged_count
                      << " atoms from " << new_atoms.size() <<std::endl;
          }
        std::set<int> charged_ids;
        srand(time(NULL));
        while (charged_ids.size() < charged_count)
          {
            charged_ids.insert(rand() % new_atoms.size());
          }
        for (auto it = charged_ids.begin(); it != charged_ids.end(); ++it)
          {
            new_atoms[*it].q = -o.bead_charge;
          }
      }

    // Results of MMT addition
    #ifdef DETAILED_OUTPUT
        std::cout << "**********\n";
        std::cout << "MMT addition in structure.hpp:\n";
        std::cout << "atoms_count = " << new_atoms.size() << std::endl;
        std::cout << "bonds_count = " << new_bonds.size() << std::endl;
        std::cout << "**********\n";
    #endif

    // Append created structure into the existing structure
    for (size_t idx = 0; idx < new_bonds.size(); ++idx)
      {
        new_bonds[idx].atom_one += this->_atoms_count;
        new_bonds[idx].atom_two += this->_atoms_count;
        this->_bonds[this->_bonds_count + 1 + idx] = new_bonds[idx];
      }
    this->_bonds_count += new_bonds.size();
    for (size_t idx = 0; idx < new_atoms.size(); ++idx)
      {
        this->_atoms[this->_atoms_count + 1 + idx] = new_atoms[idx];
      }
    this->_atoms_count += new_atoms.size();
    return true;
}


#endif  // STRUCTURE_ADD_MMT_CIRCULAR_HPP
