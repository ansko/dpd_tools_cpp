#ifndef STRUCTURE_ADD_MMT_PERIODIC_HPP
#define STRUCTURE_ADD_MMT_PERIODIC_HPP


#include "structures.hpp"


// Add an isolated MMT platelet
bool Structure::add_mmt_periodic(OptionsParser &o, float z, size_t charged_count)
{
    #ifdef DETAILED_OUTPUT  // Print parameters of MMT addition
        std::cout << "**********\n"
                  << "Structure.hpp add_mmt_periodic input parameters:\n"
                  << "platelet_edge = " << o.platelet_edge
                  << "\nbead_radius_clay = " << o.lj_bead_radius_clay
                  << "\natom_type = " << o.mmt_atom_type
                  << "\nmmt_edge_bond_type = " << o.mmt_edge_bond_type
                  << "\nmmt_diagonal_bond_type = " << o.mmt_diagonal_bond_type
                  << "\nz = " << z
                  << "\ncharged_count = " << charged_count
                  << "\nbead_charge = " << o.bead_charge
                  << "\n**********\n";
    #endif

    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    size_t new_idx = 1;

    // Add atoms and bonds top-bottom
    // bond length = reduction_k * 2 * o.lj_bead_radius_clay ...
    // but only along x and y
    float reduction_k(1.0);  // TODO: introduce !=1 values
    std::map<int, std::map<int, std::map<std::string, size_t> > > atom_ids;
    for (int idxx = 0; idxx < (int)o.platelet_edge; ++idxx)
      {
        atom_ids[idxx] = std::map<int, std::map<std::string, size_t> >();
        float dx = reduction_k * idxx * 2 * o.lj_bead_radius_clay
                   - (float)o.platelet_edge * o.lj_bead_radius_clay
                   + o.lj_bead_radius_clay;
        for (int idxy = 0; idxy < (int)o.platelet_edge; ++idxy)
          {
            atom_ids[idxx][idxy] = std::map<std::string, size_t>();
            float dy = reduction_k * idxy * 2 * o.lj_bead_radius_clay
                       - (float)o.platelet_edge * o.lj_bead_radius_clay
                       + o.lj_bead_radius_clay;
            new_atoms.push_back(Atom(dx, dy, z - o.lj_bead_radius_clay,
                0, o.mmt_atom_type, 0, 0, 0, "filler"));
            atom_ids[idxx][idxy]["bottom"] = new_idx++;
            new_atoms.push_back(Atom(dx, dy, z + o.lj_bead_radius_clay,
                0, o.mmt_atom_type, 0, 0, 0, "filler"));
            atom_ids[idxx][idxy]["top"] = new_idx++;
            // Only bonds top-bottom
            new_bonds.push_back(Bond(o.mmt_edge_bond_type, new_idx - 2,
                new_idx - 1));
          }
      }

    // Add other bonds
    for (auto idx_x = 0; idx_x < (int)o.platelet_edge; ++idx_x)
      {
        for (auto idx_y = 0; idx_y < (int)o.platelet_edge; ++idx_y)
          {
            // Edge bonds
            size_t other_idx_x_minus;
            size_t other_idx_y_minus;
            size_t other_idx_x_plus;
            size_t other_idx_y_plus;
            if (idx_x == 0)
              {
                other_idx_x_minus = o.platelet_edge - 1;
                other_idx_x_plus = 1;
              }
            else if (idx_x == o.platelet_edge - 1)
              {
                other_idx_x_minus = o.platelet_edge - 2;
                other_idx_x_plus = 0;
              }
            else
              {
                other_idx_x_minus = idx_x - 1;
                other_idx_x_plus = idx_x + 1;
              }
            if (idx_y == 0)
              {
                other_idx_y_minus = o.platelet_edge - 1;
                other_idx_y_plus = 1;
              }
            else if (idx_y == o.platelet_edge - 1)
              {
                other_idx_y_minus = o.platelet_edge - 2;
                other_idx_y_plus = 0;
              }
            else
              {
                other_idx_y_minus = idx_y - 1;
                other_idx_y_plus = idx_y + 1;
              }
            if (atom_ids.find(other_idx_x_minus) != atom_ids.end()
                && atom_ids[other_idx_x_minus].find(idx_y) !=
                   atom_ids[other_idx_x_minus].end())
              {
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[other_idx_x_minus][idx_y]["top"],
                    atom_ids[idx_x][idx_y]["top"]));
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[other_idx_x_minus][idx_y]["bottom"],
                    atom_ids[idx_x][idx_y]["bottom"]));
              }
            if (atom_ids[idx_x].find(other_idx_y_minus) !=
                atom_ids[idx_x].end())
              {
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[idx_x][other_idx_y_minus]["top"],
                    atom_ids[idx_x][idx_y]["top"]));
                new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                    atom_ids[idx_x][other_idx_y_minus]["bottom"],
                    atom_ids[idx_x][idx_y]["bottom"]));
              }

            // Diagonal bonds
            if (atom_ids.find(other_idx_x_plus) != atom_ids.end()
                && atom_ids[other_idx_x_plus].find(other_idx_y_plus) !=
                   atom_ids[other_idx_x_plus].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[other_idx_x_plus][other_idx_y_plus]["top"],
                    atom_ids[idx_x][idx_y]["bottom"]));
              }
            if (atom_ids.find(other_idx_x_plus) != atom_ids.end()
                && atom_ids[other_idx_x_plus].find(other_idx_y_minus) !=
                   atom_ids[other_idx_x_plus].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[other_idx_x_plus][other_idx_y_minus]["top"],
                    atom_ids[idx_x][idx_y]["bottom"]));
              }
            if (atom_ids.find(other_idx_x_minus) != atom_ids.end()
                && atom_ids[other_idx_x_minus].find(other_idx_y_plus) !=
                   atom_ids[other_idx_x_minus].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[other_idx_x_minus][other_idx_y_plus]["top"],
                    atom_ids[idx_x][idx_y]["bottom"]));
              }
            if (atom_ids.find(other_idx_x_minus) != atom_ids.end()
                && atom_ids[other_idx_x_minus].find(other_idx_y_minus) !=
                   atom_ids[other_idx_x_minus].end())
              {
                new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                    atom_ids[other_idx_x_minus][other_idx_y_minus]["top"],
                    atom_ids[idx_x][idx_y]["bottom"]));
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
        srand(time(NULL));
        std::set<int> charged_ids;
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

    // Append created structure to the existing structure
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


#endif  // STRUCTURE_ADD_MMT_PERIODIC_HPP
