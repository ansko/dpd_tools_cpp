#include "structures.hpp"


// Add an "infinite" (because of periodic boundary conditions) MMT platelet
// into existing Structure
bool
Structure::add_mmt_circular(AddMmtCircularParameters &parameters)
{
    // Unpacking
    float x = parameters.x;
    float y = parameters.y;
    float z = parameters.z;
    size_t charged_count = parameters.charged_count;
    float lj_bead_radius_clay(parameters.lj_bead_radius_clay);
    float bead_charge(parameters.bead_charge);
    float platelet_closing(parameters.platelet_closing);
    size_t platelet_radius(parameters.platelet_radius);
    size_t mmt_atom_type(parameters.mmt_atom_type);
    size_t mmt_edge_bond_type(parameters.mmt_edge_bond_type);
    size_t mmt_diagonal_bond_type(parameters.mmt_diagonal_bond_type);

    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    size_t new_idx = 1;

    // Add atoms and bonds between top atom and neighboring bottom atom
    std::map<int, std::map<int, std::map<std::string, size_t> > > atom_ids;
    for (int idxx = -(int)platelet_radius + 1; idxx < (int)platelet_radius;
         ++idxx)
      {
        atom_ids[idxx] = std::map<int, std::map<std::string, size_t> >();
        float dx = idxx * 2 * lj_bead_radius_clay * platelet_closing;
        for (int idxy = -(int)platelet_radius + 1;
             idxy < (int)platelet_radius; ++idxy)
          {
            if (pow(platelet_radius, 2) < pow(idxx, 2) + pow(idxy, 2))
              {
                continue;
              }
            atom_ids[idxx][idxy] = std::map<std::string, size_t>();
            float dy = idxy * 2 * lj_bead_radius_clay * platelet_closing;
            new_atoms.push_back(Atom(x + dx, y + dy,
                z - lj_bead_radius_clay*platelet_closing,
                0, mmt_atom_type, 0, 0, 0, "filler"));
            atom_ids[idxx][idxy]["bottom"] = new_idx;
            new_atoms.push_back(Atom(x + dx, y + dy,
                z + lj_bead_radius_clay*platelet_closing,
                0, mmt_atom_type, 0, 0, 0, "filler"));
            atom_ids[idxx][idxy]["top"] = new_idx + 1;
            new_bonds.push_back(Bond(mmt_edge_bond_type, new_idx,
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
                new_bonds.push_back(Bond(mmt_edge_bond_type,
                    atom_ids[itx->first - 1][ity->first]["top"],
                    atom_ids[itx->first][ity->first]["top"]));
                new_bonds.push_back(Bond(mmt_edge_bond_type,
                    atom_ids[itx->first - 1][ity->first]["bottom"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids[itx->first].find(ity->first - 1) !=
                atom_ids[itx->first].end())
              {
                new_bonds.push_back(Bond(mmt_edge_bond_type,
                    atom_ids[itx->first][ity->first - 1]["top"],
                    atom_ids[itx->first][ity->first]["top"]));
                new_bonds.push_back(Bond(mmt_edge_bond_type,
                    atom_ids[itx->first][ity->first - 1]["bottom"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            // Diagonal bonds
            if (atom_ids.find(itx->first - 1) != atom_ids.end()
                && atom_ids[itx->first - 1].find(ity->first - 1) !=
                   atom_ids[itx->first - 1].end())
              {
                new_bonds.push_back(Bond(mmt_diagonal_bond_type,
                    atom_ids[itx->first - 1][ity->first - 1]["top"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids.find(itx->first - 1) != atom_ids.end()
                && atom_ids[itx->first - 1].find(ity->first + 1) !=
                   atom_ids[itx->first - 1].end())
              {
                new_bonds.push_back(Bond(mmt_diagonal_bond_type,
                    atom_ids[itx->first - 1][ity->first + 1]["top"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids.find(itx->first + 1) != atom_ids.end()
                && atom_ids[itx->first + 1].find(ity->first - 1) !=
                   atom_ids[itx->first + 1].end())
              {
                new_bonds.push_back(Bond(mmt_diagonal_bond_type,
                    atom_ids[itx->first + 1][ity->first - 1]["top"],
                    atom_ids[itx->first][ity->first]["bottom"]));
              }
            if (atom_ids.find(itx->first + 1) != atom_ids.end()
                && atom_ids[itx->first + 1].find(ity->first + 1) !=
                   atom_ids[itx->first + 1].end())
              {
                new_bonds.push_back(Bond(mmt_diagonal_bond_type,
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
            new_atoms[*it].q = -bead_charge;
          }
      }

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
