#include "structures.hpp"


// Add an "infinite" (because of periodic boundary conditions) MMT platelet
// into existing Structure
bool
Structure::add_mmt_circular_fcc(AddMmtCircularParameters &parameters)
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
    int idx_min_x((-float(platelet_radius) - 1));
    int idx_max_x(platelet_radius + 1);
    int idx_min_y((-float(platelet_radius) - 1));
    int idx_max_y(platelet_radius + 1);

    float R(platelet_radius * 2 * lj_bead_radius_clay * platelet_closing);
    float bond_len(2 * lj_bead_radius_clay * platelet_closing);
    float bond_len2(bond_len * bond_len);
    float tet_height = bond_len * sqrt(2.0 / 3.0);  // half of thickness

    for (int idxx = idx_min_x; idxx < idx_max_x; ++idxx)
      {
        for (int idxy = idx_min_y; idxy < idx_max_y; ++idxy)
          {
            float dx = idxx * bond_len + idxy * bond_len * cos(M_PI/3);
            float dy = idxy * bond_len * sin(M_PI/3);

            if (sqrt(dx*dx + dy*dy) > R)
              {
                continue;
              }

            new_atoms.push_back(Atom(x + dx, y + dy, z - tet_height,
                                     0, mmt_atom_type, 0, 0, 0, "filler"));
            new_atoms.push_back(Atom(x + dx, y + dy, z + tet_height,
                                     0, mmt_atom_type, 0, 0, 0, "filler"));

            float dx_small (2.0/3.0 * sqrt(3) / 2 * cos(M_PI/6));
            float dy_small (2.0/3.0 * sqrt(3) / 2 * sin(M_PI/6));
            if (sqrt((dx + dx_small)*(dx + dx_small)
                      + (dy + dy_small)*(dy + dy_small)) > R)
              {
                continue;
              }


            new_atoms.push_back(Atom(x + dx + dx_small, y + dy + dy_small, z,
                                     0, mmt_atom_type, 0, 0, 0, "filler"));
            new_idx += 3;
          }
      }

    std::cout << new_idx << std::endl;

    for (size_t aidx1 = 0; aidx1 < new_atoms.size(); ++aidx1)
        for (size_t aidx2 = aidx1 + 1; aidx2 < new_atoms.size(); ++aidx2)
          {
            float dx(std::fabs(new_atoms[aidx1].x - new_atoms[aidx2].x));
            float dy(std::fabs(new_atoms[aidx1].y - new_atoms[aidx2].y));
            float dz(std::fabs(new_atoms[aidx1].z - new_atoms[aidx2].z));
            float dr2(dx*dx + dy*dy + dz*dz);
            if (dr2 < 1.05 * bond_len2)
              {
                bool found(false);
                for (Bond &b : new_bonds)
                  {
                    if((b.atom_one == aidx1 && b.atom_two == aidx2)
                        || (b.atom_two == aidx1 && b.atom_one == aidx2))
                      {
                        found = true;
                        continue;
                      }
                  }
                if (!found)
                  {
                    new_bonds.push_back(Bond(1, aidx1, aidx2));
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
        new_bonds[idx].atom_one += this->_atoms_count + 1;
        new_bonds[idx].atom_two += this->_atoms_count + 1;
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
