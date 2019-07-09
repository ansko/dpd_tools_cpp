#ifndef STRUCTURE_ADD_POLYMER
#define STRUCTURE_ADD_POLYMER


#include "structures.hpp"


// Add polymer to the system
bool Structure::add_polymer(OptionsParser &o)
{
    float polymer_bond_length(o.get<float>("polymer_bond_length"));
    float lj_bead_radius_soft(o.get<float>("lj_bead_radius_soft"));
    float too_close_threshold_mmt(o.get<float>("too_close_threshold_mmt"));
    float too_close_threshold_soft(o.get<float>("too_close_threshold_soft"));
    size_t polymer_atom_type(o.get<size_t>("polymer_atom_type"));
    size_t polymer_bond_type(o.get<size_t>("polymer_bond_type"));
    size_t polymerization(o.get<size_t>("polymerization"));

    #ifdef DETAILED_OUTPUT  // Print parameters of polymer addition
        std::cout << "**********\n"
                  << "Structure.hpp add_polymer input parameters:\n"
                  << "\npolymer_bond_length = " << polymer_bond_length
                  << "\nlj_bead_radius = " << lj_bead_radius_soft
                  << "\ntoo_close_threshold_mmt = " << too_close_threshold_mmt
                  << "\ntoo_close_threshold_soft = " << too_close_threshold_soft
                  << "\npolymer_atom_type = " << polymer_atom_type
                  << "\npolymer_bond_type = " << polymer_bond_type
                  << "\npolymerization = " << polymerization << std::endl;
    #endif

    srand(time(NULL));
    float lx = this->xhi - this->xlo;
    float ly = this->yhi - this->ylo;
    float lz = this->zhi - this->zlo;
    size_t fails_done = 0;
    // TODO adjust fails allowed
    size_t fails_allowed = 500;
    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    float close_r_sq_mmt = pow(too_close_threshold_mmt, 2)
        * pow(lj_bead_radius_soft, 2);
    float close_r_sq_soft = pow(too_close_threshold_soft, 2)
        * pow(lj_bead_radius_soft, 2);
    while (fails_done < fails_allowed && new_atoms.size() != polymerization)
      {
        float x;
        float y;
        float z;
        if (new_atoms.size() == 0)
          {
            float x_coeff = (float)(rand()) / (float)(RAND_MAX);// - 0.5;
            float y_coeff = (float)(rand()) / (float)(RAND_MAX);// - 0.5;
            float z_coeff = (float)(rand()) / (float)(RAND_MAX);// - 0.5;
            x = this->xlo + lx * x_coeff;
            y = this->ylo + ly * y_coeff;
            z = this->zlo + lz * z_coeff;
          }
        else
          {
            x = new_atoms[new_atoms.size() - 1].x;
            y = new_atoms[new_atoms.size() - 1].y;
            z = new_atoms[new_atoms.size() - 1].z;
          }
        float theta_coeff = (float)(rand()) / (float)(RAND_MAX);
        float phi_coeff = (float)(rand()) / (float)(RAND_MAX);
        float theta = theta_coeff * M_PI;
        float phi = phi_coeff * 2*M_PI;
        float r = polymer_bond_length * lj_bead_radius_soft;
        x += r * sin(theta) * cos(phi);
        y += r * sin(theta) * sin(phi);
        z += r * cos(theta);
        bool is_close = false;
        for (auto & atom : this->_atoms)
          {
            float dx = std::min(lx - std::fabs(atom.second.x - x),
                                std::fabs(atom.second.x - x));
            float dy = std::min(ly - std::fabs(atom.second.y - y),
                                std::fabs(atom.second.y - y));
            float dz = std::min(lz - std::fabs(atom.second.z - z),
                                std::fabs(atom.second.z - z));
            float dr = dx*dx + dy*dy + dz*dz;
            if ((dr < close_r_sq_mmt && atom.second.phase == "filler")
                || dr < close_r_sq_soft)
              {
                is_close = true;
                break;
              }
          }
        if (is_close)
          {
            fails_done++;
            continue;
          }
        for (auto & atom : new_atoms)
          {
            float dx = std::min(lx - std::fabs(atom.x - x), std::fabs(atom.x - x));
            float dy = std::min(ly - std::fabs(atom.y - y), std::fabs(atom.y - y));
            float dz = std::min(lz - std::fabs(atom.z - z), std::fabs(atom.z - z));
            float dr = dx*dx + dy*dy + dz*dz;
            if (dr < close_r_sq_soft)
              {
                is_close = true;
                break;
              }
          }
        if (is_close)
          {
            fails_done++;
            continue;
          }
        new_atoms.push_back(Atom(x, y, z, 0, polymer_atom_type, 0, 0, 0,
            "polymer"));
      }
    if (new_atoms.size() != polymerization)
      {
        return false;
      }
    for (size_t idx = 0; idx < new_atoms.size(); ++idx)
      {
        this->_atoms[this->_atoms_count + 1 + idx] = new_atoms[idx];
        if (idx < new_atoms.size() - 1)
          {
            this->_bonds[this->_bonds_count + 1 + idx] = Bond(
                polymer_bond_type,
                this->_atoms_count + 1 + idx,
                this->_atoms_count + 1 + idx + 1); 
          }
      }
    this->_atoms_count += polymerization;
    this->_bonds_count += polymerization - 1;

    return true;
}


#endif  // STRUCTURE_ADD_POLYMER