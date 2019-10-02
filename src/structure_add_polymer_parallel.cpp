#include "structures.hpp"


// Parallel version of polymer addition.
// Crops whole box into nx*ny*nz parallelepipeds and inserts
// polymer chains into the parallelepiped that is defined by:
//     0 <= idx_x < nx
//     0 <= idx_y < ny
//     0 <= idx_z < nz
bool Structure::add_polymer_parallel(AddPolymerParallelParameters &parameters)
{
    // Unpacking parameters
    size_t idx_x = parameters.idx_x;
    size_t nx = parameters.nx;
    size_t idx_y = parameters.idx_y;
    size_t ny = parameters.ny;
    size_t idx_z = parameters.idx_z;
    size_t nz = parameters.nz;

    float polymer_bond_length = parameters.polymer_bond_length;
    float lj_bead_radius_soft = parameters.lj_bead_radius_soft;
    float too_close_threshold_mmt = parameters.too_close_threshold_mmt;
    float too_close_threshold_soft = parameters.too_close_threshold_soft;
    size_t polymer_atom_type = parameters.polymer_atom_type;
    size_t polymer_bond_type = parameters.polymer_bond_type;
    size_t polymerization = parameters.polymerization;

    srand(time(NULL));
    float lx = (this->xhi - this->xlo) / nx;
    float ly = (this->yhi - this->ylo) / ny;
    float lz = (this->zhi - this->zlo) / nz;
    float my_xlo = this->xlo + lx * idx_x;
    float my_xhi = this->xlo + lx * (idx_x + 1);
    float my_ylo = this->ylo + ly * idx_y;
    float my_yhi = this->ylo + ly * (idx_y + 1);
    float my_zlo = this->zlo + lz * idx_z;
    float my_zhi = this->zlo + lz * (idx_z + 1);
    size_t fails_done = 0;

    // TODO adjust fails allowed
    size_t fails_allowed = 500;
    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    float close_r_sq_mmt = pow(too_close_threshold_mmt, 2)
        * pow(lj_bead_radius_soft, 2);
    float close_r_sq_soft = pow(too_close_threshold_soft, 2)
        * pow(lj_bead_radius_soft, 2);

    float r = polymer_bond_length * lj_bead_radius_soft;

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
            x = my_xlo + lx * x_coeff;
            y = my_ylo + ly * y_coeff;
            z = my_zlo + lz * z_coeff;
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
    if (new_atoms.size() != polymerization)  // incomplete chain
      {
        return false;
      }

    this->_mtx_polymers_addition.lock();

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

    this->_mtx_polymers_addition.unlock();

    return true;
}
