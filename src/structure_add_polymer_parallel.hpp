#ifndef STRUCTURE_ADD_POLYMER_PARALLEL
#define STRUCTURE_ADD_POLYMER_PARALLEL


#include "structures.hpp"


// Parallel version of polymer addition.
// Crops whole box into nx*ny*nz parallelepipeds and inserts
// polymer chains into the parallelepiped that is defined by:
//     0 <= idx_x < nx
//     0 <= idx_y < ny
//     0 <= idx_z < nz
bool Structure::add_polymer_parallel(OptionsParser &o, size_t idx_x, size_t nx,
size_t idx_y, size_t ny, size_t idx_z, size_t nz)
{
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

    #ifdef DETAILED_OUTPUT  // Print parameters of polmyer addition
        std::cout << "**********\n"
                  << "Structure.hpp add_polymer input parameters:\n"
                  << "polymer_bond_length = " << o.polymer_bond_length
                  << "\nlj_bead_radius_soft = " << o.lj_bead_radius_soft
                  << "\ntoo_close_threshold_mmt = " << o.too_close_threshold_mmt
                  << "too_close_threshold_soft = " << o.too_close_threshold_soft
                  << "polymer_atom_type = " << o.polymer_atom_type
                  << "polymer_bond_type = " << o.polymer_bond_type
                  << "polymerization = " << o.polymerization << std::endl;
    #endif

    // TODO adjust fails allowed
    size_t fails_allowed = 500;
    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    float close_r_sq_mmt = pow(o.too_close_threshold_mmt, 2)
        * pow(o.lj_bead_radius_soft, 2);
    float close_r_sq_soft = pow(o.too_close_threshold_soft, 2)
        * pow(o.lj_bead_radius_soft, 2);
    while (fails_done < fails_allowed && new_atoms.size() != o.polymerization)
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
        float r = o.polymer_bond_length * o.lj_bead_radius_soft;
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
        new_atoms.push_back(Atom(x, y, z, 0, o.polymer_atom_type, 0, 0, 0,
            "polymer"));
      }
    if (new_atoms.size() != o.polymerization)  // incomplete chain
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
                o.polymer_bond_type,
                this->_atoms_count + 1 + idx,
                this->_atoms_count + 1 + idx + 1); 
          }
      }
    this->_atoms_count += o.polymerization;
    this->_bonds_count += o.polymerization - 1;

    this->_mtx_polymers_addition.unlock();

    return true;
}


#endif  // STRUCTURE_ADD_POLYMER_PARALLEL
