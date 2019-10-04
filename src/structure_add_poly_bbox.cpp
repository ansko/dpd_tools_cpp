#include "structures.hpp"


/* Parallel version of polymer addition:
   add polymer into bbox, which may be equal to the whole cell or
   to its part; when inserting constructing polymer into exsisting
   structure, the mutex is locked.
*/
bool
Structure::add_poly_bbox(BBox &bbox, AddPolymerParallelParameters &parameters)
{
    float polymer_bond_length = parameters.polymer_bond_length;
    float lj_bead_radius_soft = parameters.lj_bead_radius_soft;
    float too_close_threshold_mmt = parameters.too_close_threshold_mmt;
    float too_close_threshold_soft = parameters.too_close_threshold_soft;
    size_t polymer_atom_type = parameters.polymer_atom_type;
    size_t polymer_bond_type = parameters.polymer_bond_type;
    size_t polymerization = parameters.polymerization;

    srand(time(NULL));

    float xlo = bbox.xlo;
    float xhi = bbox.xhi;
    float ylo = bbox.ylo;
    float yhi = bbox.yhi;
    float zlo = bbox.zlo;
    float zhi = bbox.zhi;
    float lx = xhi - xlo;
    float ly = yhi - ylo;
    float lz = zhi - zlo;

    size_t fails_done = 0;

    // TODO adjust fails allowed
    size_t fails_allowed = 500;
    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    float close_r_sq_mmt(pow(too_close_threshold_mmt * lj_bead_radius_soft, 2));
    float close_r_sq_soft(pow(too_close_threshold_soft * lj_bead_radius_soft, 2));

    float r = polymer_bond_length * lj_bead_radius_soft;  // random move step

    while (fails_done < fails_allowed && new_atoms.size() != polymerization)
      {
        float x;
        float y;
        float z;
        if (new_atoms.size() == 0)
          {
            // Randomly creating first atom of the chain
            float x_coeff = (float)(rand()) / (float)(RAND_MAX);
            float y_coeff = (float)(rand()) / (float)(RAND_MAX);
            float z_coeff = (float)(rand()) / (float)(RAND_MAX);
            x = xlo + lx * x_coeff;
            y = ylo + ly * y_coeff;
            z = zlo + lz * z_coeff;
          }
        else
          {
            // Move in a random direction from the last atom appened
            x = new_atoms[new_atoms.size() - 1].x;
            y = new_atoms[new_atoms.size() - 1].y;
            z = new_atoms[new_atoms.size() - 1].z;
            float theta_coeff = (float)(rand()) / (float)(RAND_MAX);
            float phi_coeff = (float)(rand()) / (float)(RAND_MAX);
            float theta = theta_coeff * M_PI;
            float phi = phi_coeff * 2*M_PI;
            x += r * sin(theta) * cos(phi);
            y += r * sin(theta) * sin(phi);
            z += r * cos(theta);
          }

        bool is_close = false;
        for (auto & atom : this->_atoms)
          {
            float dx = std::fabs(atom.second.x - x);
            float dy = std::fabs(atom.second.y - y);
            float dz = std::fabs(atom.second.z - z);
            dx = std::min(lx - dx, dx);
            dy = std::min(ly - dy, dy);
            dz = std::min(lz - dz, dz);
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
            float dx = std::fabs(atom.x - x);
            float dy = std::fabs(atom.y - y);
            float dz = std::fabs(atom.z - z);
            dx = std::min(lx - dx, dx);
            dy = std::min(ly - dy, dy);
            dz = std::min(lz - dz, dz);
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
    if (new_atoms.size() != polymerization)  // incomplete chain constructed
      {
        return false;
      }

    // Insert constructed polymer chain into the exsising structure
    this->_mtx_polymers_addition.lock();

    for (size_t idx = 0; idx < new_atoms.size(); ++idx)
      {
        this->_atoms[this->_atoms_count + 1 + idx] = new_atoms[idx];
        if (idx < new_atoms.size() - 1)
          {
            Bond b(polymer_bond_type, this->_atoms_count + 1 + idx,
                   this->_atoms_count + 1 + idx + 1); 
            this->_bonds[this->_bonds_count + 1 + idx] = b;
          }
      }
    this->_atoms_count += polymerization;
    this->_bonds_count += polymerization - 1;

    this->_mtx_polymers_addition.unlock();

    return true;
}
