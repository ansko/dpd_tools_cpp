#include "structures.hpp"


// Add modifier into the space sconstrained by z = top and z = bottom
bool Structure::add_modifier_gallery(AddModifierGalleryParameters &parameters)
{
    // Unpacking
    float top = parameters.top;
    float bottom = parameters.bottom;
    float modifier_head_tail_bond_length
        = parameters.modifier_head_tail_bond_length;
    float modifier_tail_tail_bond_length
        = parameters.modifier_tail_tail_bond_length;
    float lj_bead_radius_soft = parameters.lj_bead_radius_soft;
    float lj_bead_radius_clay = parameters.lj_bead_radius_clay;
    float too_close_threshold_mmt = parameters.too_close_threshold_mmt;
    float too_close_threshold_soft = parameters.too_close_threshold_soft;
    float bead_charge = parameters.bead_charge;
    float platelet_closing = parameters.platelet_closing;
    size_t tail_length = parameters.tail_length;
    size_t modifier_head_atom_type = parameters.modifier_head_atom_type;
    size_t modifier_tail_atom_type = parameters.modifier_tail_atom_type;
    size_t head_tail_type = parameters.head_tail_type;
    size_t tail_tail_type = parameters.tail_tail_type;
    size_t platelet_radius = parameters.platelet_radius;
    std::string mode = parameters.mode;

    srand(time(NULL));
    float lx = this->xhi - this->xlo;
    float ly = this->yhi - this->ylo;
    float lz = this->zhi - this->zlo;
    float r_platelet = platelet_radius*2;
    float r_platelet_lj = (2 * platelet_radius * platelet_closing)
        * lj_bead_radius_clay;
    float interlayer = top - bottom;
    size_t fails_done = 0;

    // TODO adjust fails allowed
    size_t fails_allowed = 100;
    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    float close_r_sq_mmt = pow(too_close_threshold_mmt, 2)
        * pow(lj_bead_radius_soft, 2);
    float close_r_sq_soft = pow(too_close_threshold_soft, 2)
        * pow(lj_bead_radius_soft, 2);
    while (fails_done < fails_allowed && new_atoms.size() != 1 + tail_length)
      {
        float x;
        float y;
        float z;
        if (new_atoms.size() == 0)
          {
            if (mode == "isolated")
              {
                float alpha = 2*M_PI * (float)(rand()) / (float)(RAND_MAX);
                float z_coeff = (float)(rand()) / (float)(RAND_MAX);
                x = this->xlo + lx/2 + r_platelet_lj * cos(alpha);
                y = this->ylo + ly/2 + r_platelet_lj * sin(alpha);
                z = bottom + interlayer * z_coeff;
              }
            else
              {
                float x_coeff = (float)(rand()) / (float)(RAND_MAX);
                float y_coeff = (float)(rand()) / (float)(RAND_MAX);
                float z_coeff = (float)(rand()) / (float)(RAND_MAX);
                x = this->xlo + lx * x_coeff;
                y = this->ylo + ly * y_coeff;
                z = bottom + interlayer * z_coeff;
              }
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
        float r;
        if (new_atoms.size() == 0)
          {
            r = modifier_head_tail_bond_length * lj_bead_radius_soft;
          }
        else
          {
            r = modifier_tail_tail_bond_length * lj_bead_radius_soft;
          }
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
                || (dr < close_r_sq_soft && atom.second.phase != "filler"))
              {
                is_close = true;
                break;
              }
          }
        if (is_close)
          {
            fails_done ++;
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
            fails_done ++;
            continue;
          }
        float q = 0;
        size_t type = modifier_tail_atom_type;
        if (new_atoms.size() == 0)
          {
            q = bead_charge;
            type = modifier_head_atom_type;
          }
        new_atoms.push_back(Atom(x, y, z, q, type, 0, 0, 0, "modifier"));
      }
    if (new_atoms.size() != 1 + tail_length)
      {
        return false;
      }
    for (size_t idx = 0; idx < new_atoms.size(); ++idx)
      {
        this->_atoms[this->_atoms_count + 1 + idx] = new_atoms[idx];
        if (idx == 0)
          {
            this->_bonds[this->_bonds_count + 1] = Bond(
                head_tail_type,
                this->_atoms_count + 1 + idx,
                this->_atoms_count + 1 + idx + 1); 
          }
        else if (idx < new_atoms.size() - 1)
          {
            this->_bonds[this->_bonds_count + 1 + idx] = Bond(
                tail_tail_type,
                this->_atoms_count + 1 + idx,
                this->_atoms_count + 1 + idx + 1); 
          }
      }
    this->_atoms_count += tail_length + 1;
    this->_bonds_count += tail_length;

    return true;
}
