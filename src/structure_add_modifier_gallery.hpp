#ifndef STRUCTURE_ADD_MODIFIER_GALLERY
#define STRUCTURE_ADD_MODIFIER_GALLERY


#include "structures.hpp"


// Add modifier into the space sconstrained by z = top and z = bottom
bool Structure::add_modifier_gallery(OptionsParser &o, float top, float bottom)
{
    #ifdef DETAILED_OUTPUT  // Parameters of modifier addition
        std::cout << "**********\n"
                  << "Structure.hpp add_modifier_gallery input parameters:\n"
                  << "tail_length = " << o.tail_length
                  << "\nmodifier_head_tail_bond_length = "
                      << o.modifier_head_tail_bond_length
                  << "\nmodifier_tail_tail_bond_length = "
                      << o.modifier_tail_tail_bond_length
                  << "\nlj_bead_radius_soft = " << o.lj_bead_radius_soft
                  << "\ntoo_close_threshold_mmt = " << o.too_close_threshold_mmt
                  << "\ntoo_close_threshold_soft = " << o.too_close_threshold_soft
                  << "\nmodifier_head_atom_type = " << o.modifier_head_atom_type
                  << "\nmodifier_tail_atom_type = " << o.modifier_tail_atom_type
                  << "\nbead_charge = " << o.bead_charge
                  << "\nhead_tail_type = " << o.head_tail_type
                  << "\ntail_tail_type = " << o.tail_tail_type
                  << "\ntop = " << top
                  << "\nbottom = " << bottom
                  << "\n**********\n";
    #endif

    srand(time(NULL));
    float lx = this->xhi - this->xlo;
    float ly = this->yhi - this->ylo;
    float lz = this->zhi - this->zlo;
    float r_platelet = o.platelet_radius*2;
    float interlayer = top - bottom;
    size_t fails_done = 0;

    // TODO adjust fails allowed
    size_t fails_allowed = 100;
    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;
    float close_r_sq_mmt = pow(o.too_close_threshold_mmt, 2)
        * pow(o.lj_bead_radius_soft, 2);
    float close_r_sq_soft = pow(o.too_close_threshold_soft, 2)
        * pow(o.lj_bead_radius_soft, 2);
    while (fails_done < fails_allowed && new_atoms.size() != 1 + o.tail_length)
      {
        float x;
        float y;
        float z;
        if (new_atoms.size() == 0)
          {
            float alpha = 2*M_PI * (float)(rand()) / (float)(RAND_MAX);
            //float x_coeff = (float)(rand()) / (float)(RAND_MAX) - 0.5;
            //float y_coeff = (float)(rand()) / (float)(RAND_MAX) - 0.5;
            float z_coeff = (float)(rand()) / (float)(RAND_MAX);
            //x = this->xlo + lx/2 + lx/2 * o.planar_expansion_coeff * x_coeff;
            //y = this->ylo + ly/2 + ly/2 * o.planar_expansion_coeff * y_coeff;
            //x = this->xlo + lx/2 + lx/4 * x_coeff;  // planar coeff == 2
            //y = this->ylo + ly/2 + ly/4 * y_coeff;
            x = this->xlo + lx/2 + r_platelet * cos(alpha);
            y = this->ylo + ly/2 + r_platelet * sin(alpha);
            z = bottom + interlayer * z_coeff;
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
            r = o.modifier_head_tail_bond_length * o.lj_bead_radius_soft;
          }
        else
          {
            r = o.modifier_tail_tail_bond_length * o.lj_bead_radius_soft;
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
        size_t type = o.modifier_tail_atom_type;
        if (new_atoms.size() == 0)
          {
            q = o.bead_charge;
            type = o.modifier_head_atom_type;
          }
        new_atoms.push_back(Atom(x, y, z, q, type, 0, 0, 0, "modifier"));
      }
    if (new_atoms.size() != 1 + o.tail_length)
      {
        return false;
      }
    for (size_t idx = 0; idx < new_atoms.size(); ++idx)
      {
        this->_atoms[this->_atoms_count + 1 + idx] = new_atoms[idx];
        if (idx == 0)
          {
            this->_bonds[this->_bonds_count + 1] = Bond(
                o.head_tail_type,
                this->_atoms_count + 1 + idx,
                this->_atoms_count + 1 + idx + 1); 
          }
        else if (idx < new_atoms.size() - 1)
          {
            this->_bonds[this->_bonds_count + 1 + idx] = Bond(
                o.tail_tail_type,
                this->_atoms_count + 1 + idx,
                this->_atoms_count + 1 + idx + 1); 
          }
      }
    this->_atoms_count += o.tail_length + 1;
    this->_bonds_count += o.tail_length;

    return true;
}


#endif  // STRUCTURE_ADD_MODIFIER_GALLERY
