#include "structures.hpp"


// Add modifier into the space sconstrained by z = top and z = bottom
bool Structure::add_modifier_gallery(AddModifierGalleryParameters &parameters)
{
    // Unpacking passed paramters:
    float top = parameters.top;
    float bottom = parameters.bottom;
    float ht_len(parameters.modifier_head_tail_bond_length);
    float tt_len(parameters.modifier_tail_tail_bond_length);
    float lj_br_soft(parameters.lj_bead_radius_soft);
    float lj_br_clay(parameters.lj_bead_radius_clay);
    float bead_charge = parameters.bead_charge;
    float platelet_closing = parameters.platelet_closing;
    size_t tail_length = parameters.tail_length;
    size_t modifier_head_atom_type = parameters.modifier_head_atom_type;
    size_t modifier_tail_atom_type = parameters.modifier_tail_atom_type;
    size_t head_tail_type = parameters.head_tail_type;
    size_t tail_tail_type = parameters.tail_tail_type;
    size_t platelet_radius = parameters.platelet_radius;
    std::string mode = parameters.mode;
    // Calculate minimal distances based on parameters
    float r2_min_mmt(pow(parameters.too_close_threshold_mmt, 2)
                     * pow(lj_br_soft, 2));
    float r2_min_soft(pow(parameters.too_close_threshold_soft, 2)
                      * pow(lj_br_soft, 2));

    // Start algorithm
    srand(time(NULL));
    float lx = this->xhi - this->xlo;
    float ly = this->yhi - this->ylo;
    float lz = this->zhi - this->zlo;
    float r_platelet_lj((2 * platelet_radius * lj_br_clay * platelet_closing));
    float interlayer = top - bottom;
    size_t fails_done = 0;

    // TODO adjust fails allowed
    size_t fails_allowed = 100;
    std::vector<Atom> new_atoms;
    std::vector<Bond> new_bonds;

    // Start molecule construction
    while (fails_done < fails_allowed && new_atoms.size() != 1 + tail_length)
      {
        float x;
        float y;
        float z;
        if (new_atoms.size() == 0)
          {
            if (mode == "isolated")
              {
                float r_coeff(float(rand()) / float(RAND_MAX));
                float tmp_r(r_platelet_lj * r_coeff);
                float alpha = 2*M_PI * (float)(rand()) / (float)(RAND_MAX);
                float z_coeff = (float)(rand()) / (float)(RAND_MAX);
                x = this->xlo + lx/2 + tmp_r * cos(alpha);
                y = this->ylo + ly/2 + tmp_r * sin(alpha);
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
            r = ht_len * lj_br_soft;
          }
        else
          {
            r = tt_len * lj_br_soft;
          }
        x += r * sin(theta) * cos(phi);
        y += r * sin(theta) * sin(phi);
        z += r * cos(theta);
        bool is_close = false;
        // Check whether molecule fits into existing structure
        for (auto & atom : this->_atoms)
          {
            float dx(std::fabs(atom.second.x - x));
            float dy(std::fabs(atom.second.y - y));
            float dz(std::fabs(atom.second.z - z));
            dx = std::min(lx - dx, dx);
            dy = std::min(ly - dy, dy);
            dz = std::min(lz - dz, dz);
            float dr = dx*dx + dy*dy + dz*dz;
            if ((dr < r2_min_mmt && atom.second.phase == "filler")
                || (dr < r2_min_soft && atom.second.phase != "filler"))
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
        // Check whether new molecule atom fits existing new molecule atoms
        for (auto & atom : new_atoms)
          {
            float dx(std::fabs(atom.x - x));
            float dy(std::fabs(atom.y - y));
            float dz(std::fabs(atom.z - z));
            dx = std::min(lx - dx, dx);
            dy = std::min(ly - dy, dy);
            dz = std::min(lz - dz, dz);
            float dr = dx*dx + dy*dy + dz*dz;
            if (dr < r2_min_soft)
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
