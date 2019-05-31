#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP


#include <fstream>
#include <unistd.h>

#include <map>
#include <mutex>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "options_common.hpp"
#include "options_periodic.hpp"
#include "options_isolated.hpp"
#include "options_isolated_parallel.hpp"


// A simple container class for atom properties
struct Atom
{
public:
    float x;
    float y;
    float z;
    float q;
    size_t type;
    size_t nx;
    size_t ny;
    size_t nz;
    std::string phase;

    Atom() {};

    Atom(float x, float y, float z, float q, size_t type,
        int nx=0, int ny=0, int nz=0, std::string phase="undefined")
    : x(x), y(y), z(z), q(q), type(type), nx(nx), ny(ny), nz(nz), phase(phase)
    {};
};


// A simple container class for bond properties
struct Bond
{
public:
    size_t type;
    size_t atom_one;
    size_t atom_two;

    Bond() {};

    Bond(size_t type, size_t atom_one, size_t atom_two)
    : type(type), atom_one(atom_one), atom_two(atom_two)
    {};
};


// A very big class for storing and creating composite structure
class Structure
{
public:
    float xlo;
    float xhi;
    float ylo;
    float yhi;
    float zlo;
    float zhi;

    Structure()
    : _atoms_count(0), _bonds_count(0)
    {};

    std::map<int, Atom> atoms()
      {
        return this->_atoms;
      }

    std::map<int, Bond> bonds()
      {
        return this->_bonds;
      }

    // Add an infinite (because of periodic boundary conditions) MMT platelet
    bool add_mmt_circular(OptionsIsolated &o,
        float x=0, float y=0, float z=0, size_t charged_count=0)
      {
        #ifdef DETAILED_OUTPUT  // Print parameters of MMT additions
        {
            std::cout << "**********\n"
                << "\nStructure.hpp add_mmt_circular input parameters:\n"
                << "\nplatelet_radius = " << o.platelet_radius
                << "\nbead_radius_clay = " << o.lj_bead_radius_clay
                << "\natom_type = " << o.mmt_atom_type
                << "\nmmt_edge_bond_type = " << o.mmt_edge_bond_type
                << "\nmmt_diagonal_bond_type = " << o.mmt_diagonal_bond_type
                << "\nx = " << x
                << "\ny = " << y
                << "\nz = " << z
                << "\ncharged_count = " << charged_count
                << "\nbead_charge = " << o.bead_charge
                << "\n**********\n";
          }
        #endif

        std::vector<Atom> new_atoms;
        std::vector<Bond> new_bonds;
        size_t new_idx = 1;

        // Add atoms and bonds top-bottom
        std::map<int, std::map<int, std::map<std::string, size_t> > > atom_ids;
        for (int idxx = -(int)o.platelet_radius + 1; idxx < (int)o.platelet_radius;
            ++idxx)
          {
            atom_ids[idxx] = std::map<int, std::map<std::string, size_t> >();
            float dx = idxx * 2 * o.lj_bead_radius_clay;
            for (int idxy = -(int)o.platelet_radius + 1;
                idxy < (int)o.platelet_radius; ++idxy)
              {
                if (pow(o.platelet_radius, 2) < pow(idxx, 2) + pow(idxy, 2))
                  {
                    continue;
                  }
                atom_ids[idxx][idxy] = std::map<std::string, size_t>();
                float dy = idxy * 2 * o.lj_bead_radius_clay;
                new_atoms.push_back(Atom(x + dx, y + dy, z - o.lj_bead_radius_clay,
                    0, o.mmt_atom_type, 0, 0, 0, "filler"));
                atom_ids[idxx][idxy]["bottom"] = new_idx;
                new_atoms.push_back(Atom(x + dx, y + dy, z + o.lj_bead_radius_clay,
                    0, o.mmt_atom_type, 0, 0, 0, "filler"));
                atom_ids[idxx][idxy]["top"] = new_idx + 1;
                new_bonds.push_back(Bond(o.mmt_edge_bond_type, new_idx,
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
                    new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                        atom_ids[itx->first - 1][ity->first]["top"],
                        atom_ids[itx->first][ity->first]["top"]));
                    new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                        atom_ids[itx->first - 1][ity->first]["bottom"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                if (atom_ids[itx->first].find(ity->first - 1) !=
                           atom_ids[itx->first].end())
                  {
                    new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                        atom_ids[itx->first][ity->first - 1]["top"],
                        atom_ids[itx->first][ity->first]["top"]));
                    new_bonds.push_back(Bond(o.mmt_edge_bond_type,
                        atom_ids[itx->first][ity->first - 1]["bottom"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                // Diagonal bonds
                if (atom_ids.find(itx->first - 1) != atom_ids.end()
                    && atom_ids[itx->first - 1].find(ity->first - 1) !=
                           atom_ids[itx->first - 1].end())
                  {
                    new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                        atom_ids[itx->first - 1][ity->first - 1]["top"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                if (atom_ids.find(itx->first - 1) != atom_ids.end()
                    && atom_ids[itx->first - 1].find(ity->first + 1) !=
                           atom_ids[itx->first - 1].end())
                  {
                    new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                        atom_ids[itx->first - 1][ity->first + 1]["top"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                if (atom_ids.find(itx->first + 1) != atom_ids.end()
                    && atom_ids[itx->first + 1].find(ity->first - 1) !=
                           atom_ids[itx->first + 1].end())
                  {
                    new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
                        atom_ids[itx->first + 1][ity->first - 1]["top"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }

                if (atom_ids.find(itx->first + 1) != atom_ids.end()
                    && atom_ids[itx->first + 1].find(ity->first + 1) !=
                           atom_ids[itx->first + 1].end())
                  {
                    new_bonds.push_back(Bond(o.mmt_diagonal_bond_type,
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
            while (charged_ids.size() < charged_count)
              {
                charged_ids.insert(rand() % new_atoms.size());
              }
            for (auto it = charged_ids.begin(); it != charged_ids.end(); ++it)
              {
                new_atoms[*it].q = -o.bead_charge;
              }
          }
        #ifdef DETAILED_OUTPUT  // Results of MMT addition
          {
            std::cout << "**********\n";
            std::cout << "MMT addition in structure.hpp:\n";
            std::cout << "atoms_count = " << new_atoms.size() << std::endl;
            std::cout << "bonds_count = " << new_bonds.size() << std::endl;
            std::cout << "**********\n";
          }
        #endif

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
      };


    // Add an isolated MMT platelet
    bool add_mmt_periodic(OptionsPeriodic &o, float z=0, size_t charged_count=0)
      {
        #ifdef DETAILED_OUTPUT  // Print parameters of MMT addition
          {
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
          }
        #endif

        std::vector<Atom> new_atoms;
        std::vector<Bond> new_bonds;
        size_t new_idx = 1;

        // Add atoms and bonds top-bottom
        std::map<int, std::map<int, std::map<std::string, size_t> > > atom_ids;
        for (int idxx = 0; idxx < (int)o.platelet_edge; ++idxx)
          {
            atom_ids[idxx] = std::map<int, std::map<std::string, size_t> >();
            float dx = idxx * 2 * o.lj_bead_radius_clay
                       - (float)o.platelet_edge * o.lj_bead_radius_clay
                       + o.lj_bead_radius_clay;
            for (int idxy = 0; idxy < (int)o.platelet_edge; ++idxy)
              {
                atom_ids[idxx][idxy] = std::map<std::string, size_t>();
                float dy = idxy * 2 * o.lj_bead_radius_clay
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
        #ifdef DETAILED_OUTPUT  // Results of MMT addition
          {
            std::cout << "**********\n";
            std::cout << "MMT addition in structure.hpp:\n";
            std::cout << "atoms_count = " << new_atoms.size() << std::endl;
            std::cout << "bonds_count = " << new_bonds.size() << std::endl;
            std::cout << "**********\n";
          }
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
      };


    // Add modifier into the space sconstrained by z = top and z = bottom
    bool add_modifier_gallery(OptionsCommon &o, float top, float bottom)
      {
        #ifdef DETAILED_OUTPUT  // Parameters of modifier addition
          {
            std::cout << "**********\n"
                << "Structure.hpp add_modifier_gallery input parameters:\n"
                << "tail_length = " << o.tail_length
                // << "\nplanar_expansion_coeff = " << o.planar_expansion_coeff
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
          }
        #endif

        srand(time(NULL));
        float lx = this->xhi - this->xlo;
        float ly = this->yhi - this->ylo;
        float lz = this->zhi - this->zlo;
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
        while (fails_done < fails_allowed
               && new_atoms.size() != 1 + o.tail_length)
          {
            float x;
            float y;
            float z;
            if (new_atoms.size() == 0)
              {
                float x_coeff = (float)(rand()) / (float)(RAND_MAX) - 0.5;
                float y_coeff = (float)(rand()) / (float)(RAND_MAX) - 0.5;
                float z_coeff = (float)(rand()) / (float)(RAND_MAX);
                //x = this->xlo + lx/2 + lx/2 * o.planar_expansion_coeff * x_coeff;
                //y = this->ylo + ly/2 + ly/2 * o.planar_expansion_coeff * y_coeff;
                x = this->xlo + lx/2 + lx/4 * x_coeff;  // planar coeff == 2
                y = this->ylo + ly/2 + ly/4 * y_coeff;
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
                float dx = std::min(lx - std::fabs(atom.x - x),
                                    std::fabs(atom.x - x));
                float dy = std::min(ly - std::fabs(atom.y - y),
                                    std::fabs(atom.y - y));
                float dz = std::min(lz - std::fabs(atom.z - z),
                                    std::fabs(atom.z - z));
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


    // Add polymer to the system
    bool add_polymer(OptionsCommon &o)
      {
        #ifdef DETAILED_OUTPUT  // Print parameters of polymer addition
        if (false)
          {
            std::cout << "**********\n"
                << "Structure.hpp add_polymer input parameters:\n"
                << "\npolymer_bond_length = " << o.polymer_bond_length
                << "\nlj_bead_radius = " << o.lj_bead_radius_soft
                << "\ntoo_close_threshold_mmt = " << o.too_close_threshold_mmt
                << "\ntoo_close_threshold_soft = " << o.too_close_threshold_soft
                << "\npolymer_atom_type = " << o.polymer_atom_type
                << "\npolymer_bond_type = " << o.polymer_bond_type
                << "\npolymerization = " << o.polymerization << std::endl;
          }
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
        float close_r_sq_mmt = pow(o.too_close_threshold_mmt, 2)
                               * pow(o.lj_bead_radius_soft, 2);
        float close_r_sq_soft = pow(o.too_close_threshold_soft, 2)
                                * pow(o.lj_bead_radius_soft, 2);
        while (fails_done < fails_allowed
               && new_atoms.size() != o.polymerization)
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
                float dx = std::min(lx - std::fabs(atom.x - x),
                                    std::fabs(atom.x - x));
                float dy = std::min(ly - std::fabs(atom.y - y),
                                    std::fabs(atom.y - y));
                float dz = std::min(lz - std::fabs(atom.z - z),
                                    std::fabs(atom.z - z));
                float dr = dx*dx + dy*dy + dz*dz;
                if (dr < close_r_sq_soft)
                  {
                    is_close = true;
                    //std::cout << "close 2\n";
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
        if (new_atoms.size() != o.polymerization)
          {
            return false;
          }
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
        return true;
      }


    // Parallel version of polymer addition.
    // Crops whole box into nx*ny*nz parallelepipeds and inserts
    // polymer chains into the parallelepiped that is defined by:
    //     0 <= idx_x < nx
    //     0 <= idx_y < ny
    //     0 <= idx_z < nz
    bool add_polymer_parallel(OptionsIsolatedParallel &o,
        size_t idx_x, size_t nx, size_t idx_y, size_t ny, size_t idx_z, size_t nz)
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
        if (false)
          {
            std::cout << "**********\n"
                << "Structure.hpp add_polymer input parameters:\n"
                << "polymer_bond_length = " << o.polymer_bond_length
                << "\nlj_bead_radius_soft = " << o.lj_bead_radius_soft
                << "\ntoo_close_threshold_mmt = " << o.too_close_threshold_mmt
                << "too_close_threshold_soft = " << o.too_close_threshold_soft
                << "polymer_atom_type = " << o.polymer_atom_type
                << "polymer_bond_type = " << o.polymer_bond_type
                << "polymerization = " << o.polymerization << std::endl;
          }
        #endif

        // TODO adjust fails allowed
        size_t fails_allowed = 500;
        std::vector<Atom> new_atoms;
        std::vector<Bond> new_bonds;
        float close_r_sq_mmt = pow(o.too_close_threshold_mmt, 2)
                               * pow(o.lj_bead_radius_soft, 2);
        float close_r_sq_soft = pow(o.too_close_threshold_soft, 2)
                                * pow(o.lj_bead_radius_soft, 2);
        while (fails_done < fails_allowed
               && new_atoms.size() != o.polymerization)
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
                float dx = std::min(lx - std::fabs(atom.x - x),
                                    std::fabs(atom.x - x));
                float dy = std::min(ly - std::fabs(atom.y - y),
                                    std::fabs(atom.y - y));
                float dz = std::min(lz - std::fabs(atom.z - z),
                                    std::fabs(atom.z - z));
                float dr = dx*dx + dy*dy + dz*dz;
                if (dr < close_r_sq_soft)
                  {
                    is_close = true;
                    //std::cout << "close 2\n";
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

private:
    size_t _atoms_count;
    size_t _bonds_count;
    std::map<int, Atom> _atoms;
    std::map<int, Bond> _bonds;

    std::mutex _mtx_polymers_addition;
};


#endif  // STRUCTURE_HPP
