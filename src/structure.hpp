#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP


#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>


class Options;


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

    bool add_mmt_circular(size_t platelet_radius, float bead_radius,
        size_t atom_type, size_t bond_type_edge, size_t bond_type_diagonal,
        float x=0, float y=0, float z=0,
        size_t charged_count=0, float bead_charge=0)
      {
        #ifdef DEBUG
        std::cout << "**********\n";
        std::cout << "Structure.hpp add_mmt_circular input parameters:\n";
        std::cout << "platelet_radius = " << platelet_radius << std::endl;
        std::cout << "bead_radius = " << bead_radius << std::endl;
        std::cout << "atom_type = " << atom_type << std::endl;
        std::cout << "bond_type_edge = " << bond_type_edge << std::endl;
        std::cout << "bond_type_diagonal = " << bond_type_diagonal << std::endl;
        std::cout << "x = " << x << std::endl;
        std::cout << "y = " << y << std::endl;
        std::cout << "z = " << z << std::endl;
        std::cout << "charged_count = " << charged_count << std::endl;
        std::cout << "bead_charge = " << bead_charge << std::endl;
        std::cout << "**********\n";
        #endif
        std::vector<Atom> new_atoms;
        std::vector<Bond> new_bonds;
        size_t new_idx = 1;
        // Add atoms and bonds top-bottom
        std::map<int, std::map<int, std::map<std::string, size_t> > > atom_ids;
        for (int idxx = -(int)platelet_radius + 1; idxx < (int)platelet_radius;
            ++idxx)
          {
            atom_ids[idxx] = std::map<int, std::map<std::string, size_t> >();
            float dx = idxx * 2 * bead_radius;
            for (int idxy = -(int)platelet_radius + 1;
                idxy < (int)platelet_radius; ++idxy)
              {
                if (pow(platelet_radius, 2) < pow(idxx, 2) + pow(idxy, 2))
                  {
                    continue;
                  }
                atom_ids[idxx][idxy] = std::map<std::string, size_t>();
                float dy = idxy * 2 * bead_radius;
                new_atoms.push_back(Atom(x + dx, y + dy, z - bead_radius,
                    0, atom_type, 0, 0, 0, "filler"));
                atom_ids[idxx][idxy]["bottom"] = new_idx;
                new_atoms.push_back(Atom(x + dx, y + dy, z + bead_radius,
                    0, atom_type, 0, 0, 0, "filler"));
                atom_ids[idxx][idxy]["top"] = new_idx + 1;
                new_bonds.push_back(Bond(bond_type_edge, new_idx, new_idx + 1));
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
                    new_bonds.push_back(Bond(bond_type_edge,
                        atom_ids[itx->first - 1][ity->first]["top"],
                        atom_ids[itx->first][ity->first]["top"]));
                    new_bonds.push_back(Bond(bond_type_edge,
                        atom_ids[itx->first - 1][ity->first]["bottom"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                if (atom_ids[itx->first].find(ity->first - 1) !=
                           atom_ids[itx->first].end())
                  {
                    new_bonds.push_back(Bond(bond_type_edge,
                        atom_ids[itx->first][ity->first - 1]["top"],
                        atom_ids[itx->first][ity->first]["top"]));
                    new_bonds.push_back(Bond(bond_type_edge,
                        atom_ids[itx->first][ity->first - 1]["bottom"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                // Diagonal bonds
                if (atom_ids.find(itx->first - 1) != atom_ids.end()
                    && atom_ids[itx->first - 1].find(ity->first - 1) !=
                           atom_ids[itx->first - 1].end())
                  {
                    new_bonds.push_back(Bond(bond_type_diagonal,
                        atom_ids[itx->first - 1][ity->first - 1]["top"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                if (atom_ids.find(itx->first - 1) != atom_ids.end()
                    && atom_ids[itx->first - 1].find(ity->first + 1) !=
                           atom_ids[itx->first - 1].end())
                  {
                    new_bonds.push_back(Bond(bond_type_diagonal,
                        atom_ids[itx->first - 1][ity->first + 1]["top"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }
                if (atom_ids.find(itx->first + 1) != atom_ids.end()
                    && atom_ids[itx->first + 1].find(ity->first - 1) !=
                           atom_ids[itx->first + 1].end())
                  {
                    new_bonds.push_back(Bond(bond_type_diagonal,
                        atom_ids[itx->first + 1][ity->first - 1]["top"],
                        atom_ids[itx->first][ity->first]["bottom"]));
                  }

                if (atom_ids.find(itx->first + 1) != atom_ids.end()
                    && atom_ids[itx->first + 1].find(ity->first + 1) !=
                           atom_ids[itx->first + 1].end())
                  {
                    new_bonds.push_back(Bond(bond_type_diagonal,
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
                new_atoms[*it].q = bead_charge;
              }
          }
        #ifdef DEBUG
        std::cout << "**********\n";
        std::cout << "MMT addition in structure.hpp:\n";
        std::cout << "atoms_count = " << new_atoms.size() << std::endl;
        std::cout << "bonds_count = " << new_bonds.size() << std::endl;
        std::cout << "**********\n";
        #endif
        // Append created structure to the existing structure
        std::cout << "Old atoms count\n" << this->_atoms_count << "\n";
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

    bool add_modifier_gallery(float bead_charge, size_t tail_length,
        size_t head_type, size_t tail_type, size_t head_tail_type,
        size_t tail_tail_type, float planar_expansion_coeff,
        float lj_bead_radius, float top, float bottom)
      {
        srand(time(NULL));
        float lx = this->xhi - this->xlo;
        float ly = this->yhi - this->ylo;
        float lz = this->zhi - this->zlo;
        float interlayer = top - bottom;
        size_t fails_done = 0;
        size_t fails_allowed = 100;
        std::vector<Atom> new_atoms;
        std::vector<Bond> new_bonds;
        while (fails_done < fails_allowed
            && new_atoms.size() != 1 + tail_length)
          {
            float x;
            float y;
            float z;
            if (new_atoms.size() == 0)
              {
                float x_coeff = (float)(rand()) / (float)(RAND_MAX) - 0.5;
                float y_coeff = (float)(rand()) / (float)(RAND_MAX) - 0.5;
                float z_coeff = (float)(rand()) / (float)(RAND_MAX) - 0.5;
                x = this->xlo + lx/2 + lx/2 / planar_expansion_coeff * x_coeff;
                y = this->ylo + ly/2 + ly/2 / planar_expansion_coeff * y_coeff;
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
            // TODO bond length
            x += 2.25 * lj_bead_radius * sin(theta) * cos(phi);
            y += 2.25 * lj_bead_radius * sin(theta) * sin(phi);
            z += 2.25 * lj_bead_radius * cos(theta);
            bool is_close = false;
            for (auto & atom : this->_atoms)
              {
                float dx = atom.second.x - x;
                float dy = atom.second.y - y;
                float dz = atom.second.z - z;
                float dr = dx*dx + dy*dy + dz*dz;
                // TODO threshold
                if (dr < 3.5 * pow(lj_bead_radius, 2))
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
                float dx = atom.x - x;
                float dy = atom.y - y;
                float dz = atom.z - z;
                float dr = dx*dx + dy*dy + dz*dz;
                // TODO threshold
                if (dr < 3.5 * pow(lj_bead_radius, 2))
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
            size_t type = tail_type;
            if (new_atoms.size() == 0)
              {
                q = bead_charge;
                type = head_type;
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
                this->_bonds[this->_bonds_count + 1] = Bond(head_tail_type,
                    this->_atoms_count + 1 + idx,
                    this->_atoms_count + 1 + idx + 1); 
              }
            else if (idx < new_atoms.size() - 1)
              {
                this->_bonds[this->_bonds_count + 1 + idx] = Bond(tail_tail_type,
                    this->_atoms_count + 1 + idx,
                    this->_atoms_count + 1 + idx + 1); 
              }
          }
        this->_atoms_count += tail_length + 1;
        this->_bonds_count += tail_length;
        return true;
      }

private:
    size_t _atoms_count;
    size_t _bonds_count;
    std::map<int, Atom> _atoms;
    std::map<int, Bond> _bonds;
};


#endif  // STRUCTURE_HPP
