#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP


#include <cmath>

#include <map>
#include <mutex>
#include <set>
#include <vector>

#include "options_parser.hpp"


struct AddMmtCircularParameters
{
public:
    AddMmtCircularParameters(OptionsParser &o, float x, float y, float z,
        size_t charged_count)
    : x(x), y(y), z(z), charged_count(charged_count),
      lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay")),
      bead_charge(o.get<float>("bead_charge")),
      platelet_closing(o.get<float>("platelet_closing")),
      platelet_radius(o.get<size_t>("platelet_radius")),
      mmt_atom_type(o.get<size_t>("mmt_atom_type")),
      mmt_edge_bond_type(o.get<size_t>("mmt_edge_bond_type")),
      mmt_diagonal_bond_type(o.get<size_t>("mmt_diagonal_bond_type"))
    {
    }

    float x;               // x of the circluar platelet center
    float y;               // y of the circluar platelet center
    float z;               // z of the circluar platelet center
    size_t charged_count;  // number of platelet atoms carrying non zero charge

    float lj_bead_radius_clay;
    float bead_charge;
    float platelet_closing;
    size_t platelet_radius;
    size_t mmt_atom_type;
    size_t mmt_edge_bond_type;
    size_t mmt_diagonal_bond_type;
};


struct AddMmtPeriodicParameters
{
public:
    AddMmtPeriodicParameters(OptionsParser &o, float z=0, size_t
        charged_count=0)
    : z(z), charged_count(charged_count),
      lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay")),
      bead_charge(o.get<float>("bead_charge")),
      platelet_edge(o.get<size_t>("platelet_edge")),
      mmt_atom_type(o.get<size_t>("mmt_atom_type")),
      mmt_edge_bond_type(o.get<size_t>("mmt_edge_bond_type")),
      mmt_diagonal_bond_type(o.get<size_t>("mmt_diagonal_bond_type"))
    {}

    float z;               // z of the circluar platelet center
    size_t charged_count;  // number of platelet atoms carrying non zero charge
    float lj_bead_radius_clay;
    float bead_charge;
    size_t platelet_edge;
    size_t mmt_atom_type;
    size_t mmt_edge_bond_type;
    size_t mmt_diagonal_bond_type;
};


struct AddModifierGalleryParameters
{
public:
    AddModifierGalleryParameters (OptionsParser &o, float top, float bottom)
    : top(top), bottom(bottom),
      modifier_head_tail_bond_length(
          o.get<float>("modifier_head_tail_bond_length")),
      modifier_tail_tail_bond_length(
          o.get<float>("modifier_tail_tail_bond_length")),
      lj_bead_radius_soft(o.get<float>("lj_bead_radius_soft")),
      lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay")),
      too_close_threshold_mmt(o.get<float>("too_close_threshold_mmt")),
      too_close_threshold_soft(o.get<float>("too_close_threshold_soft")),
      bead_charge(o.get<float>("bead_charge")),
      platelet_closing(o.get<float>("platelet_closing")),
      tail_length(o.get<size_t>("tail_length")),
      modifier_head_atom_type(o.get<size_t>("modifier_head_atom_type")),
      modifier_tail_atom_type(o.get<size_t>("modifier_tail_atom_type")),
      head_tail_type(o.get<size_t>("head_tail_type")),
      tail_tail_type(o.get<size_t>("tail_tail_type")),
      platelet_radius(o.get<size_t>("platelet_radius"))
    {}

    float top;     // top z-coordintate of the gallery
    float bottom;  // bottom z-coordinate of the gallery
    float modifier_head_tail_bond_length;
    float modifier_tail_tail_bond_length;
    float lj_bead_radius_soft;
    float lj_bead_radius_clay;
    float too_close_threshold_mmt;
    float too_close_threshold_soft;
    float bead_charge;
    float platelet_closing;
    size_t tail_length;
    size_t modifier_head_atom_type;
    size_t modifier_tail_atom_type;
    size_t head_tail_type;
    size_t tail_tail_type;
    size_t platelet_radius;
};


struct AddPolymerParameters
{
public:
    AddPolymerParameters(OptionsParser &o)
    : polymer_bond_length(o.get<float>("polymer_bond_length")),
      lj_bead_radius_soft(o.get<float>("lj_bead_radius_soft")),
      too_close_threshold_mmt(o.get<float>("too_close_threshold_mmt")),
      too_close_threshold_soft(o.get<float>("too_close_threshold_soft")),
      polymer_atom_type(o.get<size_t>("polymer_atom_type")),
      polymer_bond_type(o.get<size_t>("polymer_bond_type")),
      polymerization(o.get<size_t>("polymerization"))
    {}

    float polymer_bond_length;
    float lj_bead_radius_soft;
    float too_close_threshold_mmt;
    float too_close_threshold_soft;
    size_t polymer_atom_type;
    size_t polymer_bond_type;
    size_t polymerization;
};


struct AddPolymerParallelParameters
{
public:
    AddPolymerParallelParameters(OptionsParser &o, size_t idx_x, size_t nx,
        size_t idx_y, size_t ny, size_t idx_z, size_t nz)
    : idx_x(idx_x), nx(nx), idx_y(idx_y), ny(ny), idx_z(idx_z), nz(nz),
      polymer_bond_length(o.get<float>("polymer_bond_length")),
      lj_bead_radius_soft(o.get<float>("lj_bead_radius_soft")),
      too_close_threshold_mmt(o.get<float>("too_close_threshold_mmt")),
      too_close_threshold_soft(o.get<float>("too_close_threshold_soft")),
      polymer_atom_type(o.get<size_t>("polymer_atom_type")),
      polymer_bond_type(o.get<size_t>("polymer_bond_type")),
      polymerization(o.get<size_t>("polymerization"))
    {}

    size_t idx_x;  // index along x of thread running this function
    size_t nx;     // full number of thread along x
    size_t idx_y;  // ... same for other axes ...
    size_t ny;     // ... same for other axes ...
    size_t idx_z;  // ... same for other axes ...
    size_t nz;     // ... same for other axes ...
    float polymer_bond_length;
    float lj_bead_radius_soft;
    float too_close_threshold_mmt;
    float too_close_threshold_soft;
    size_t polymer_atom_type;
    size_t polymer_bond_type;
    size_t polymerization;
};


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

    Atom(float x, float y, float z, float q, size_t type, int nx=0, int ny=0,
         int nz=0, std::string phase="undefined")
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


// A big class to store vectors of atoms and bonds in the box
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

    std::map<int, Atom> atoms() { return this->_atoms; }
    std::map<int, Bond> bonds() { return this->_bonds; }

    bool add_mmt_circular(AddMmtCircularParameters &parameters);
    bool add_mmt_periodic(AddMmtPeriodicParameters &parameters);
    bool add_modifier_gallery(AddModifierGalleryParameters &parameters);
    bool add_polymer(AddPolymerParameters &paramters);
    bool add_polymer_parallel(AddPolymerParallelParameters &parameters);

private:
    size_t _atoms_count;
    size_t _bonds_count;
    std::map<int, Atom> _atoms;
    std::map<int, Bond> _bonds;

    std::mutex _mtx_polymers_addition;  // for parallel algorithms
};


#endif  // STRUCTURES_HPP
