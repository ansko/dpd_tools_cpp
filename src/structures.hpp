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
    : o(o), x(x), y(y), z(z), charged_count(charged_count)
    {}

    OptionsParser o;
    float x;               // x of the circluar platelet center
    float y;               // y of the circluar platelet center
    float z;               // z of the circluar platelet center
    size_t charged_count;  // number of platelet atoms carrying non zero charge
};


struct AddMmtPeriodicParameters
{
public:
    AddMmtPeriodicParameters(OptionsParser &o, float z=0, size_t
        charged_count=0)
    : o(o), z(z), charged_count(charged_count)
    {}

    OptionsParser o;
    float z;               // z of the circluar platelet center
    size_t charged_count;  // number of platelet atoms carrying non zero charge
};


struct AddModifierGalleryParameters
{
public:
    AddModifierGalleryParameters (OptionsParser &o, float top, float bottom)
    : o(o), top(top), bottom(bottom)
    {}

    OptionsParser o;
    float top;     // top z-coordintate of the gallery
    float bottom;  // bottom z-coordinate of the gallery
};


struct AddPolymerParallelParameters
{
public:
    AddPolymerParallelParameters(OptionsParser &o, size_t idx_x, size_t nx,
        size_t idx_y, size_t ny, size_t idx_z, size_t nz)
    : o(o), idx_x(idx_x), nx(nx), idx_y(idx_y), ny(ny), idx_z(idx_z), nz(nz)
    {}

    OptionsParser o;
    size_t idx_x;  // index along x of thread running this function
    size_t nx;     // full number of thread along x
    size_t idx_y;  // ... same for other axes ...
    size_t ny;
    size_t idx_z;
    size_t nz;
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
    bool add_polymer(OptionsParser &o);
    bool add_polymer_parallel(AddPolymerParallelParameters &parameters);

private:
    size_t _atoms_count;
    size_t _bonds_count;
    std::map<int, Atom> _atoms;
    std::map<int, Bond> _bonds;

    std::mutex _mtx_polymers_addition;  // for parallel algorithms
};


#endif  // STRUCTURES_HPP
