#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP


#include <cmath>

#include <map>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include "options_parser.hpp"
#include "parameters_classes.hpp"


// Data structure for storing bounding box
struct BBox
{
    BBox(const float xlo, const float xhi,
         const float ylo, const float yhi,
         const float zlo, const float zhi)
      : xlo(xlo), xhi(xhi), ylo(ylo), yhi(yhi), zlo(zlo), zhi(zhi)
    {};

    const float xlo = 0.111;
    const float xhi = 0.111;
    const float ylo = 0.111;
    const float yhi = 0.111;
    const float zlo = 0.111;
    const float zhi = 0.111;
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

    bool add_mmt_circular_fcc(AddMmtCircularParameters &parameters);

    bool add_mmt_periodic(AddMmtPeriodicParameters &parameters);
    bool add_modifier_gallery(AddModifierGalleryParameters &parameters);
    bool add_polymer_parallel(AddPolymerParallelParameters &parameters);

    bool add_poly_bbox(BBox &bbox, AddPolymerParallelParameters &parameters);

private:
    size_t _atoms_count;
    size_t _bonds_count;
    std::map<int, Atom> _atoms;
    std::map<int, Bond> _bonds;

    std::mutex _mtx_polymers_addition;  // for parallel algorithms
};


#endif  // STRUCTURES_HPP include guard
