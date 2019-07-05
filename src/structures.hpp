#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP


#include <mutex>
#include <set>
#include <vector>


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

    bool add_mmt_circular(OptionsParser &o,
        float x=0, float y=0, float z=0, size_t charged_count=0);
    bool add_mmt_periodic(OptionsParser &o, float z=0, size_t charged_count=0);
    bool add_modifier_gallery(OptionsParser &o, float top, float bottom);
    bool add_polymer(OptionsParser &o);
    bool add_polymer_parallel(OptionsParser &o,
        size_t idx_x, size_t nx, size_t idx_y, size_t ny, size_t idx_z, size_t nz);

private:
    size_t _atoms_count;
    size_t _bonds_count;
    std::map<int, Atom> _atoms;
    std::map<int, Bond> _bonds;

    std::mutex _mtx_polymers_addition;  // for parallel algorithms
};


#endif  // STRUCTURES_HPP
