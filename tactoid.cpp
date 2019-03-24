#define _USE_MATH_DEFINES


#include <cmath>
#include <iostream>


#include "src/options.hpp"
#include "src/structure.hpp"
#include "src/write_data.hpp"


int main()
{
    bool verbose = true;

    // Parse options
    Options o("options");
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "Parsed options:\n";
    o.print_options(true);
    std::cout << "**********\n";
    #endif

    // Compute general parameters
    float real_bead_radius = std::cbrt(o.lx * o.ly * (o.lz - o.mmt_real_thickness)
        / o.md_soft_atoms / (4/3 * M_PI));
    float real_cutoff = std::cbrt(o.dpd_rho * 4/3 * M_PI) * real_bead_radius;
    float lj_interlayer = o.real_interlayer * o.lj_bead_radius / real_bead_radius;
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "General parameters:\n";
    std::cout << "real_bead_radius = " << real_bead_radius << std::endl;
    std::cout << "real_cutoff = " << real_cutoff << std::endl;
    std::cout << "lj_interlayer = " << lj_interlayer << std::endl;
    std::cout << "**********\n";
    #endif

    // Compute modifiers count per one lamella:
    float real_mmt_area = M_PI * pow(o.platelet_radius* real_bead_radius, 2);
    size_t charged_count = 0.015 * real_mmt_area * o.stacking;
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "Derived parameters:\n";
    std::cout << "real_mmt_area = " << real_mmt_area << std::endl;
    std::cout << "charged_count (per stack) = " << charged_count << std::endl;
    std::cout << "**********\n";
    #endif

    // Compute box size
    float min_height = o.real_interlayer * (o.stacking - 1)
        + 4 * real_bead_radius * o.stacking;
    float min_width = 2 * o.platelet_radius * real_bead_radius
        * o.planar_expansion_coeff;
    float cube_edge = std::max(cube_edge, min_height)
        * o.lj_bead_radius / real_bead_radius;
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "Box computations:\n";
    std::cout << "min_height = " << min_height << std::endl;
    std::cout << "min_width = " << min_width << std::endl;
    std::cout << "cube_edge (in lj units) = " << cube_edge << std::endl;
    std::cout << "**********\n";
    #endif

    // Adjust polymers count:
    float free_volume = pow(cube_edge, 3)
        -4 * 4*pow(o.platelet_radius, 2) * pow(o.lj_bead_radius, 3)
        -charged_count * (1 + o.tail_length) * 4*pow(o.lj_bead_radius, 3);
    float polymer_volume = o.polymerization * pow(real_bead_radius, 3);
    size_t polymers_count = o.dpd_rho * free_volume / polymer_volume;
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "Polymers computations:\n";
    std::cout << "polymers_count = " << polymers_count << std::endl;
    std::cout << "**********\n";
    #endif

    // Start main algorithm
    Structure s;

    // Set calcualted box
    s.xlo = -cube_edge;
    s.xhi = cube_edge;
    s.ylo = -cube_edge;
    s.yhi = cube_edge;
    s.zlo = -cube_edge;
    s.zhi = cube_edge;

    // Add MMT stack
    size_t half_stacking = o.stacking / 2;
    float dz = lj_interlayer / 2 + 2 * o.lj_bead_radius;
    if (o.stacking % 2 != 0)
      {
        bool status = s.add_mmt_circular(o.platelet_radius, o.lj_bead_radius,
            o.mmt_atom_type, o.mmt_edge_bond_type, o.mmt_diagonal_bond_type,
            0, 0, 0, 0, 0);
        dz = lj_interlayer + 4 * o.lj_bead_radius;
      }
    for (size_t idx = 0; idx < half_stacking; ++idx)
      {
        bool status = s.add_mmt_circular(o.platelet_radius, o.lj_bead_radius,
            o.mmt_atom_type, o.mmt_edge_bond_type, o.mmt_diagonal_bond_type,
            0, 0, dz, 0, 0);
        status = s.add_mmt_circular(o.platelet_radius, o.lj_bead_radius,
            o.mmt_atom_type, o.mmt_edge_bond_type, o.mmt_diagonal_bond_type,
            0, 0, -dz, 0, 0);
        dz += lj_interlayer + 4 * o.lj_bead_radius;
      }

    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "MMT addition:\n";
    std::cout << "atoms_count = " << s.atoms().size() << std::endl;
    std::cout << "**********\n";
    #endif

    write_data("mmt_test.data", s);

    // Add modifiers
    std::vector<std::pair<float, float> > galleries;
    if (o.stacking % 2 == 0)
      {
        for (size_t idx = 0; idx < half_stacking; ++idx)
          {
            float top = lj_interlayer / 2
                + idx * (lj_interlayer + 4 * o.lj_bead_radius);
            galleries.push_back(std::pair<float, float>(top - lj_interlayer, top));
          }
      }
    else
      {
        for (size_t idx = 0; idx < half_stacking; ++idx)
          {
            float top = lj_interlayer + 2 * o.lj_bead_radius
                + idx * (lj_interlayer + 4 * o.lj_bead_radius);
            galleries.push_back(std::pair<float, float>(top - lj_interlayer, top));
          }
      }
    size_t modifiers_done = 0;
    size_t modifiers_fails_done = 0;
    size_t modifiers_fails_allowed = charged_count;
    while (modifiers_done < charged_count
        && modifiers_fails_done < modifiers_fails_allowed)
      {
        #ifdef DEBUG
        std::cout << "**********\n";
        std::cout << "Modifier addition: "
                  <<  modifiers_done << " of " << charged_count
                  << "; fais: " << modifiers_fails_done
                  << " of " << modifiers_fails_allowed << "\n";
        std::cout << "**********\n";
        #endif
        size_t idx = rand() % galleries.size();
        float bottom = galleries[idx].first;
        float top = galleries[idx].second;
        bool status = s.add_modifier_gallery(o.bead_charge, o.tail_length,
            o.modifier_head_atom_type, o.modifier_tail_atom_type,
            o.head_tail_type, o.tail_tail_type, o.planar_expansion_coeff,
            o.lj_bead_radius, top, bottom);
        if (status)
          {
            modifiers_done++;
          }
        else
          {
            modifiers_fails_done ++;
          }
      }

    write_data("mmt_mod_test.data", s);

    return 0;
}
/*
def make_tactoid():
    top = lj_interlayer/2
    bottom = -lj_interlayer/2
    #print('atoms in mmt', mmt_atoms_count)
    modifiers_done = 0
    fails_done = 0
    fails_allowed = charged_count * 10
    while modifiers_done < charged_count and fails_done < fails_allowed:
        status, new_structure = between_clays(
            top=top, bottom=bottom,
            structure=structure,
            head_charge=bead_charge,
            tail_length=tail_length,
            bead_radius=lj_bead_radius,
            head_atom_type=2,
            tail_atom_type=3,
            head_tail_bond_type=3,
            tail_tail_bond_type=4)
        if status:
            modifiers_done += 1
            structure = new_structure
        else:
            fails_done += 1
    print('Modifiers done after all:', modifiers_done)
    print('charged', len(structure['atoms']))

    polymers_done = 0
    polymers_fails_done = 0
    polymers_fails_allowed = 10000 #10 * polymers_count
    while (polymers_done < polymers_count and
        polymers_fails_done < polymers_fails_allowed):
        old_structure = copy.deepcopy(structure)
        if polymerization > 1:
            status, new_structure = random_chain(
                polymerization=polymerization,
                bead_radius=lj_bead_radius,
                atom_type=4,
                bond_type=5,
                start_atom_id=len(structure['atoms']) + 1,
                start_bond_id=len(structure['bonds']) + 1,
                structure=structure)
        else:
            status, new_structure = polymer_chain(
                polymerization=polymerization,
                bead_radius=bead_radius,
                atom_type=4,
                start_atom_id=len(structure['atoms']) + 1,
                structure=structure)
        if status:
            structure = new_structure
            polymers_done += 1
            #if polymers_done % polymers_count // 10 == 0:
        else:
            structure = copy.deepcopy(old_structure)
            polymers_fails_done += 1
        print('polymers: {0}, fails: {1}/{2}'.format(
            polymers_done, polymers_fails_done, polymers_fails_allowed))
    print('Polymers done after all:', polymers_done)

    # Summary
    print('Structure:\n\t{0} atoms\n\t{1} bonds'.format(
        len(structure['atoms']), len(structure['bonds'])))
    print('Cell: volume={0}, dpd_density={1}'.format(
        cube_edge**3, len(structure['atoms'])/cube_edge**3))
    write_data('tactoid.data', structure)
*/
