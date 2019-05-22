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
      {
        std::cout << "**********\n";
        std::cout << "Parsed options:\n";
        o.print_options(true);
        std::cout << "**********\n";
      }
    #endif

    // Compute general parameters
    // 1.35
    //float real_bead_radius = std::cbrt(o.lx * o.ly * (o.lz - o.mmt_real_thickness)
    //    / o.md_soft_atoms / (4/3 * M_PI));
    float real_bead_radius = 1.35;
    float real_cutoff = std::cbrt(o.dpd_rho * 4/3 * M_PI) * real_bead_radius;
    float lj_interlayer = o.real_interlayer * o.lj_bead_radius / real_bead_radius;
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "General parameters:\n";
        std::cout << "real_bead_radius = " << real_bead_radius << std::endl;
        std::cout << "real_cutoff = " << real_cutoff << std::endl;
        std::cout << "lj_interlayer = " << lj_interlayer << std::endl;
        std::cout << "**********\n";
      }
    #endif


    // Compute modifiers count per one lamella:
    float real_mmt_area = M_PI * pow(o.platelet_radius * real_bead_radius, 2);
    //size_t charged_count = 0.015 * real_mmt_area * o.stacking;
    size_t charged_count = 20;
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Derived parameters:\n";
        std::cout << "real_mmt_area = " << real_mmt_area << std::endl;
        std::cout << "charged_count (per stack) = " << charged_count << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // Compute box size
    float min_height = o.real_interlayer * o.stacking  // top and bottom
        + 4 * real_bead_radius * o.stacking;
    float min_width = 2*o.platelet_radius * 2*real_bead_radius
        * std::max(float(1.25), o.planar_expansion_coeff);
    float cube_edge = std::max(min_width, min_height)
        * o.lj_bead_radius / real_bead_radius;
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Box computations:\n";
        std::cout << "min_height = " << min_height << std::endl;
        std::cout << "min_width = " << min_width << std::endl;
        std::cout << "cube_edge (in lj units) = " << cube_edge << std::endl;
        std::cout << "mmt real radius: " << 2*o.lj_bead_radius *o.platelet_radius
            << std::endl;
        std::cout << "radii: " << o.lj_bead_radius << " " << real_bead_radius
            << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // Adjust polymers count:
    float free_volume = pow(cube_edge, 3)
        -4 * 4*pow(o.platelet_radius, 2) * pow(o.lj_bead_radius, 3)
        -charged_count * (1 + o.tail_length) * 4*pow(o.lj_bead_radius, 3);
    float polymer_volume = o.polymerization * pow(real_bead_radius, 3);
    size_t polymers_count = o.dpd_rho * free_volume / polymer_volume;
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Polymers computations:\n";
        std::cout << "Free volume = " << free_volume << std::endl;
        std::cout << "Polymer volume = " << polymer_volume << std::endl;
        std::cout << "polymers_count = " << polymers_count << std::endl;
        std::cout << "**********\n";
      }
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
        bool status = s.add_mmt_circular(o, 0, 0, 0, charged_count);
        dz = lj_interlayer + 4 * o.lj_bead_radius;
      }
    for (size_t idx = 0; idx < half_stacking; ++idx)
      {
        bool status = s.add_mmt_circular(o, 0, 0, dz, charged_count);
        status = s.add_mmt_circular(o, 0, 0, -dz, charged_count);
        dz += lj_interlayer + 4 * o.lj_bead_radius;
      }
    size_t mmt_atoms = s.atoms().size();
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "MMT addition:\n";
        std::cout << "atoms_count = " << s.atoms().size() << std::endl;
        std::cout << "**********\n";

        write_data("isolated_mmt.data", s);
      }
    #endif

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
    while (modifiers_done < charged_count * o.stacking
        && modifiers_fails_done < modifiers_fails_allowed)
      {
        #ifdef DEBUG
          {
            std::cout << "**********\n";
            std::cout << "Modifier addition: "
                      <<  modifiers_done << " of " << charged_count * o.stacking
                      << "; fais: " << modifiers_fails_done
                      << " of " << modifiers_fails_allowed << "\n";
            std::cout << "**********\n";
          }
        #endif
        size_t idx = rand() % galleries.size();
        float bottom = galleries[idx].first;
        float top = galleries[idx].second;
        bool status = s.add_modifier_gallery(o, top, bottom);
        if (status)
          {
            modifiers_done++;
          }
        else
          {
            modifiers_fails_done ++;
          }
      }
    size_t modifier_atoms = s.atoms().size() - mmt_atoms;
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Modifier addition after all: "
                  <<  modifiers_done << " of " << charged_count
                  << "; fais: " << modifiers_fails_done
                  << " of " << modifiers_fails_allowed << "\n";
        std::cout << "**********\n";

        write_data("isolated_mmt_mod.data", s);
      }
    #endif

    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    size_t polymers_fails_allowed = polymers_count * o.polymerization;
    while (polymers_done < polymers_count
        && polymers_fails_done < polymers_fails_allowed)
      {
        #ifdef DEBUG
          {
            std::cout << "**********\n";
            std::cout << "Polymer addition: "
                      <<  polymers_done << " of " << polymers_count
                      << "; fais: " << polymers_fails_done
                      << " of " << polymers_fails_allowed << "\n";
            std::cout << "**********\n";
          }
        #endif
        bool status = s.add_polymer(o);
        if (status)
          {
            polymers_done++;
          }
        else
          {
            polymers_fails_done ++;
          }
      }
    size_t polymer_atoms = s.atoms().size() - mmt_atoms - modifier_atoms;
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Polymers addition after all: "
                  <<  polymers_done << " of " << polymers_count
                  << "; fais: " << polymers_fails_done
                  << " of " << polymers_fails_allowed << "\n";
        std::cout << "**********\n";
      }
    #endif

    write_data("isolated_mmt_mod_poly.data", s);

    #ifdef DEBUG
      {
        std::cout << "Structure created:\n"
            << "\tAtoms: " << s.atoms().size() << std::endl
            << "\tBonds: " << s.bonds().size() << std::endl
            << "\tBox side: " << cube_edge << std::endl
           << "\tDPD_rho: " << s.atoms().size() / pow(cube_edge, 3) << std::endl;

        std::cout << "MMT: 1 - " << mmt_atoms << std::endl
            << "modifier: " << mmt_atoms + 1 << " - " << mmt_atoms + modifier_atoms
            << std::endl
            << "polymer: " << mmt_atoms + modifier_atoms + 1 << " - "
                << mmt_atoms + modifier_atoms + polymer_atoms << std::endl;
      }
    #endif

    return 0;
}
