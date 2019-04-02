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
    float real_mmt_area = 4 * pow(o.platelet_edge * real_bead_radius, 2);
    // CEC = 92.6
    // cell_square = 400
    // cell_mass = 720
    // nx*ny = 18
    //     coeff = CEC * cell_mass
    // 108 
    // 0.015 = 108 / ly*lx
    // 1/720 * 2/3 = 93.6 meq/100g
    size_t charged_count = 0.015 * real_mmt_area;
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "Derived parameters:\n";
    std::cout << "real_mmt_area = " << real_mmt_area << std::endl;
    std::cout << "charged_count (per stack == at all) = " << charged_count
        << std::endl;
    std::cout << "**********\n";
    #endif

    // Compute box size
    float xy_size = o.platelet_edge * 2*o.lj_bead_radius;
    float z_size = 2*o.lj_bead_radius + o.real_interlayer / o.lj_bead_radius;
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "Box computations:\n";
    std::cout << "xy_size (in lj units) = " << xy_size << std::endl;
    std::cout << "z_size (in lj units) = " << z_size << std::endl;
    std::cout << "**********\n";
    #endif

    // Adjust polymers count:
    float free_volume = pow(xy_size, 3) * z_size
        -4 * pow(o.platelet_edge, 2) * pow(o.lj_bead_radius, 3)
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
    s.xlo = -xy_size/2;
    s.xhi = xy_size/2;
    s.ylo = -xy_size/2;
    s.yhi = xy_size/2;
    s.zlo = -z_size/2;
    s.zhi = z_size/2;

    // Add single MMT platelet
    bool status = s.add_mmt_periodic(o, 0, charged_count);
    size_t mmt_atoms = s.atoms().size();
    #ifdef DEBUG
    std::cout << "**********\n";
    std::cout << "MMT addition:\n";
    std::cout << "atoms_count = " << s.atoms().size() << std::endl;
    std::cout << "**********\n";
    #endif

    write_data("periodic_mmt.data", s);

    // Add modifiers
    size_t modifiers_done = 0;
    size_t modifiers_fails_done = 0;
    size_t modifiers_fails_allowed = charged_count;
    while (modifiers_done < charged_count * o.stacking
        && modifiers_fails_done < modifiers_fails_allowed)
      {
        #ifdef DEBUG
        std::cout << "**********\n";
        std::cout << "Modifier addition: "
                  <<  modifiers_done << " of " << charged_count * o.stacking
                  << "; fais: " << modifiers_fails_done
                  << " of " << modifiers_fails_allowed << "\n";
        std::cout << "**********\n";
        #endif
        bool status = s.add_modifier_gallery(o, xy_size/2, -xy_size/2);
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
    std::cout << "**********\n";
    std::cout << "Modifier addition after all: "
              <<  modifiers_done << " of " << charged_count
              << "; fais: " << modifiers_fails_done
              << " of " << modifiers_fails_allowed << "\n";
    std::cout << "**********\n";
    #endif

    write_data("periodic_mmt_mod.data", s);

    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    // TODO
    size_t polymers_fails_allowed = std::min(size_t(10000),
        polymers_count * o.polymerization);
    polymers_count = std::min(size_t(10000), polymers_count);
    while (polymers_done < polymers_count
        && polymers_fails_done < polymers_fails_allowed)
      {
        #ifdef DEBUG
        std::cout << "**********\n";
        std::cout << "Polymer addition: "
                  <<  polymers_done << " of " << polymers_count
                  << "; fais: " << polymers_fails_done
                  << " of " << polymers_fails_allowed << "\n";
        std::cout << "**********\n";
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
    std::cout << "**********\n";
    std::cout << "Polymers addition after all: "
              <<  polymers_done << " of " << polymers_count
              << "; fais: " << polymers_fails_done
              << " of " << polymers_fails_allowed << "\n";
    std::cout << "**********\n";
    #endif

    write_data("periodic_mmt_mod_poly.data", s);

    std::cout << "Structure created:\n"
        << "\tAtoms: " << s.atoms().size() << std::endl
        << "\tBonds: " << s.bonds().size() << std::endl
        << "\tBox side (x,y): " << xy_size << std::endl
        << "\tBox side (z): " << z_size << std::endl
        << "\tDPD_rho: " << s.atoms().size() / pow(xy_size, 2) * z_size
        << std::endl;

    std::cout << "MMT: 1 - " << mmt_atoms << std::endl
        << "modifier: " << mmt_atoms + 1 << " - " << mmt_atoms + modifier_atoms
        << std::endl
        << "polymer: " << mmt_atoms + modifier_atoms + 1 << " - "
            << mmt_atoms + modifier_atoms + polymer_atoms << std::endl;

    std::cout << s.atoms()[1].phase << " " << (s.atoms()[1].phase == "filler")
        << std::endl;

    return 0;
}
