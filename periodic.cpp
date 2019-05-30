#define _USE_MATH_DEFINES  // define some mathematical constants (M_PI, ...)


#include <cmath>
#include <iostream>


#include "src/options.hpp"
#include "src/structure.hpp"
#include "src/write_data.hpp"


int main()
{
    bool verbose = true;

    // Parse options
    Options o("options_periodic");
    #ifdef DEBUG
      {
        // Output for checking wether parsed options are read correctly
        std::cout << "**********\n";
        std::cout << "Parsed options:\n";
        o.print_options(true);
        std::cout << "**********\n";
      }
    #endif

    // Check options:
    if (o.mmt_real_thickness <= 0
        || o.real_interlayer <= 0
        || o.platelet_edge <= 0)
      {
        std::cout << "Parsed options problem: mmt geometry\n"
          << "\tthickness: " << o.mmt_real_thickness << " [Angstoms]\n"
          << "\tinterlayer: " << o.real_interlayer << " [Angstoms]\n"
          << "\tsize: " << o.platelet_edge << " [beads]\n";
        return 0;
      }
    if (o.dpd_rho <= 0 || o.lj_bead_radius <= 0)
      {
        std::cout << "Parsed options problem: DPD parameters\n"
            << "\tnumeric DPD density: " << o.dpd_rho
            << "\n\tbead radius: " << o.lj_bead_radius << " [Angstroms]\n";
        return 0;
      }
    if (o.tail_length <= 0)
      {
        std::cout << "Parsed options problem: modifier\n"
          << "\ttail lenght: " << o.tail_length << " [beads]\n";
        return 0;
      }
    if (o.polymerization <= 0)
      {
        std::cout << "Parsed options problem: polymer\n"
            << "\tpolymerization: " << o.polymerization << " [beads]\n";
        return 0;
      }

    // Compute general parameters

    // In MD modeling:
    // 48420 atoms = 9 * 5380 =
    //   = 9 * (720 clay atoms + 840 mod atoms + 3820 poly atoms)
    //   = 9 * (720 clay atoms + 12*70 mod atoms + 10*382 poly atoms)
    //                                                 /\
    //                                                 20 monomers
    // This transforma into:
    // N soft beads = 9 * (12 * 3 mod beads + 10 * 20 poly beads) = 
    //   = 9 * (36 mod beads + 200 poly beads) = 2124 soft beads
    // -> real_bead_volume = 80 * 94 * (64 - mmt_real_thickness) / 2124

    float real_bead_volume = 80 * 94 * (64 - o.mmt_real_thickness) / 2124;
    float real_bead_radius = std::cbrt(real_bead_volume *3/4 / M_PI);
    float real_cutoff = std::cbrt(o.dpd_rho * 4/3 * M_PI) * real_bead_radius;
    float lj_interlayer = o.real_interlayer * o.lj_bead_radius / real_bead_radius;
    #ifdef DEBUG
      {
        // Output of general parameters used in the further computations
        std::cout << "**********\n";
        std::cout << "General parameters:\n";
        std::cout << "real_bead_radius = " << real_bead_radius << std::endl;
        std::cout << "real_cutoff = " << real_cutoff << std::endl;
        std::cout << "lj_interlayer = " << lj_interlayer << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // Compute modifiers count per one lamella:
    float real_mmt_area = pow(o.platelet_edge * 2*real_bead_radius, 2);
    float exchange_surface_density = 0.015;  // see magic_constants.md
    //size_t charged_count = exchange_surface_density * real_mmt_area;
    size_t charged_count = o.modifiers_count_preset;   // TODO this is not ok
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Derived parameters:\n";
        std::cout << "real_mmt_area = " << real_mmt_area << std::endl;
        std::cout << "charged_count = " << charged_count
            << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // MMT atoms count should be multiple of that of the modifiers
    size_t mmt_atoms_count = o.platelet_edge*o.platelet_edge * 2;
    float max_exchanged_fraction = 0.25; // TODO: get real value
        // (though 0.25 seems quite ok).
    if (mmt_atoms_count * max_exchanged_fraction < charged_count)
      {
        std::cout << "Too high CEC in periodic:\n"
            << "\tMMT atoms: " << mmt_atoms_count
            << "\n\tcharged_count: " << charged_count
            << "\n\tratio: " << float(charged_count) / float(mmt_atoms_count)
            << "\n\t while max allowed ratio is: " << max_exchanged_fraction
            << std::endl;
        return 0;
      }

    // Compute box size
    float xy_size = o.platelet_edge * 2*o.lj_bead_radius;
    float z_size = 2*o.lj_bead_radius + 2* o.real_interlayer / o.lj_bead_radius;
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Box computations:\n";
        std::cout << "xy_size (in lj units) = " << xy_size << std::endl;
        std::cout << "z_size (in lj units) = " << z_size << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // Adjust polymers count:
    float bv = 4/3*M_PI * pow(o.lj_bead_radius, 3);  // bead volume
    float lj_cell_volume = pow(xy_size, 2) * z_size;
    size_t all_beads_count = o.dpd_rho * lj_cell_volume;
    size_t mmt_beads_count = 2 * pow(o.platelet_edge, 2);
    size_t single_mod_beads_count = 1 + o.tail_length;
    size_t all_mods_beads_count = single_mod_beads_count * charged_count;
    size_t single_polymer_beads_count = o.polymerization;
    size_t polymers_count = round((all_beads_count
        -mmt_beads_count - all_mods_beads_count) / o.polymerization);

    /*
    std::cout << "all_beads_count: " << all_beads_count
        << "\nmmt_beads_count: " << mmt_beads_count
        << "\nsingle_mod_beads_count: " << single_mod_beads_count
        << "\nall_mods_beads_count: " << all_mods_beads_count
        << "\nsingle_polymer_beads_count: " << single_polymer_beads_count
        << "\npolymers_count: " << polymers_count << std::endl;

    size_t num_density = lj_cell_volume / 
        (polymers_count * o.polymerization
         +all_mods_beads_count + mmt_beads_count);
    std::cout << "num_density: " << num_density << std::endl;
    std::cout << lj_cell_volume / all_beads_count << std::endl;

    return 0;*/


    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "Polymers computations:\n";
        std::cout << "polymers_count = " << polymers_count << std::endl;
        std::cout << "**********\n";
      }
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
    // TODO: adjustable bead charge
    bool status = s.add_mmt_periodic(o, 0, charged_count);
    size_t mmt_atoms = s.atoms().size();
    #ifdef DEBUG
      {
        std::cout << "**********\n";
        std::cout << "MMT addition:\n";
        std::cout << "atoms_count = " << s.atoms().size() << std::endl;
        std::cout << "**********\n";
      }
    #endif

    #ifdef PARTIAL_DATAFILES
      {
        //write_data("periodic_mmt.data", s);
      }
    #endif

    // Add modifiers
    // structure.hpp, add_modifiers:
    size_t modifiers_done = 0;
    size_t modifiers_fails_done = 0;
    size_t modifiers_fails_allowed = charged_count;
    while (modifiers_done < charged_count
           && modifiers_fails_done < modifiers_fails_allowed)
      {
        #ifdef DEBUG
          {
            std::cout << "**********\n";
            std::cout << "Modifier addition: "
                      <<  modifiers_done << " of " << charged_count
                      << "; fais: " << modifiers_fails_done
                      << " of " << modifiers_fails_allowed << "\n";
            std::cout << "**********\n";
          }
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
      {
        std::cout << "**********\n";
        std::cout << "Modifier addition after all: "
                  <<  modifiers_done << " of " << charged_count
                  << "; fais: " << modifiers_fails_done
                  << " of " << modifiers_fails_allowed << "\n";
        std::cout << "**********\n";
      }
    #endif

    #ifdef PARTIAL_DATAFILES
      {
        //write_data("periodic_mmt_mod.data", s);
      }
    #endif

    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    // TODO
    size_t polymers_fails_allowed = std::min(size_t(10000),
        polymers_count * o.polymerization);
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

    std::string data_out_fname("periodic_mmt");
    data_out_fname += "_edge" + std::to_string(o.platelet_edge);
    data_out_fname += "_mod_n" + std::to_string(modifiers_done);
    data_out_fname += "_tail" + std::to_string(o.tail_length);
    data_out_fname += "_poly_p" + std::to_string(o.polymerization);
    data_out_fname += "_n" + std::to_string(polymers_done);
    data_out_fname += ".data";

    //write_data("isolated_mmt_mod_poly.data", s);
    write_data(data_out_fname, s);

    #ifdef DEBUG
      {
        std::cout << "Structure created:\n"
            << "\tAtoms: " << s.atoms().size() << std::endl
            << "\tBonds: " << s.bonds().size() << std::endl
            << "\tBox side (x,y): " << xy_size << std::endl
            << "\tBox side (z): " << z_size << std::endl
            << "\tDPD_rho: " << s.atoms().size() / pow(xy_size, 2) / z_size
            << std::endl;

        std::cout << "Indices:\n"
            << "MMT: 1 - "  << mmt_atoms << std::endl
            << "modifier: " << mmt_atoms + 1
                << " - " << mmt_atoms + modifier_atoms
                << std::endl
            << "polymer: " << mmt_atoms + modifier_atoms + 1
                << " - " << mmt_atoms + modifier_atoms + polymer_atoms
                << std::endl;

        std::cout << "real_bead_radius " << real_bead_radius << std::endl;

        // In composites studied my md:
        // lx = 90, ly = 8, cell_mass = 720 (40 atoms), n_cells = 162
        double mmt_real_surface = pow(xy_size * real_bead_radius, 2);
        double mmt_real_mass = 720*162 * mmt_real_surface / 90/80;
        double mmt_real_bead_mass = mmt_real_mass / mmt_atoms;
        std::cout << "mmt real bead mass " << mmt_real_bead_mass << std::endl;
      }
    #endif  // DEBUG

    // Check numerical density:
    std::cout << "Resulting numerical density: "
        << s.atoms().size() / (s.xhi-s.xlo) / (s.yhi-s.ylo) / (s.zhi-s.zlo)
        << std::endl;

    return 0;
}
