#define _USE_MATH_DEFINES  // define some mathematical constants (M_PI, ...)


#include <cmath>
#include <iostream>


#include "src/options_parser.hpp"
#include "src/write_data.hpp"


int main()
{
    bool verbose = false;
    #ifdef DETAILED_OUTPUT
        verbose = true;
    #endif

    // Parse options
    OptionsParser o("options_periodic");

    float mmt_real_thickness(o.get<float>("mmt_real_thickness"));
    float real_interlayer(o.get<float>("real_interlayer"));
    float dpd_rho(o.get<float>("dpd_rho"));
    float lj_bead_radius_soft(o.get<float>("lj_bead_radius_soft"));
    float lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay"));
    float real_r_c(o.get<float>("real_r_c"));
    size_t platelet_edge(o.get<size_t>("platelet_edge"));
    size_t tail_length(o.get<size_t>("tail_length"));
    size_t polymerization(o.get<size_t>("polymerization"));
    size_t modifiers_count_preset(o.get<size_t>("modifiers_count_preset"));

    #ifdef DETAILED_OUTPUT  // Print all read options
      {
        std::cout << "**********\n";
        std::cout << "Parsed options:\n";
        o.print_options(true);
        std::cout << "**********\n";
      }
    #endif

    // Check options:
    if (mmt_real_thickness <= 0
        || real_interlayer <= 0
        || platelet_edge <= 0)
      {
        std::cout << "Parsed options problem: mmt geometry\n"
          << "\tthickness: " << mmt_real_thickness << " [Angstoms]\n"
          << "\tinterlayer: " << real_interlayer << " [Angstoms]\n"
          << "\tsize: " << platelet_edge << " [beads]\n";
        return 0;
      }

    if (dpd_rho <= 0
        || lj_bead_radius_soft <= 0
        || lj_bead_radius_clay <= 0
        || real_r_c <= 0)
      {
        std::cout << "Parsed options problem: DPD parameters\n"
            << "\tnumeric DPD density: " << dpd_rho
            << "\n\tbead radius soft: "
                << lj_bead_radius_clay << " [Angstroms]"
            << "\n\tbead radius clay: " << lj_bead_radius_soft
                << " [Angstroms]\n"
            << "\n\treal_r_c: " << real_r_c << " [Angstroms]\n";
        return 0;
      }
    if (tail_length <= 0)
      {
        std::cout << "Parsed options problem: modifier\n"
          << "\ttail lenght: " << tail_length << " [beads]\n";
        return 0;
      }
    if (polymerization <= 0)
      {
        std::cout << "Parsed options problem: polymer\n"
            << "\tpolymerization: " << polymerization << " [beads]\n";
        return 0;
      }


    // Compute modifiers count per one lamella:
    float mmt_real_edge = real_r_c * 2*lj_bead_radius_clay * platelet_edge;
    float real_mmt_area = pow(mmt_real_edge, 2);
    float exchange_surface_density = 0.015;  // see magic_constants.md
    size_t charged_count = modifiers_count_preset;   // TODO this is not ok

    #ifdef DETAILED_OUTPUT  // Print information about modifier
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
    size_t mmt_atoms_count = platelet_edge * platelet_edge * 2;
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
    float xy_size = platelet_edge * 2*lj_bead_radius_clay;
    float z_size = 4 * lj_bead_radius_clay + real_interlayer / real_r_c;

    #ifdef DETAILED_OUTPUT  // Print information about box
      {
        std::cout << "**********\n";
        std::cout << "Box computations:\n";
        std::cout << "xy_size (in lj units) = " << xy_size << std::endl;
        std::cout << "z_size (in lj units) = " << z_size << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // Adjust polymers count:
    float lj_cell_volume = pow(xy_size, 2) * z_size;
    size_t all_beads_count = dpd_rho * lj_cell_volume;
    size_t mmt_beads_count = 2 * pow(platelet_edge, 2);
    size_t single_mod_beads_count = 1 + tail_length;
    size_t all_mods_beads_count = single_mod_beads_count * charged_count;
    size_t single_polymer_beads_count = polymerization;
    size_t polymers_count = round((all_beads_count
        -mmt_beads_count - all_mods_beads_count) / polymerization);

    #ifdef DETAILED_OUTPUT  // Print information about polymers
      {
        std::cout << "**********\n";
        std::cout << "Polymers computations:\n";
        std::cout << "polymers_count = " << polymers_count << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // *************************************************************************
    //                             Start main algorithm
    // *************************************************************************
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
    AddMmtPeriodicParameters parameters(o, 0, charged_count);
    bool status = s.add_mmt_periodic(parameters);
    size_t mmt_atoms = s.atoms().size();

    #ifdef DETAILED_OUTPUT  // Print information about mmt addition status
      {
        std::cout << "**********\n";
        std::cout << "MMT addition:\n";
        std::cout << "atoms_count = " << s.atoms().size() << std::endl;
        std::cout << "**********\n";
      }
    #endif

    #ifdef PARTIAL_DATAFILES
      {
        std::string data_out_fname("incomplete_periodic_mmt");
        data_out_fname += "_r" + std::to_string(platelet_radius);
        data_out_fname += "_n" + std::to_string(stacking);
        data_out_fname += ".data";
        std::cout << "Writing incomplete data into " << data_out_fname
            << std::endl;
        write_data(data_out_fname, s);
      }
    #endif

    std::cout << "7\n";

    // Add modifiers
    size_t modifiers_done = 0;
    size_t modifiers_fails_done = 0;
    size_t modifiers_fails_allowed = charged_count;
    while (modifiers_done < charged_count
           && modifiers_fails_done < modifiers_fails_allowed)
      {
        #ifdef DETAILED_OUTPUT  // Print information about last step
          {
            std::cout << "**********\n";
            std::cout << "Modifier addition: "
                      <<  modifiers_done << " of " << charged_count
                      << "; fais: " << modifiers_fails_done
                      << " of " << modifiers_fails_allowed << "\n";
            std::cout << "**********\n";
          }
        #endif

        AddModifierGalleryParameters parameters(o, xy_size/2, -xy_size/2);
        bool status = s.add_modifier_gallery(parameters);
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

    std::cout << "8\n";

    #ifdef DETAILED_OUTPUT  // Print overall information about modifier addition
      {
        std::cout << "**********\n";
        std::cout << "Modifier addition after all: "
                  <<  modifiers_done << " of " << charged_count
                  << "; fais: " << modifiers_fails_done
                  << " of " << modifiers_fails_allowed << "\n";
        std::cout << "**********\n";
      }
    #endif

    #ifdef PARTIAL_DATAFILES  // Output incomplete data into file
      {
        std::string data_out_fname("incomplete_peridic_mmt");
        data_out_fname += "_r" + std::to_string(platelet_radius);
        data_out_fname += "_n" + std::to_string(stacking);
        data_out_fname += "_mod_n" + std::to_string(modifiers_done);
        data_out_fname += "_tail" + std::to_string(tail_length);
        data_out_fname += ".data";
        std::cout << "Writing incomplete data into " << data_out_fname
            << std::endl;
        write_data(data_out_fname, s);
      }
    #endif

    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    size_t polymers_fails_allowed = std::min(size_t(10000),
        polymers_count * polymerization);
    while (polymers_done < polymers_count
        && polymers_fails_done < polymers_fails_allowed)
      {
        #ifdef DETAILED_OUTPUT  // Print information about last step (sometimes)
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

    #ifdef DETAILED_OUTPUT  // Print overall information about polymer addition
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
    data_out_fname += "_edge" + std::to_string(platelet_edge);
    data_out_fname += "_mod_n" + std::to_string(modifiers_done);
    data_out_fname += "_tail" + std::to_string(tail_length);
    data_out_fname += "_poly_p" + std::to_string(polymerization);
    data_out_fname += "_n" + std::to_string(polymers_done);
    data_out_fname += ".data";
    write_data(data_out_fname, s);

    #ifdef GENERAL_OUTPUT  // Print summary structural information
      {
        std::cout << "Structure created:"
            << "\n\tAtoms: " << s.atoms().size() << std::endl
            << "\n\tBonds: " << s.bonds().size() << std::endl
            << "\n\tBox side (x,y): " << xy_size << std::endl
            << "\n\tBox side (z): " << z_size << std::endl
            << std::endl;
        std::cout << "Atoms:"
            << "\n\tMMT:      1 - "  << mmt_atoms
            << "\n\tmodifier: " << mmt_atoms + 1
                << " - " << mmt_atoms + modifier_atoms
            << "\n\tpolymer:  " << mmt_atoms + modifier_atoms + 1
                << " - " << mmt_atoms + modifier_atoms + polymer_atoms
                << std::endl;

        std::cout << "Resulting structure written into:\n\t"
            << data_out_fname << std::endl;

        float DPD_rho = float(s.atoms().size())
            / (s.xhi-s.xlo) / (s.yhi-s.ylo) / (s.zhi-s.zlo);
        float CEC(93 * 90*84/real_mmt_area * modifiers_done/108);
        std::cout << "Physical parameters:\n"
            << "\tDPD_rho (numerical): " << DPD_rho
            << "\n\tCEC: " << CEC << std::endl;
      }
    #endif

    return 0;
}
