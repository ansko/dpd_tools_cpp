#define _USE_MATH_DEFINES  // define some mathematical constants (M_PI, ...)


#include <cmath>
#include <iostream>

#include "src/options_parser.hpp"
#include "src/write_data.hpp"


int main()
{
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

    // Adjust polymers count:
    float lj_cell_volume = pow(xy_size, 2) * z_size;
    size_t all_beads_count = dpd_rho * lj_cell_volume;
    size_t mmt_beads_count = 2 * pow(platelet_edge, 2);
    size_t single_mod_beads_count = 1 + tail_length;
    size_t all_mods_beads_count = single_mod_beads_count * charged_count;
    size_t single_polymer_beads_count = polymerization;
    size_t polymers_count = round((all_beads_count
        -mmt_beads_count - all_mods_beads_count) / polymerization);

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

    // Add modifiers
    size_t modifiers_done = 0;
    size_t modifiers_fails_done = 0;
    size_t modifiers_fails_allowed = charged_count;
    while (modifiers_done < charged_count
           && modifiers_fails_done < modifiers_fails_allowed)
      {
        AddModifierGalleryParameters parameters(o, xy_size/2, -xy_size/2,
                                                "periodic");
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

    if (modifiers_done != charged_count)
      {
        std::cout << "Failed to add all modifiers!\n";
        return 0;
      }

    size_t modifier_atoms = s.atoms().size() - mmt_atoms;
    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    size_t polymers_fails_allowed = std::min(size_t(10000),
        polymers_count * polymerization);
    while (polymers_done < polymers_count
        && polymers_fails_done < polymers_fails_allowed)
      {
        AddPolymerParameters parameters(o);
        bool status = s.add_polymer(parameters);
        if (status)
          {
            polymers_done++;
          }
        else
          {
            polymers_fails_done++;
          }
      }

    if (polymers_done != polymers_count)
      {
        std::cout << "Failed to add all polymers!\n";
        return 0;
      }

    size_t polymer_atoms = s.atoms().size() - mmt_atoms - modifier_atoms;

    std::string data_out_fname("periodic_mmt");
    data_out_fname += "_edge" + std::to_string(platelet_edge);
    data_out_fname += "_mod_n" + std::to_string(modifiers_done);
    data_out_fname += "_tail" + std::to_string(tail_length);
    data_out_fname += "_poly_p" + std::to_string(polymerization);
    data_out_fname += "_n" + std::to_string(polymers_done);
    data_out_fname += ".data";
    write_data(data_out_fname, s);

    std::cout << "Success!\n" << data_out_fname << std::endl;

    return 0;
}
