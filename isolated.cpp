#define _USE_MATH_DEFINES


#include <cmath>
#include <iostream>


#include "src/options_isolated.hpp"
#include "src/structure.hpp"
#include "src/write_data.hpp"


int main()
{
    bool verbose = false;
    #ifdef DETAILED_OUTPUT
        verbose = true;
    #endif

    // Parse options
    OptionsIsolated o("options_isolated");
    #ifdef DETAILED_OUTPUT  // Print all read options
      {
        std::cout << "**********\n";
        std::cout << "Parsed options:\n";
        o.print_options(true);
        std::cout << "**********\n";
      }
    #endif


    // Compute modifiers count per one lamella:
    float real_mmt_bead_d(o.real_r_c * 2*o.lj_bead_radius_clay);
    float real_mmt_area = M_PI * pow(o.platelet_radius * real_mmt_bead_d, 2);
    size_t charged_count;
    std::string modifier_count_taken_from;
    if (!o.modifiers_count_preset)
      {
        charged_count = 0.015 * real_mmt_area * o.stacking;
        modifier_count_taken_from = "CEC";
      }
    else
      {
        charged_count = o.modifiers_count_preset;
        modifier_count_taken_from = "options";
      }
    #ifdef DETAILED_OUTPUT  // Print information about modifier
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
        + 2 * real_mmt_bead_d * o.stacking;
    float min_width = 2*o.platelet_radius * real_mmt_bead_d
        * float(o.planar_expansion_coeff);
    float cube_edge = std::max(min_width, min_height) / o.real_r_c;
    #ifdef DETAILED_OUTPUT  // Print information about box
      {
        std::cout << "**********\n";
        std::cout << "Box computations:\n";
        std::cout << "min_height = " << min_height << std::endl;
        std::cout << "min_width = " << min_width << std::endl;
        std::cout << "cube_edge (in lj units) = " << cube_edge << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // Adjust polymers count:
    float lj_cell_volume = pow(cube_edge, 3);
    size_t all_beads_count = o.dpd_rho * lj_cell_volume;
    size_t mmt_beads_count = 2 * size_t(M_PI * pow(o.platelet_radius, 2));
    size_t single_mod_beads_count = 1 + o.tail_length;
    size_t all_mods_beads_count = single_mod_beads_count * charged_count;
    size_t single_polymer_beads_count = o.polymerization;
    size_t polymers_count = round((all_beads_count
        -mmt_beads_count - all_mods_beads_count) / o.polymerization);

    #ifdef DETAILED_OUTPUT  // Print information about polymers
      {
        std::cout << "**********\n";
        std::cout << "Polymers computations:\n";
        std::cout << "\tpolymers_count = " << polymers_count << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // *************************************************************************
    //                             Start main algorithm
    // *************************************************************************
    Structure s;

    // Set calcualted size to the box
    s.xlo = -cube_edge/2;
    s.xhi = cube_edge/2;
    s.ylo = -cube_edge/2;
    s.yhi = cube_edge/2;
    s.zlo = -cube_edge/2;
    s.zhi = cube_edge/2;

    // Add MMT stack
    size_t half_stacking = o.stacking / 2;
    size_t part_charged_count = charged_count / o.stacking;
    size_t charged_mmt_done = 0;
    float lj_interlayer = o.real_interlayer / o.real_r_c;
    float dz = lj_interlayer / 2 + 2*o.lj_bead_radius_clay;
    if (o.stacking % 2 != 0)
      {
        bool status = s.add_mmt_circular(o, 0, 0, 0, part_charged_count);
        charged_mmt_done += part_charged_count;
        dz = lj_interlayer + 4 * o.lj_bead_radius_clay;
      }
    for (size_t idx = 0; idx < half_stacking; ++idx)
      {
        size_t part_charged_count1 = part_charged_count;
        size_t part_charged_count2 = part_charged_count;
        if (idx == half_stacking - 1)
          {
            size_t charged_mmt_left = charged_count - charged_mmt_done;
            if (charged_mmt_left % 2 == 0)
              {
                part_charged_count1 = charged_mmt_left / 2;
                part_charged_count2 = charged_mmt_left / 2;
              }
            else
              {
                part_charged_count1 = charged_mmt_left / 2;
                part_charged_count2 = charged_mmt_left - part_charged_count1;
              }
          }
        bool status = s.add_mmt_circular(o, 0, 0, dz, part_charged_count1);
        charged_mmt_done += part_charged_count1;
        status = s.add_mmt_circular(o, 0, 0, -dz, part_charged_count2);
        dz += lj_interlayer + 4 * o.lj_bead_radius_clay;
        charged_mmt_done += part_charged_count2;
      }
    size_t mmt_atoms = s.atoms().size();
    #ifdef DETAILED_OUTPUT  // Print information about mmt addition status
      {
        std::cout << "**********\n";
        std::cout << "MMT addition:\n";
        std::cout << "atoms_count = " << s.atoms().size() << std::endl;
        std::cout << "\ncharged_mmt_done: " << charged_mmt_done;
        std::cout << "**********\n";
      }
    #endif
    #ifdef PARTIAL_DATAFILES  // Output incomplete data into file
      {
        std::string data_out_fname("incomplete_isolated_mmt");
        data_out_fname += "_r" + std::to_string(o.platelet_radius);
        data_out_fname += "_n" + std::to_string(o.stacking);
        data_out_fname += ".data";
        std::cout << "Writing incomplete data into " << data_out_fname
            << std::endl;
        write_data(data_out_fname, s);
      }
    #endif

    // Add modifiers
    std::vector<std::pair<float, float> > galleries;
    if (o.stacking % 2 == 0)
      {
        for (size_t idx = 0; idx < half_stacking; ++idx)
          {
            float top = lj_interlayer / 2
                + idx * (lj_interlayer + 4 * o.lj_bead_radius_clay);
            galleries.push_back(std::pair<float, float>(top - lj_interlayer, top));
          }
      }
    else
      {
        for (size_t idx = 0; idx < half_stacking; ++idx)
          {
            float top = lj_interlayer + 2 * o.lj_bead_radius_clay
                + idx * (lj_interlayer + 4 * o.lj_bead_radius_clay);
            galleries.push_back(std::pair<float, float>(top - lj_interlayer, top));
          }
      }
    size_t modifiers_done = 0;
    size_t modifiers_fails_done = 0;
    size_t modifiers_fails_allowed = charged_count * 10000;
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

    if (modifiers_done != charged_count)
      {
        std::cout << "Failed to add all modifiers required\n";
        return 0;
      }

    size_t modifier_atoms = s.atoms().size() - mmt_atoms;

    #ifdef DETAILED_OUTPUT  // Print overall information about modifier addition
      {
        std::cout << "**********\n";
        std::cout << "Modifier addition after all: "
            << "modifier_count_taken_from: " << modifier_count_taken_from
            <<  modifiers_done << " of " << charged_count
            << "; fais: " << modifiers_fails_done
            << " of " << modifiers_fails_allowed << "\n";
        std::cout << "**********\n";
      }
    #endif

    #ifdef PARTIAL_DATAFILES  // Output incomplete data into file
      {
        std::string data_out_fname("incomplete_isolated_mmt");
        data_out_fname += "_r" + std::to_string(o.platelet_radius);
        data_out_fname += "_n" + std::to_string(o.stacking);
        data_out_fname += "_mod_n" + std::to_string(modifiers_done);
        data_out_fname += "_tail" + std::to_string(o.tail_length);
        data_out_fname += ".data";
        std::cout << "Writing incomplete data into " << data_out_fname
            << std::endl;
        write_data(data_out_fname, s);
      }
    #endif

    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    size_t polymers_fails_allowed = polymers_count * o.polymerization;
    while (polymers_done < polymers_count
        && polymers_fails_done < polymers_fails_allowed)
      {
        #ifdef DETAILED_OUTPUT  // Print information about last step (sometimes)
        if (polymers_done % (size_t(polymers_count / 10)) == 0)
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
    if (polymers_done < polymers_count)
      {
        std::cout << "Failed to add all polymers required\n";
        return 0;
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

    // Write resulting datafile
    std::string data_out_fname("isolated_mmt");
    data_out_fname += "_r" + std::to_string(o.platelet_radius);
    data_out_fname += "_n" + std::to_string(o.stacking);
    data_out_fname += "_mod_n" + std::to_string(modifiers_done);
    data_out_fname += "_tail" + std::to_string(o.tail_length);
    data_out_fname += "_poly_p" + std::to_string(o.polymerization);
    data_out_fname += "_n" + std::to_string(polymers_done);
    data_out_fname += ".data";
    write_data(data_out_fname, s);

    #ifdef GENERAL_OUTPUT  // Print summary structural information
      {
        std::cout << "Structure created:"
            << "\n\tAtoms: " << s.atoms().size()
            << "\n\tBonds: " << s.bonds().size()
            << "\n\tCubic box side: " << cube_edge
            << std::endl;
        std::cout << "Atoms:"
            << "\n\tMMT:      1 - " << mmt_atoms
            << "\n\tmodifier: " << mmt_atoms + 1 << " - "
                << mmt_atoms + modifier_atoms
                << " (" << modifiers_done << " molecules)"
                << " (" << modifier_count_taken_from << ")"
            << "\n\tpolymer:  " << mmt_atoms + modifier_atoms + 1 << " - "
                << mmt_atoms + modifier_atoms + polymer_atoms
                << " (" << polymers_done << " molecules)" << std::endl;

        std::cout << "Resulting structure written into:\n\t"
            << data_out_fname << std::endl;

        float DPD_rho = float(s.atoms().size())
            / (s.xhi-s.xlo) / (s.yhi-s.ylo) / (s.zhi-s.zlo);
        float CEC(93 * 90*84/real_mmt_area * modifiers_done/108 / o.stacking);
        std::cout << "Physical parameters:\n"
            << "\tDPD_rho (numerical): " << DPD_rho
            << "\n\tCEC: " << CEC << std::endl;
      }
    #endif

    return 0;
}
