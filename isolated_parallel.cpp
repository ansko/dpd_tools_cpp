#define _USE_MATH_DEFINES


#include <cmath>

#include <array>
#include <iostream>
#include <thread>

#include "src/options_parser.hpp"
#include "src/write_data.hpp"


// Convert (idx_x, idx_y, idx_z) tuple into a single index
size_t three2one(size_t idx_x, size_t nx, size_t idx_y, size_t ny, size_t idx_z,
size_t nz)
{
    return idx_x * ny*nz + idx_y * nz + idx_z;
}


// Convert a single index into a (idx_x, idx_y, idx_z) tuple
std::array<size_t, 3> one2three(size_t idx1d, size_t nx, size_t ny, size_t nz)
{
    size_t idx_x, idx_y, idx_z;
    idx_x = idx1d / (ny*nz);
    idx1d -= idx_x * ny*nz;
    idx_y = idx1d / nz;
    idx_z = idx1d - idx_y * nz;
    return std::array<size_t, 3>({idx_x, idx_y, idx_z});
}


// Parallel polymers addition
void threading_function(size_t thread_idx, Structure &s,
    OptionsParser &o,
    size_t nx, size_t ny, size_t nz, size_t polymers_count)
{
    std::array<size_t, 3> idcs3d = one2three(thread_idx, nx, ny, nz);
    size_t polymerization(o.get<float>("polymerization"));
    size_t idx_x = idcs3d[0];
    size_t idx_y = idcs3d[1];
    size_t idx_z = idcs3d[2];
    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    size_t polymers_fails_allowed = polymers_count * polymerization;

    while (polymers_done < polymers_count
           && polymers_fails_done < polymers_fails_allowed)
      {
        #ifdef GENERAL_OUTPUT  // Information about thread's last step
        if (polymers_done % (size_t(polymers_count / 10)) == 0)
          {
            std::cout << thread_idx << " "
                << polymers_done << " " << polymers_count << std::endl;
          }
        #endif

        bool status = s.add_polymer_parallel(o, idx_x, nx, idx_y, ny, idx_z, nz);
        if (status)
          {
            polymers_done++;
          }
        else
          {
            polymers_fails_done ++;
          }
      }
}


int main()
{
    bool verbose = false;
    #ifdef DETAILED_OUTPUT
        verbose = true;
    #endif

    // Parse options
    OptionsParser o("options_isolated");
    #ifdef DETAILED_OUTPUT  // Print parsed options
      {
        std::cout << "**********\n";
        std::cout << "Parsed options:\n";
        o.print_options(true);
        std::cout << "**********\n";
      }
    #endif

    // Compute general parameters
    float platelet_closing(o.get<float>("platelet_closing"));
    float real_r_c(o.get<float>("real_r_c"));
    float lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay"));
    float dpd_rho(o.get<float>("dpd_rho"));
    float planar_expansion_coeff(o.get<float>("planar_expansion_coeff"));
    size_t platelet_radius(o.get<size_t>("platelet_radius"));
    size_t stacking(o.get<size_t>("stacking"));
    size_t modifiers_count_preset(o.get<size_t>("modifiers_count_preset"));
    size_t tail_length(o.get<size_t>("tail_length"));
    size_t polymerization(o.get<size_t>("polymerization"));
    size_t threads_nx(o.get<size_t>("threads_nx"));
    size_t threads_ny(o.get<size_t>("threads_ny"));
    size_t threads_nz(o.get<size_t>("threads_nz"));

    // Compute modifiers count per one lamella:
    float real_mmt_bead_d(real_r_c * 2*lj_bead_radius_clay * platelet_closing);
    float real_mmt_area = M_PI * pow(platelet_radius * real_mmt_bead_d, 2);

    size_t charged_count;
    std::string modifier_count_taken_from;
    if (!modifiers_count_preset)
      {
        charged_count = 0.015 * real_mmt_area * stacking;
        modifier_count_taken_from = "CEC";
      }
    else
      {
        charged_count = modifiers_count_preset;
        modifier_count_taken_from = "options";
      }

    #ifdef DETAILED_OUTPUT  // Print parameters of MMT and modifier charges
      {
        std::cout << "**********\n";
        std::cout << "Derived parameters:\n";
        std::cout << "real_mmt_area = " << real_mmt_area << std::endl;
        std::cout << "charged_count (per stack) = " << charged_count << std::endl;
        std::cout << "**********\n";
      }
    #endif

    // Compute interlayer based on modifiers size and count
    float lj_mmt_area = real_mmt_area / real_r_c / real_r_c;
    float lj_volume = (tail_length + 1) * charged_count / dpd_rho;
    float lj_interlayer = lj_volume / lj_mmt_area;
    float real_interlayer(lj_interlayer * real_r_c);

    // Compute box size
    float min_height = real_interlayer * stacking
        + 2 * real_mmt_bead_d * stacking;
    float min_width = 2*platelet_radius * real_mmt_bead_d
        * float(planar_expansion_coeff);
    float cube_edge = std::max(min_width, min_height) / real_r_c;

    #ifdef DETAILED_OUTPUT  // Print box infomation
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
    size_t all_beads_count = dpd_rho * lj_cell_volume;
    size_t mmt_beads_count = 2 * size_t(M_PI * pow(platelet_radius, 2));
    size_t single_mod_beads_count = 1 + tail_length;
    size_t all_mods_beads_count = single_mod_beads_count * charged_count;
    size_t single_polymer_beads_count = polymerization;
    size_t polymers_count = round((all_beads_count
        -mmt_beads_count - all_mods_beads_count) / polymerization);

    #ifdef DETAILED_OUTPUT  // Print polymer addition parameters
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
    s.xlo = -cube_edge/2;
    s.xhi = cube_edge/2;
    s.ylo = -cube_edge/2;
    s.yhi = cube_edge/2;
    s.zlo = -cube_edge/2;
    s.zhi = cube_edge/2;

    // Add MMT stack
    size_t half_stacking = stacking / 2;
    size_t part_charged_count = charged_count / stacking;
    size_t charged_mmt_done = 0;
    float dz = lj_interlayer / 2 + 2*lj_bead_radius_clay;
    if (stacking % 2 != 0)
      {
        bool status = s.add_mmt_circular(o, 0, 0, 0, part_charged_count);
        charged_mmt_done += part_charged_count;
        dz = lj_interlayer + 4 * lj_bead_radius_clay;
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
        dz += lj_interlayer + 4 * lj_bead_radius_clay;
        charged_mmt_done += part_charged_count2;
      }
    size_t mmt_atoms = s.atoms().size();

    #ifdef DETAILED_OUTPUT  // Print results of MMT addition
      {
        std::cout << "**********\n";
        std::cout << "MMT addition:\n";
        std::cout << "atoms_count = " << s.atoms().size() << std::endl;
        std::cout << "**********\n";
      }
    #endif

    #ifdef PARTIAL_DATAFILES
      {
        std::string data_out_fname("incomlete_parallel_isolated_mmt");
        data_out_fname += "_r" + std::to_string(platelet_radius);
        data_out_fname += "_n" + std::to_string(stacking);
        data_out_fname += ".data";
        write_data(data_out_fname, s);
      }
    #endif

    // Add modifiers
    std::vector<std::pair<float, float> > galleries;
    if (stacking % 2 == 0)
      {
        for (size_t idx = 0; idx < half_stacking; ++idx)
          {
            float top = lj_interlayer / 2
                + idx * (lj_interlayer + 4 * lj_bead_radius_clay);
            galleries.push_back(std::pair<float, float>(top - lj_interlayer, top));
          }
      }
    else
      {
        for (size_t idx = 0; idx < half_stacking; ++idx)
          {
            float top = lj_interlayer + 2 * lj_bead_radius_clay
                + idx * (lj_interlayer + 4 * lj_bead_radius_clay);
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
        std::cout << "Failed to add modifier\n";
        return 0;
      }

    size_t modifier_atoms = s.atoms().size() - mmt_atoms;

    #ifdef DETAILED_INFORMATION  // Print results of modifiers addition
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
        std::string data_out_fname("incomplete_parallel_isolated_mmt");
        data_out_fname += "_r" + std::to_string(platelet_radius);
        data_out_fname += "_n" + std::to_string(stacking);
        data_out_fname += "_mod_n" + std::to_string(modifiers_done);
        data_out_fname += "_tail" + std::to_string(tail_length);
        data_out_fname += ".data";
        write_data(data_out_fname, s);
      }
    #endif

    // Parallelism with polymers starts here
    size_t nx = threads_nx;
    size_t ny = threads_ny;
    size_t nz = threads_nz;

    std::vector<std::thread> my_threads;
    for (size_t idx_x = 0; idx_x < nx; ++idx_x)
        for (size_t idx_y = 0; idx_y < ny; ++idx_y)
            for (size_t idx_z = 0; idx_z < nz; ++idx_z)
              {
                size_t thread_idx(three2one(idx_x, nx, idx_y, ny, idx_z, nz));
                my_threads.emplace_back(
                    threading_function,
                    thread_idx,
                    std::ref(s), std::ref(o),
                    nx, ny, nz, size_t(polymers_count / (nx*ny*nz)));
              }
    for (auto &th : my_threads)
      {
        th.join();
      }

    // Parallel version add polymers by parts all of which are size_t,
    // so their sum may differ from the real value of polymers_count.
    // This part adds left chains in not-parallel regime,
    // since this count should be small.
    size_t polymer_atoms = s.atoms().size() - mmt_atoms - modifier_atoms;
    size_t polymers_done = polymer_atoms / polymerization;
    if (polymers_done < polymers_count)
      {
        std::cout << "Not all polymers were added parallely, correcting...\n";
        size_t polymers_fails_done(0);
        size_t polymers_fails_allowed((polymers_count - polymers_done) * 10);
        while (polymers_done < polymers_count
            && polymers_fails_done < polymers_fails_allowed)
          {
            #ifdef DETAILED_OUTPUT  // Print information about last step
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
      }

    #ifdef DETAILED_OUTPUT  // Print overall information about polymer addition
      {
        std::cout << "**********\n";
        std::cout << "Polymers addition after all: "
                  <<  polymers_done << " of " << polymers_count << "\n";
        std::cout << "**********\n";
      }
    #endif

    std::string data_out_fname("parallel_isolated_mmt");
    data_out_fname += "_r" + std::to_string(platelet_radius);
    data_out_fname += "_n" + std::to_string(stacking);
    data_out_fname += "_mod_n" + std::to_string(modifiers_done);
    data_out_fname += "_tail" + std::to_string(tail_length);
    data_out_fname += "_poly_p" + std::to_string(polymerization);
    data_out_fname += "_n" + std::to_string(polymers_done);
    data_out_fname += ".data";
    write_data(data_out_fname, s);

    #ifdef GENERAL_OUTPUT
      {
        std::cout << "Structure created:\n"
            << "\tAtoms: " << s.atoms().size() << std::endl
            << "\tBonds: " << s.bonds().size() << std::endl
            << "\tCubic box side: " << cube_edge << std::endl;
        std::cout << "Atoms:"
            << "\n\tMMT:      1 - " << mmt_atoms
            << "\n\tmodifier: " << mmt_atoms + 1 << " - "
                << mmt_atoms + modifier_atoms
                << " (" << modifiers_done << " molecules)"
                << " (" << modifier_count_taken_from << ")"
            << "\n\tpolymer:  " << mmt_atoms + modifier_atoms + 1 << " - "
                << s.atoms().size()
                << " (" << polymers_done << " molecules)" << std::endl;
        std::cout << "Resulting structure written into:\n\t"
            << data_out_fname << std::endl;
        std::cout << "Calculated lj interlayer = " << lj_interlayer << std::endl;
      }
    #endif

    #ifdef GENERAL_OUTPUT  // Print DPD_rho and CEC
      {
        float DPD_rho = float(s.atoms().size())
            / (s.xhi-s.xlo) / (s.yhi-s.ylo) / (s.zhi-s.zlo);
        float CEC(93 * 90*84/real_mmt_area * modifiers_done/108 / stacking);
        std::cout << "Physical parameters:\n"
            << "\tDPD_rho (numerical): " << DPD_rho
            << "\n\tCEC: " << CEC << std::endl;
      }
    #endif

    return 0;
}
