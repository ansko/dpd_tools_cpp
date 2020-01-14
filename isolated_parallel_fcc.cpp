#define _USE_MATH_DEFINES  // M_PI from cmath define support


#include <cmath>

#include <algorithm>
#include <array>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

#include "src/options_parser.hpp"
#include "src/write_data.hpp"


std::mutex mmm;

// Parallel polymers addition:
// performs splitting initial cell into the 3d array of small boxes
// which are independent from each other
void
par_poly_add(Structure &s, OptionsParser &o, BBox bbox, size_t polymers_count)
{
    size_t polymerization(o.get<float>("polymerization"));
    size_t polymers_done = 0;
    size_t polymers_fails_done = 0;
    size_t polymers_fails_allowed = polymers_count * polymerization;

    while (polymers_done < polymers_count
           && polymers_fails_done < polymers_fails_allowed)
      {
        AddPolymerParallelParameters parameters(o);
        bool status = s.add_poly_bbox(bbox, parameters);
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
    // Some option from MD
    float md_mmt_area(7500);
    float md_cec(92.6);
    float md_mods_count(108);

    // Parse options
    OptionsParser o("options_isolated");

    // Compute general parameters
    float platelet_closing(o.get<float>("platelet_closing"));
    float real_r_c(o.get<float>("real_r_c"));
    float lj_bead_radius_clay(o.get<float>("lj_bead_radius_clay"));
    float dpd_rho(o.get<float>("dpd_rho"));
    float planar_expansion_coeff(o.get<float>("planar_expansion_coeff"));
    size_t platelet_radius(o.get<size_t>("platelet_radius"));
    size_t stacking(o.get<size_t>("stacking"));
    //size_t modifiers_count_preset(o.get<size_t>("modifiers_count_preset"));
    size_t tail_length(o.get<size_t>("tail_length"));
    size_t polymerization(o.get<size_t>("polymerization"));
    size_t threads_nx(o.get<size_t>("threads_nx"));
    size_t threads_ny(o.get<size_t>("threads_ny"));
    size_t threads_nz(o.get<size_t>("threads_nz"));

    // Compute modifiers count per one lamella:
    float real_mmt_bead_d(real_r_c * 2 * lj_bead_radius_clay * platelet_closing);
    float real_mmt_area(M_PI * pow(platelet_radius * real_mmt_bead_d, 2));

    size_t charged_count;
    std::string modifier_count_taken_from;

    try
      {
        charged_count = o.get<size_t>("modifiers_count_preset");
        modifier_count_taken_from = "options";
      }
    catch(char const *ex)
      {
        charged_count = md_mods_count * (real_mmt_area * stacking / md_mmt_area);
        modifier_count_taken_from = "CEC";
      }

    // Compute interlayer based on modifiers size and count
    float mmt_R_lj(platelet_radius * 2 * lj_bead_radius_clay * platelet_closing);

    float lj_mmt_area = mmt_R_lj * mmt_R_lj * M_PI;
    float lj_all_mmt_area = lj_mmt_area * stacking;

    // This volume should be reserved for modifiers
    float lj_volume = (tail_length + 1) * charged_count / dpd_rho;

    float lj_interlayer = lj_volume / lj_all_mmt_area;

    // Periodicity of MMT (with intercalated modifier) along x-axis
    float mmt_z_period = lj_interlayer + 6 * lj_bead_radius_clay * platelet_closing;

    // Compute box size
    float min_height = mmt_z_period * (stacking + 1);
    float min_width =  2*platelet_radius * lj_bead_radius_clay * platelet_closing
        * float(planar_expansion_coeff);
    float cube_edge = std::max(min_width, min_height);

    // Adjust polymers count:
    float lj_cell_volume = pow(cube_edge, 3);
    size_t all_beads_count = dpd_rho * lj_cell_volume;
    size_t mmt_beads_count = 2 * size_t(M_PI * pow(platelet_radius, 2));
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
    s.xlo = -cube_edge/2;
    s.xhi = cube_edge/2;
    s.ylo = -cube_edge/2;
    s.yhi = cube_edge/2;
    s.zlo = -cube_edge/2;
    s.zhi = cube_edge/2;
    float lx(s.xhi - s.xlo);
    float ly(s.yhi - s.ylo);
    float lz(s.zhi - s.zlo);

    // Add MMT stack
    size_t half_stacking = stacking / 2;
    size_t part_charged_count = charged_count / stacking;
    size_t charged_mmt_done = 0;
    float dz = lj_interlayer / 2 + 2*lj_bead_radius_clay;
    if (stacking % 2 != 0)
      {
        AddMmtCircularParameters parameters(o, 0, 0, 0, part_charged_count);
        bool status = s.add_mmt_circular_fcc(parameters);
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
        AddMmtCircularParameters parameters_top(o, 0, 0, dz,
                                                part_charged_count1);
        bool status = s.add_mmt_circular_fcc(parameters_top);
        charged_mmt_done += part_charged_count1;
        AddMmtCircularParameters parameters_bot(o, 0, 0, -dz,
                                                part_charged_count2);
        status = s.add_mmt_circular_fcc(parameters_bot);
        dz += lj_interlayer + 4 * lj_bead_radius_clay;
        charged_mmt_done += part_charged_count2;
      }
    size_t mmt_atoms = s.atoms().size();

    // Add modifiers
    std::vector<std::pair<float, float> > galleries;
    if (stacking % 2 == 0)
      {
        for (size_t idx = 0; idx < half_stacking + 1; ++idx)
          {
            float top;
            float bottom;
            if (idx != half_stacking)
              {
                top = lj_interlayer / 2 + idx * mmt_z_period;
                bottom = top - lj_interlayer;
              }
            else
              {
                top = idx * mmt_z_period;
                bottom = top - lj_interlayer / 2;
              }
            galleries.push_back(std::pair<float, float>(bottom, top));
            if (idx != 0)
              {
                galleries.push_back(std::pair<float, float>(-top, -bottom));
              }
          }
      }
    else
      {
        for (size_t idx = 0; idx < half_stacking + 1; ++idx)
          {
            float top;
            float bottom;


            if (idx != half_stacking)
              {
                top = lj_interlayer + 2 * lj_bead_radius_clay;
                bottom = top - lj_interlayer;
              }
            else
              {
                top = lj_interlayer / 2 + 2 * lj_bead_radius_clay;
                bottom = top - lj_interlayer / 2;
              }
            galleries.push_back(std::pair<float, float>(bottom, top));
            galleries.push_back(std::pair<float, float>(-top, -bottom));
          }
      }

    size_t modifiers_done = 0;
    size_t modifiers_fails_done = 0;
    size_t modifiers_fails_allowed = charged_count * 10000;
    while (modifiers_done < charged_count
        && modifiers_fails_done < modifiers_fails_allowed)
      {
        // Chose gallery idx; since outer galleries are smaller,
        // the probability of chosnig them should be increased:

        std::vector<size_t> mod_idcs(galleries.size() + 4);
        std::vector<size_t> chosen;
        std::iota (std::begin(mod_idcs), std::end(mod_idcs), 0);
        mod_idcs[mod_idcs.size() - 4] = 0;
        mod_idcs[mod_idcs.size() - 3] = 0;
        mod_idcs[mod_idcs.size() - 2] = galleries.size() - 1;
        mod_idcs[mod_idcs.size() - 1] = galleries.size() - 1;

        std::sample(mod_idcs.begin(), mod_idcs.end(), std::back_inserter(chosen),
                1, std::mt19937{std::random_device{}()});
        size_t idx = mod_idcs[chosen[0]];

        float bottom = galleries[idx].first;
        float top = galleries[idx].second;
        AddModifierGalleryParameters parameters(o, top, bottom, "isolated");
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
        std::cout << "Failed to add all modifiers\n";
        return 0;
      }

    size_t modifier_atoms = s.atoms().size() - mmt_atoms;

    // Parallelism with polymers starts here
    size_t nx = threads_nx;
    size_t ny = threads_ny;
    size_t nz = threads_nz;
    size_t poly_count_sub(polymers_count / (nx * ny * nz));
    std::vector<std::thread> my_threads;
    for (size_t idx_x = 0; idx_x < nx; ++idx_x)
        for (size_t idx_y = 0; idx_y < ny; ++idx_y)
            for (size_t idx_z = 0; idx_z < nz; ++idx_z)
              {
                float cur_xlo(s.xlo + lx / nx * idx_x);
                float cur_xhi(s.xlo + lx / nx * (idx_x + 1));
                float cur_ylo(s.ylo + ly / ny * idx_y);
                float cur_yhi(s.ylo + ly / ny * (idx_y + 1));
                float cur_zlo(s.zlo + lz / nz * idx_z);
                float cur_zhi(s.zlo + lz / nz * (idx_z + 1));
                BBox bbox(cur_xlo, cur_xhi, cur_ylo, cur_yhi, cur_zlo, cur_zhi);
                my_threads.emplace_back(par_poly_add, std::ref(s), std::ref(o),
                                        bbox, poly_count_sub);
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
            AddPolymerParallelParameters parameters(o);
            BBox bbox(s.xlo, s.xhi, s.ylo, s.yhi, s.zlo, s.zhi);
            bool status = s.add_poly_bbox(bbox, parameters);
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

    if (polymers_done != polymers_count)
      {
        std::cout << "Failed to add all polymers!\n";
        return 0;
      }

    float CEC(md_cec * (md_mmt_area / real_mmt_area / stacking)
                     * (charged_count / md_mods_count));

    std::string data_out_fname("parallel_isolated_mmt");
    data_out_fname += "_r" + std::to_string(platelet_radius);
    data_out_fname += "_n" + std::to_string(stacking);
    data_out_fname += "_mod_n" + std::to_string(modifiers_done);
    data_out_fname += "_tail" + std::to_string(tail_length);
    data_out_fname += "_poly_p" + std::to_string(polymerization);
    data_out_fname += "_n" + std::to_string(polymers_done);
    data_out_fname += "_CEC" + std::to_string(int(CEC));
    data_out_fname += ".data";
    write_data(data_out_fname, s);

    std::cout << "Success! CEC = " << CEC << std::endl << data_out_fname
              << std::endl;

    return 0;
}
