#ifndef OPTIONS_ISOLATED_HPP
#define OPTIONS_ISOLATED_HPP


#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "options_common.hpp"


class OptionsIsolated : public OptionsCommon
{
public:
    size_t platelet_radius = 0;
    size_t stacking = 0;
    float planar_expansion_coeff = -1;

    OptionsIsolated(std::string options_fname)
        : OptionsCommon(options_fname)
      {
        std::string option_name;
        std::ifstream ofs(options_fname);
        while (ofs >> option_name)
          {
            if (option_name == "planar_expansion_coeff")
                ofs >> this->planar_expansion_coeff;
            else if (option_name == "dpd_rho")
                ofs >> this->dpd_rho;
            else if (option_name == "platelet_radius")
                ofs >> this->platelet_radius;
            else if (option_name == "stacking")
                ofs >> this->stacking;
          }
      };

    void print_options(bool verbose=false)
      {
        OptionsCommon::print_options(verbose);

        std::cout << "planar_expansion_coeff = " << this->planar_expansion_coeff;
        if (verbose)
            std::cout << " Cell increase coeff in xy plane around MMT platelet\n";
        else
            std::cout << std::endl;
        std::cout << "platelet_radius = " << this->platelet_radius;
        if (verbose)
            std::cout << " Radius of MMT platelet (in beads)\n";
        else
            std::cout << std::endl;
        std::cout << "stacking = " << this->stacking;
        if (verbose)
            std::cout << " MMT platelets count in the stack\n";
        else
            std::cout << std::endl;
      }
};


#endif  // OPTIONS_ISOLATED_HPP
