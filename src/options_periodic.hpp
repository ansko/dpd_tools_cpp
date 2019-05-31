#ifndef OPTIONS_PERIODIC_HPP
#define OPTIONS_PERIODIC_HPP


#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "options_common.hpp"

class OptionsPeriodic : public OptionsCommon
{
public:
    size_t platelet_edge = 0;

    OptionsPeriodic(std::string options_fname)
        : OptionsCommon(options_fname)
      {
        std::string option_name;
        std::ifstream ofs(options_fname);
        while (ofs >> option_name)
          {
            if (option_name == "platelet_edge")
                ofs >> this->platelet_edge;
          }
      };

    void print_options(bool verbose=false)
      {
        OptionsCommon::print_options(verbose);

        std::cout << "platelet_edge = " << this->platelet_edge;
        if (verbose)
            std::cout << " MMT platelet edge in beads\n";
        else
            std::cout << std::endl;
      }
};


#endif  // OPTIONS_PERIODIC_HPP
