#ifndef OPTIONS_ISOLATED_PARALLEL_HPP
#define OPTIONS_ISOLATED_PARALLEL_HPP


#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "options_isolated.hpp"


class OptionsIsolatedParallel : public OptionsIsolated
{
public:
    size_t threads_nx = 0;
    size_t threads_ny = 0;
    size_t threads_nz = 0;

    OptionsIsolatedParallel(std::string options_fname)
        : OptionsIsolated(options_fname)
      {
        std::string option_name;
        std::ifstream ofs(options_fname);
        while (ofs >> option_name)
          {
            if (option_name == "threads_nx")
                ofs >> this->threads_nx;
            else if (option_name == "threads_ny")
                ofs >> this->threads_ny;
            else if (option_name == "threads_nz")
                ofs >> this->threads_nz;
          }
      };

    void print_options(bool verbose=false)
      {
        OptionsIsolated::print_options(verbose);

        std::cout << "threads_nx = " << this->threads_nx;
        if (verbose)
            std::cout << " Threads along x\n";
        else
            std::cout << std::endl;
        std::cout << "threads_ny = " << this->threads_ny;
        if (verbose)
            std::cout << " Threads along y\n";
        else
            std::cout << std::endl;
        std::cout << "threads_nz = " << this->threads_nz;
        if (verbose)
            std::cout << " Threads along z\n";
        else
            std::cout << std::endl;
      }
};


#endif  // OPTIONS_ISOLATED_PARALLEL_HPP
