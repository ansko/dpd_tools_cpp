#ifndef OPTIONS_PARSER_HPP
#define OPTIONS_PARSER_HPP


#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>


// Read options from a specified file, where they are stored as
// OPTION_NAME VALUE
// and store read values in a map<option_name, option_value>.
class OptionsParser
{
public:
    OptionsParser(std::string options_fname)
      {
        std::string option_name;
        std::string buffer_line;
        std::ifstream ofs(options_fname);
        while (std::getline(ofs, buffer_line))
          {
            if (buffer_line.find("#") == std::string::npos)
              {
                std::stringstream ss(buffer_line);
                ss >> option_name >> this->options_[option_name];
              }
          }
      };

    void print_options(bool verbose=false)
      {
        for (auto it = this->options_.begin(); it != this->options_.end(); ++it)
          {
            std::cout << it->first << " " << it->second << std::endl;
          }
      }

    template <typename T> T get(std::string option_name)
      {
        return T(this->options_[option_name]);
      };

private:
    std::map<std::string, float> options_;
};


#endif  // OPTIONS_PARSER_HPP
