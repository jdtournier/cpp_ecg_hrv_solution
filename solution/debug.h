#pragma once

#include <iostream>
#include <string>

#define VAR(x) std::cerr << __FILE__ << ", in " << __func__ << "(), line " << __LINE__  << ": " << #x << " = " << x << "\n"

namespace debug {

  inline bool verbose_mode = false;

  void set_verbose_mode (int& argc, char* argv[]);

  inline void log (const std::string& text)
  {
    if (verbose_mode)
      std::cerr << "[DEBUG] " << text << "\n";
  }

}
