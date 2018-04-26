#ifndef MESSAGE_H_
#define MESSAGE_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <string>
#include <ostream>

///Namespace for all stdout related output to the user
namespace msg {

    const int width=60;

    void print_welcome() {
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "#==========================================================================#" << std::endl;
        std::cout << "#   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.   #" << std::endl;
        std::cout << "#                                                                          #" << std::endl;
        std::cout << "#                            -----------------                             #" << std::endl;
        std::cout << "#                ViennaTS - The Vienna Topography Simulator                #" << std::endl;
        std::cout << "#                            -----------------                             #" << std::endl;
        std::cout << "#                                                                          #" << std::endl;
        std::cout << "#    Contact:         viennats@iue.tuwien.ac.at                            #" << std::endl;
        std::cout << "#                                                                          #" << std::endl;
        std::cout << "#    License:         MIT (X11), see file LICENSE in the base directory    #" << std::endl;
        std::cout << "#==========================================================================#" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }

    void print_help_extended() {
      std::cout << "USAGE:" << std::endl;
      std::cout << "SIMULATION: viennats parameterfile" << std::endl;
      std::cout << "CONVERSION: viennats [-ls2vtk | -ls2dx | -p | -p2o] FILES" << std::endl;
      std::cout << "Options:" << std::endl;
      std::cout << "    -ls2vtk    Convert LVST file(s) to VTK file(s)." << std::endl;
      std::cout << "    -ls2dx     Convert LVST file(s) to DX file(s)." << std::endl;
      std::cout << "    -p         Print LVST file(s) to the console." << std::endl;
      std::cout << "    -p2o       Print LVST file(s) to TXT file(s)." << std::endl;
      std::cout << "NOTE: File paths and names are taken from FILES." << std::endl;
      std::cout << "Other:" << std::endl;
      std::cout << "   --help      Shows (extended) help." << std::endl;
      std::cout << "   --version   Shows version." << std::endl;

    }

    void print_help() {
      std::cout << "USAGE: viennats [option [files...]| parameterfile]" << std::endl;
      std::cout << "Try viennats --help for more information." << std::endl;
    }

    void print_version() {
      std::cout << "viennats " << VIENNATS_VERSION << std::endl;
      std::cout << "Copyright (C) 2008-2015  Institute for Microelectronics, TU Wien" << std::endl << std::endl;
      std::cout << "Contact: viennats@iue.tuwien.ac.at" << std::endl;
      std::cout << "GitHub: https://github.com/viennats/viennats-dev" << std::endl;
      std::cout << "License: MIT (X11), see file LICENSE in the base directory" << std::endl;
      std::cout << "The software is provided \"AS IS\", without warranty of any kind," << std::endl
                << "express or implied, including but not limited to the Warranies of" << std::endl
                << "merchantability, fitness for a particular purpose and noninfringement." << std::endl;
    }

    void print_message(const std::string& s) {
        std::cout << s << std::endl << std::endl;
    }

    void print_message_2(const std::string& s) {
        std::cout << s << std::endl;
    }

    void print_message(const std::string& s, double v) {
        std::cout << s << v << std::endl << std::endl;
    }

    void print_message(const std::string& s, long double v) {
        std::cout << s << v << std::endl << std::endl;
    }

    void print_message(const std::string& s, int v) {
        std::cout << s << v << std::endl << std::endl;
    }

    void print_start(const std::string& s) {
        std::cout << std::setw(40) << std::left << s << std::flush;
    }

    void print_start(const std::string& s, int i) {
        std::cout << std::setw(40) << std::left << s << " surface " << i << std::flush << std::endl;
    }

    void print_done() {
        std::cout << "    [done]" << std::endl << std::endl;
    }

    void print_warning(const std::string& s) {
        std::cout << "    WARNING: " << s << std::endl << std::endl;
    }

    void print_error(const std::string& s) {
        std::cout << "    ERROR: " << s << std::endl << std::endl;
        abort();
    }


}


#endif /* MESSAGE_H_ */
