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
      std::cout << "USAGE: viennats OPTION FILE" << std::endl;
      std::cout << "Convert:" << std::endl;
      std::cout << "   -lvst2vtk[N]   Convert LVST FILE of dimension N to VTK FILE." << std::endl;
      std::cout << "Other:" << std::endl;
      std::cout << "  --help          Shows extended help." << std::endl;
      std::cout << "   -h             Shows help." << std::endl;
      std::cout << "NOTE:If no OPTION was entered FILE will be used as a parameter file for ViennaTS." << std::endl;
    }

    void print_help() {
      std::cout << "USAGE: viennats [OPTION] [FILE]" << std::endl;
      std::cout << "Try viennats --help for more information." << std::endl;
    }

    void print_message_2(const std::string& s) {
        std::cout << s << std::endl;
    }

    void print_message(const std::string& s) {
        std::cout << s << std::endl << std::endl;
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
