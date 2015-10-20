#ifndef PARSER_ACTORS_H_
#define PARSER_ACTORS_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <string>

namespace parser_actors {

    class assign_bool {
        bool & bo;
    public:

        assign_bool(bool& b2) : bo(b2) {}

        template <class iter>
        void operator()(const iter&  a, const iter& b ) const {
            bo=(std::string(a,b)=="true");

        }
    };



    class assign_dir {
        int & dir;
        bool & sign;
    public:

        assign_dir(int & dir2, bool& sign2) : dir(dir2), sign(sign2) {}

        template <class iter>
        void operator()(const iter&  a, const iter& b ) const {
            dir=a[1]-'x';
            sign=(a[0]=='-');
        }
    };

    class assign_input_transformation {
        std::vector<int> & dir;
        std::vector<bool> & signs;
    public:

        assign_input_transformation(std::vector<int> & dir2, std::vector<bool> & sign2) : dir(dir2), signs(sign2) {}

        template <class iter>
        void operator()(const iter&  a, const iter& b ) const {
            dir.push_back(a[1]-'x');
            signs.push_back(a[0]=='-');
        }


    };

    template <class C>
    class assign_enum {
        C& a;
        const C b;
    public:
        assign_enum(C& x, const C B): a(x), b(B) {}

        template <class iter>
        void operator()(const iter&  v, const iter& w ) const  {
            a=b;
        }
    };



}

#endif /* PARSER_ACTORS_H_ */
