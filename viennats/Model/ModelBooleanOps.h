/*
 * ModelBooleanOps.h
 *
 *  Created on: Sep 24, 2009
 *      Author: ertl
 */

#ifndef MODELBOOLEANOPS_H_
#define MODELBOOLEANOPS_H_

#include "../message.h"
#include "../parser_actors.h"
#include <string>
#include <boost/spirit/include/classic.hpp>


namespace model {
///Model for boolean operations between geometries

    class BooleanOps {

        std::string fn;

        int lvl;

        bool inv;

        bool remove;

        bool surf;

    public:

        const std::string & file_name() const {
            return fn;
        }

        const int level() const {
            return lvl;
        }

        const bool invert() const {
            return inv;
        }

        const bool remove_bottom() const {
            return remove;
        }

        const bool surface() const {
        	return surf;
        }

        BooleanOps(const std::string & Parameters):lvl(0),inv(false)  {
            using namespace boost::spirit::classic;
            using namespace parser_actors;

            remove=true;
            surf=false;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("geometry_file") >> '='  >>  '\"' >> *((~ch_p('\"'))[push_back_a(fn)]) >> '\"' >> ';') |
                            (str_p("level")  >> '='  >> int_p[assign_a(lvl)]  >> ';') |
                            (str_p("invert")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(inv)]) >> ';') |
                            (str_p("surface_geometry")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(surf)]) >> ';') |
                            (str_p("remove_bottom")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(remove)]) >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");
        }
    };


}

#endif /* MODELBOOLEANOPS_H_ */
