/*
 * ModelMask.h
 *
 *  Created on: Apr 21, 2009
 *      Author: ertl
 */

#ifndef MODELMASK_H_
#define MODELMASK_H_

#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include "../parser_actors.h"
#include <string>

namespace model {
///Add a layer to act as a mask

    class Mask {

        std::string fn;

        bool remove;

        bool surf;

        bool ignore_other_mat;

    public:

        const std::string & file_name() const {
            return fn;
        }

         const bool remove_bottom() const {
            return remove;
        }

         const bool surface() const {
            return surf;
        }

        const bool ignore_other_materials() const {
            return ignore_other_mat;
        }

        Mask(const std::string & Parameters)  {
            using namespace boost::spirit::classic;
            using namespace parser_actors;

            remove=false;
            surf=false;
            ignore_other_mat = false;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("mask_file") >> '='  >>  '\"' >> *((~ch_p('\"'))[push_back_a(fn)]) >> '\"' >> ';') |
                            (str_p("surface_geometry")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(surf)]) >> ';') |
                            (str_p("remove_bottom")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(remove)]) >> ';') |
                            (str_p("ignore_other_materials")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(ignore_other_mat)]) >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if(!b) msg::print_error("Failed interpreting model parameters!");
        }
    };


}

#endif /* MODELMASK_H_ */
