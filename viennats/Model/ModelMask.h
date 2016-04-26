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
#include <string>

namespace model {
///Add a layer to act as a mask

    class Mask {

        std::string fn;

        bool remove;

        bool surf;

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

        Mask(const std::string & Parameters)  {
            using namespace boost::spirit::classic;
            using namespace parser_actors;

            remove=true;
            surf=false;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("mask_file") >> '='  >>  '\"' >> *((~ch_p('\"'))[push_back_a(fn)]) >> '\"' >> ';') |
                            (str_p("surface_geometry")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(surf)]) >> ';') |
                            (str_p("remove_bottom")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(remove)]) >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");
        }
    };


}

#endif /* MODELMASK_H_ */
