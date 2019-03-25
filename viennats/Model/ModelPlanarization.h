/*
 * ModelPlanarization.h
 *
 *  Created on: Apr 20, 2009
 *      Author: ertl
 */

#ifndef MODELPLANARIZATION_H_
#define MODELPLANARIZATION_H_

#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include "../parser_actors.h"

namespace model {
///Chemical mechanical planarization

    class Planarization {

        double Coordinate;
        bool Fill;

    public:

        bool fill_up() const {
            return Fill;
        }

        double get_coordinate() const {
            return Coordinate;
        }



        Planarization(const std::string & Parameters) : Fill(false) {
            using namespace boost::spirit::classic;
            using namespace parser_actors;

            Fill=false;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("coordinate")  >> '='  >>  real_p[assign_a(Coordinate)] >>  ';') |
                            (str_p("fill_up")  >>  '='  >> (str_p("true") | str_p("false"))[assign_bool(Fill)] >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");
        }
    };


}

#endif /* MODELPLANARIZATION_H_ */
