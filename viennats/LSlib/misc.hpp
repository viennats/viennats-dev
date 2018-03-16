#ifndef MISC_HPP_
#define MISC_HPP_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <fstream>
#include <string>

//#include "kernel.hpp"
#include "math.hpp"



namespace lvlset {
    namespace misc {

        template <class LevelSetType>
        void PrintStatistics(const LevelSetType& l, const std::string& FileName) {
            std::ofstream f;
            if(!std::ifstream(FileName.c_str())) {
                f.open(FileName.c_str());
                f << "Number of Grid Points"    <<";";
                f << "Number of active Points"    <<";";
                f << "Overhead (number of ints)"  <<";";
                f << std::endl;
            } else {
                f.open(FileName.c_str(),std::ios_base::app);
            }
            f << l.num_pts()      <<";";
            f << l.num_active_pts()    <<";";
            f << "TODO"          <<std::endl;
            f.close();
        }

        template <class value_type> int GetStatusFromDistance(value_type value) {
            int x=static_cast<int>(value);
            if (value>=0.) {
                return (value<=static_cast<value_type>(x)+static_cast<value_type>(0.5))?x:x+1;
            } else {
                return (value>=static_cast<value_type>(x)-static_cast<value_type>(0.5))?x:x-1;
            }
        }

        template <class LevelSetType>
        std::string test(const LevelSetType& l) {

            //typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;
            typedef typename LevelSetType::value_type value_type;
            const int D=LevelSetType::dimensions;

            std::ostringstream oss;

            for (typename LevelSetType::template const_iterator_neighbor_filtered<typename LevelSetType::filter_all, 1> it(l);
                        !it.is_finished();it.next()) {

                assert(it.center().pt_id()!=LevelSetType::UNDEF_PT);            //TODO

                assert(it.center().pt_id()<LevelSetType::SEGMENT_PT);


                if (it.center().is_defined()) {
                    for (int i=0; i<2*D; ++i ) {
                        assert(it.neighbor(i).pt_id()<LevelSetType::SEGMENT_PT);
                        if (it.neighbor(i).is_defined()) {
                            if(math::abs(GetStatusFromDistance(it.center().value())-GetStatusFromDistance(it.neighbor(i).value()))>1)  {
                                oss << "The defined point " << it.center().start_indices() << " has an inconsistent defined neighbor in direction " << i << "!" << std::endl;
                                oss.precision(24);
                                oss << "Value center point: " << it.center().value() << "  Value neighbor point: " << it.neighbor(i).value() << std::endl;
                            }
                        } else {
                            if(it.neighbor(i).sign()==0)  {                     //TODO POS_SIGN
                                if (it.center().value()<value_type(-0.5)) {
                                    oss << "The defined point " << it.center().start_indices() << " has a level set value less than -0.5 but has an undefined positive neighbor in direction " << i << "!" << std::endl;
                                }
                            } else {
                                if (it.center().value()>value_type(0.5)) {
                                    oss << "The defined point " << it.center().start_indices() << " has a level set value greater than 0.5 but has an undefined negative neighbor in direction " << i << "!" << std::endl;
                                }
                            }
                        }
                    }
                } else {
                    for (int i=0; i<2*D; ++i ) {
                        assert(it.neighbor(i).pt_id()<LevelSetType::SEGMENT_PT);
                        if (!it.neighbor(i).is_defined()) {
                            if(it.center().sign()!=it.neighbor(i).sign()) {
                                oss << "The undefined run from " << it.center().start_indices() << " to " << it.center().end_indices() << " has undefined neighbor grid points of opposite sign in direction " << i << "!" << std::endl;
                            }
                        }
                    }
                }

            }

            return oss.str();
        }

        /*template <class GridTraitsType, class LevelSetTraitsType>
        void Write(const levelset<GridTraitsType, LevelSetTraitsType>& l, std::string FileName, int option, float limit) {
            //option = true : print only defined grid point
            //option = false : print also start of undefined runs

            typedef levelset<GridTraitsType, LevelSetTraitsType> LevelSetType;
            typedef typename LevelSetType::value_type value_type;
            typedef typename LevelSetType::size_type size_type;
            const int D=LevelSetType::Dimension;


            std::ofstream f(FileName.c_str());

            size_type num=0;

            for (typename LevelSetType::const_iterator_runs it(l);
                !it.IsAtEnd();it.GoToNextRun()) {
                if (option==0) if (!it.IsDefined()) continue;
                if (!l.Grid.IsPositionInfinite(it.start_coordinates())) num++;
                if (it.start_coordinates()!=it.end_coordinates()) if (!l.Grid.IsPositionInfinite(it.end_coordinates())) num++;
            }


            f<< "object 1 class array type float rank 1 shape " << D << " items "<< num <<" data follows" << std::endl;

            for (typename LevelSetType::const_iterator_runs it(l);
                    !it.IsAtEnd();it.GoToNextRun()) {
                if (option==0) if (!it.IsDefined()) continue;
                if (!l.Grid.IsPositionInfinite(it.start_coordinates())) {
                    for (int j=0;j<D;j++) f << (it.start_coordinates()[j]) << " ";
                    f << std::endl;
                }
                if (!l.Grid.IsPositionInfinite(it.end_coordinates())) {
                    if (it.start_coordinates()!=it.end_coordinates()) {
                        for (int j=0;j<D;j++) f << (it.end_coordinates()[j]) << " ";
                        f << std::endl;
                    }
                }

            }

            f << "object 2 class array type float rank 0 items "<< num<<" data follows" <<std::endl;
            for (typename LevelSetType::const_iterator_runs it(l);
                    !it.IsAtEnd();it.GoToNextRun()) {
                if (option==0) if (!it.IsDefined()) continue;

                float dist=static_cast<float>(it.value());
                dist=std::min(limit,dist);
                dist=std::max(-limit,dist);
                if (!l.Grid.IsPositionInfinite(it.start_coordinates())) f << dist << std::endl;
                if (it.start_coordinates()!=it.end_coordinates()) if (!l.Grid.IsPositionInfinite(it.end_coordinates())) f << dist << std::endl;
            }

            f << "attribute \"dep\" string \"positions\"" << std::endl;

            f << "object \"data\" class field" << std::endl;
            f << "component \"positions\" value 1" << std::endl;
            f << "component \"data\" value 2" << std::endl;
            f << "end" << std::endl;

            f.close();
        }*/


    }
}
#endif /*MISC_HPP_*/
