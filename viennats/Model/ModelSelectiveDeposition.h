/*
 * ModelSelectiveDeposition.h
 *
 *  Created on: Apr 17, 2019
 *      Author: toifl
 */

#ifndef MODELSELECTIVEDEPOSITION_H_
#define MODELSELECTIVEDEPOSITION_H_

#include <stack>
#include <vector>
#include <cassert>
#include <algorithm>

#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include "../parser_actors.h"
#include "vector.hpp"
#include "Math.h"


namespace model {

///Anisotropic wet etching model

    class SelectiveDeposition {

        lvlset::vec<double,3> direction100;
        lvlset::vec<double,3> direction010;

        //the rates have to be given in the layer order, from bottom to topmost layer
        double r100;
        double r110;
        double r111;
        double r311;

        std::vector<int> depo_possible; //TODO change to boolean

    public:

        static const bool OutputFluxes=false;
        static const bool SpatiallyEqualDistributedFlux=true;
        static const bool ReemissionIsMaterialDependent=true;	//TODO is set to true for output
        static const bool CalculateConnectivities=false;        //TODO can be useful?
        static const bool CalculateVisibilities=false;
        static const bool CalculateNormalVectors=true;


        static const int CoverageStorageSize=0;
        static const int RatesStorageSize=0;
        static const unsigned int NumberOfParticleTypes=0;
        unsigned int NumberOfParticleClusters[1];

        class ParticleType {
        public:
            double Direction[3];
            double Flux;
        };

        SelectiveDeposition(const std::string & Parameters, const int AddLayer){
            using namespace boost::spirit::classic;
            using namespace parser_actors;


            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction100")  >> '='  >> '{' >> real_p[assign_a(direction100[0])]  >> "," >> real_p[assign_a(direction100[1])] >> "," >> real_p[assign_a(direction100[2])] >> '}' >> ';') |
                            (str_p("direction010")  >> '='  >> '{' >> real_p[assign_a(direction010[0])]  >> "," >> real_p[assign_a(direction010[1])] >> "," >> real_p[assign_a(direction010[2])] >> '}' >> ';') |
                            (str_p("rates100")  >> '='  >> real_p[assign_a(r100)]  >>  ';') |
                            (str_p("rates110")  >> '='  >> real_p[assign_a(r110)]  >>  ';') |
                            (str_p("rates111")  >> '='  >> real_p[assign_a(r111)]  >>  ';') |
                            (str_p("rates311")  >> '='  >> real_p[assign_a(r311)]  >>  ';') |
                            (str_p("depo_possible")  >>  '='  >>  '{' >> ( int_p[push_back_a(depo_possible)]  % ',')>> '}'  >> ';') //TODO change to boolean parsing

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");


            // each new added layer must be depositable
            for(int i=0; i<AddLayer; ++i) depo_possible.insert(depo_possible.begin(), 1);
        }

        template<class VecType>
        void CalculateVelocity(
                double &Velocity,
                const VecType& NormalVector,
                const double *Coverages,
                const double *Rates,
                int matnum,
                bool connected,
                bool visible) const {

            assert(matnum < (int) depo_possible.size() );

            if (depo_possible[matnum] == 0) {
                Velocity=0;
                return;
            }

            lvlset::vec<double,3> nv{NormalVector[0],NormalVector[1],NormalVector[2]};

            Velocity = my::math::fourRateInterpolation<double,3>(nv, direction100, direction010, r100, r110, r111, r311);

        }

		template<class VecType>
		void CalculateVectorVelocity(
				double *Velocity,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible) const {}


        static void UpdateCoverage(double *Coverages, const double *Rates) {}

        template <class PT> static void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) {}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {}

		template <class PT, class VecType> static void ParticleReflexion(
							const PT& p,
							std::stack<PT>& particle_stack,
							const VecType& NormalVector,
							const double* Coverages,
							int Material//,
//                            int D,
//                            double dot // dot product between the incoming particle direction and the normal vector
                            ) {}

      const std::vector<int>& get_depo_possible() const{
        return depo_possible;
      }
    };

//    const unsigned int WetEtching::NumberOfParticleClusters[ConstantRates::NumberOfParticleTypes]={};

}
#endif /* MODELSELECTIVEDEPOSITION_H_ */
