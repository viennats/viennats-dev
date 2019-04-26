/*
 * ModelWetEtching.h
 *
 *  Created on: Apr 22, 2009
 *      Author: ertl
 */

#ifndef MODELWETETCHING_H_
#define MODELWETETCHING_H_

#include <stack>
#include <vector>
#include <cassert>
//#include <algorithm>
//#include <functional>

#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include "../parser_actors.h"
#include "vector.hpp"

namespace model {

///Anisotropic wet etching model

    class WetEtching {

        const double EPS = 1e-6;

        lvlset::vec<double,3> direction100;
        lvlset::vec<double,3> direction010;



        std::vector<double> r100;
        std::vector<double> r110;
        std::vector<double> r111;
        std::vector<double> r311;

        std::vector<bool> zeroVel;


    public:

        static const bool OutputFluxes=false;
        static const bool SpatiallyEqualDistributedFlux=true;
        static const bool ReemissionIsMaterialDependent=true;	//TODO is set to true for output
        static const bool CalculateConnectivities=false;        //TODO can be useful?
        static const bool CalculateVisibilities=false;
        static const bool CalculateNormalVectors=true;

        class ParticleType {
        public:
            double Direction[3];
            double Flux;
        };

        WetEtching(const std::string & Parameters) {
            using namespace boost::spirit::classic;
            using namespace parser_actors;


            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                          (str_p("direction100")  >> '='  >> '{' >> real_p[assign_a(direction100[0])]  >> "," >> real_p[assign_a(direction100[1])] >> "," >> real_p[assign_a(direction100[2])] >> '}' >> ';') |
                          (str_p("direction010")  >> '='  >> '{' >> real_p[assign_a(direction010[0])]  >> "," >> real_p[assign_a(direction010[1])] >> "," >> real_p[assign_a(direction010[2])] >> '}' >> ';') |
                          (str_p("rates100")  >>  '='  >>  '{' >> ( real_p[push_back_a(r100)]  % ',')>> '}'  >> ';') |
                          (str_p("rates110")  >>  '='  >>  '{' >> ( real_p[push_back_a(r110)]  % ',')>> '}'  >> ';') |
                          (str_p("rates111")  >>  '='  >>  '{' >> ( real_p[push_back_a(r111)]  % ',')>> '}'  >> ';') |
                          (str_p("rates311")  >>  '='  >>  '{' >> ( real_p[push_back_a(r311)]  % ',')>> '}'  >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");


            //reverse because material id is reversed w.r.t. layer id (for output)
            /*std::reverse(r100.begin(),r100.end());
            std::reverse(r110.begin(),r110.end());
            std::reverse(r111.begin(),r111.end());
            std::reverse(r311.begin(),r311.end());*/

            // find materials with no growth in any direction and store in zeroVel
            for(unsigned int i=0; i < r100.size(); ++i){
              zeroVel.push_back(false);
              if(fabs(r100[i]) < EPS)
                if(fabs(r110[i]) < EPS)
                  if(fabs(r111[i]) < EPS)
                    if(fabs(r311[i]) < EPS)
                      zeroVel[i]=true;
            }


        }

        static const int CoverageStorageSize=0;
        static const int RatesStorageSize=0;
        static const unsigned int NumberOfParticleTypes=0;
        unsigned int NumberOfParticleClusters[1];

        template<class VecType>
        void CalculateVelocity(
                double &Velocity,
                const VecType& NormalVector,
                const double *Coverages,
                const double *Rates,
                int Material,
                bool connected,
                bool visible) const {


            if (zeroVel[Material] == true) {
                Velocity=0;
                return;
            }

            lvlset::vec<double,3> nv{NormalVector[0],NormalVector[1],NormalVector[2]};

            Velocity = my::math::fourRateInterpolation<double,3>(nv, direction100, direction010, r100[Material], r110[Material], r111[Material], r311[Material]);


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
    };

//    const unsigned int WetEtching::NumberOfParticleClusters[ConstantRates::NumberOfParticleTypes]={};

}
#endif /* MODELWETETCHING_H_ */
