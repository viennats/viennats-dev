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

        lvlset::vec<double,3> directions[3];
        //lvlset::vec<double,3> direction010;
        //lvlset::vec<double,3> direction001;
        double r100;
        double r110;
        double r111;
        double r311;
        bool has_mask;

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

        WetEtching(const std::string & Parameters) : has_mask(false) {
            using namespace boost::spirit::classic;
            using namespace parser_actors;


            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction100")  >> '='  >> '{' >> real_p[assign_a(directions[0][0])]  >> "," >> real_p[assign_a(directions[0][1])] >> "," >> real_p[assign_a(directions[0][2])] >> '}' >> ';') |
                            (str_p("direction010")  >> '='  >> '{' >> real_p[assign_a(directions[1][0])]  >> "," >> real_p[assign_a(directions[1][1])] >> "," >> real_p[assign_a(directions[1][2])] >> '}' >> ';') |
                            (str_p("r100")  >> '='  >> real_p[assign_a(r100)]  >>  ';') |
                            (str_p("r110")  >> '='  >> real_p[assign_a(r110)]  >>  ';') |
                            (str_p("r111")  >> '='  >> real_p[assign_a(r111)]  >>  ';') |
                            (str_p("r311")  >> '='  >> real_p[assign_a(r311)]  >>  ';') |
                            (str_p("has_mask")  >>  '='  >> (str_p("true") | str_p("false"))[assign_bool(has_mask)] >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

            directions[0]=Normalize(directions[0]);
            directions[1]=Normalize(directions[1]-directions[0]*dot(directions[0],directions[1]));
            directions[2]=cross(directions[0], directions[1]);

            for (int i=0;i<3;++i) assert(dot(directions[i], directions[(i+1)%3])<1e-6);

        }

        static const int CoverageStorageSize=0;
        static const int RatesStorageSize=0;
        static const unsigned int NumberOfParticleTypes=0;
//        static const unsigned int NumberOfParticleClusters[NumberOfParticleTypes];
        static const unsigned int NumberOfParticleClusters[1];

        template<class VecType>
        void CalculateVelocity(
                double &Velocity,
                const VecType& NormalVector,
                const double *Coverages,
                const double *Rates,
                int Material,
                bool connected,
                bool visible) const {

            if (has_mask && Material>0) {
                Velocity=0;
                return;
            }

            lvlset::vec<double,3> N;

            for (int i=0;i<3;i++) N[i]=std::fabs(directions[i][0]*NormalVector[0]+directions[i][1]*NormalVector[1]+directions[i][2]*NormalVector[2]);
            N.reverse_sort();

            assert(std::fabs(Norm(N)-1)<1e-4);

            if (dot(N, lvlset::vec<double,3>(-1,1,2))<0) {
                Velocity=-((r100*(N[0]-N[1]-2*N[2])+r110*(N[1]-N[2])+3*r311*N[2])/N[0]);    //region A
            } else {
                Velocity=-((r111*((N[1]-N[0])*0.5+N[2])+r110*(N[1]-N[2])+1.5*r311*(N[0]-N[1]))/N[0]);//region C
            }

           /* double isotropic_rate=(Material < static_cast<int>(isotropic_rates.size()))?isotropic_rates[Material]:0;
            double directional_rate=(Material < static_cast<int>(directional_rates.size()))?directional_rates[Material]:0;
            double constant_rate=(Material < static_cast<int>(constant_rates.size()))?constant_rates[Material]:0;

            Velocity=constant_rate;

            if (connected) Velocity+=isotropic_rate;

            if (visible) {
                double dot=0.;
                for (int i=0;i<3;++i) dot-=StartDirection[i]*NormalVector[i];
                Velocity+=directional_rate*std::max(0.,dot);
            }*/

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
                            int Material,
                            int D,
                            double dot // dot product between the incoming particle direction and the normal vector
                            ) {}
    };

//    const unsigned int WetEtching::NumberOfParticleClusters[ConstantRates::NumberOfParticleTypes]={};

}
#endif /* MODELWETETCHING_H_ */
