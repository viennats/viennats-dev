#ifndef MODELCONSTANTRATES_H_
#define MODELCONSTANTRATES_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <stack>
#include <vector>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include "../parser_actors.h"
#include <iostream>

using namespace std;

namespace model {
///Model which applies a constant isotropic or anisotropic velocity

	class ConstantRates {

        double StartDirection[3];
        std::vector<double> isotropic_rates;
        std::vector<double> directional_rates;
        std::vector<double> constant_rates;
        std::vector<double> vector_rates;
				bool bidirectionalRates;


	public:

        bool MaskLayer;
        bool DirEmpty;
        static const bool OutputFluxes=false;
        static const bool SpatiallyEqualDistributedFlux=true;
	    static const bool ReemissionIsMaterialDependent=false;
	    bool CalculateConnectivities;
	    bool CalculateVisibilities;
	    bool CalculateNormalVectors;
	    bool IncludeVectorRates;

		class ParticleType {
		public:
			double Direction[3];
			double Flux;
		};

		ConstantRates(const std::string & Parameters, bool masked=false) : CalculateConnectivities(true), CalculateVisibilities(true),CalculateNormalVectors(true),IncludeVectorRates(true) {
				using namespace boost::spirit::classic;
				using namespace parser_actors;

				bidirectionalRates=false;

		    bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("constant_rates")  >> '='  >>  '{' >> (real_p[push_back_a(constant_rates)] % ',') >> '}'  >> ';') |
                            (str_p("isotropic_rates")  >> '='  >>  '{' >> (real_p[push_back_a(isotropic_rates)] % ',') >> '}'  >> ';') |
                            (str_p("directional_rates")  >> '='  >>  '{' >> (real_p[push_back_a(directional_rates)] % ',') >> '}'  >> ';') |
                            (str_p("vector_rates")  >> '='  >>  '{' >> (real_p[push_back_a(vector_rates)] % ',') >> '}'  >> ';') |
														(str_p("bidirectional_rates")  >> '='  >> ((str_p("true") | str_p("false"))[assign_bool(bidirectionalRates)]) >> ';')

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (isotropic_rates.size()==0) CalculateConnectivities=false;

            if (directional_rates.size()==0) {
                CalculateVisibilities=false;
                if (vector_rates.size()==0) {
                	CalculateNormalVectors=false;
                }
            }

            if (vector_rates.size()==0) IncludeVectorRates=false;

            DirEmpty=true;
            for(int i=0;i<3;i++) if(StartDirection[i]!=0) DirEmpty=false;
            if (!DirEmpty) {
            	double start_norm=0;
            	for (int i=0;i<3;++i) start_norm+=StartDirection[i]*StartDirection[i];
            	start_norm=std::sqrt(start_norm);
            	for (int i=0;i<3;++i) StartDirection[i]/=start_norm;
            }// else DirEmpty=0;

            MaskLayer=masked;
            if (!b) msg::print_error("Failed interpreting process parameters!");

            // ---------------------------
            //CalculateNormalVectors=true;
            // ---------------------------

		}

		static const int CoverageStorageSize=0;
		static const int RatesStorageSize=0;
		static const unsigned int NumberOfParticleTypes=0;
//		static const unsigned int NumberOfParticleClusters[NumberOfParticleTypes];
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

		    double isotropic_rate=(Material < static_cast<int>(isotropic_rates.size()))?isotropic_rates[Material]:0;
		    double directional_rate=(Material < static_cast<int>(directional_rates.size()))?directional_rates[Material]:0;
		    double constant_rate=(Material < static_cast<int>(constant_rates.size()))?constant_rates[Material]:0;

		    Velocity=constant_rate;

		    if (connected) Velocity+=isotropic_rate;

		    if ((visible) && (directional_rate!=0)) {
		        double dot=0.;
		        for (int i=0;i<3;++i) dot-=StartDirection[i]*NormalVector[i];
		        if(bidirectionalRates) Velocity+=directional_rate*std::abs(dot);
						else Velocity+=directional_rate*std::max(0.,dot);
		    }
		}

		template<class VecType>
		void CalculateVectorVelocity(
				double *Velocity,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible
				) const {}

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

//	const unsigned int ConstantRates::NumberOfParticleClusters[ConstantRates::NumberOfParticleTypes]={};

}

#endif /*MODELCONSTANTRATES_H_*/
