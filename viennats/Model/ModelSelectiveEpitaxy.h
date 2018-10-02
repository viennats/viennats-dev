#ifndef MODELSELECTIVEEPITAXY_H_
#define MODELSELECTIVEEPITAXY_H_

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
#include <iostream>

using namespace std;

namespace model {
	///Model which applies a directional rate to surface depending on the normal direction of the surface

	class SelectiveEpitaxy {

		double CrystalDirection[3];
		double rate;


	public:
		static const bool OutputFluxes=false;
		static const bool SpatiallyEqualDistributedFlux=true;
		static const bool ReemissionIsMaterialDependent=false;
		bool CalculateConnectivities = false;
		bool CalculateVisibilities = false;
		bool CalculateNormalVectors = true;
		bool IncludeVectorRates = false;

		static const int CoverageStorageSize=0;
		static const int RatesStorageSize=0;
		static const unsigned int NumberOfParticleTypes=0;
		static const unsigned int NumberOfParticleClusters[1];

		class ParticleType {
		public:
			double Direction[3];
			double Flux;
		};



		SelectiveEpitaxy(const std::string & Parameters) {
			using namespace boost::spirit::classic;

			bool b = parse(
				Parameters.begin(),
				Parameters.end(),
				*(
					(str_p("crystal_direction")  >> '='  >> '{' >> real_p[assign_a(CrystalDirection[0])]  >> "," >> real_p[assign_a(CrystalDirection[1])] >> "," >> real_p[assign_a(CrystalDirection[2])] >> '}' >> ';') |
					(str_p("rate")  >> '='  >>  '{' >> (real_p[assign_a(rate)] % ',') >> '}'  >> ';')
				),
				space_p | comment_p("//") | comment_p("/*", "*/")).full;

				if (!b) msg::print_error("Failed interpreting model parameters!");

			}



			template<class VecType>
			void CalculateVelocity(
				double &Velocity,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible) const {

						Velocity = rate * 0;

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
						int Material
					) {}
				};

			}

			#endif /*MODELSELECTIVEEPITAXY_H_*/
