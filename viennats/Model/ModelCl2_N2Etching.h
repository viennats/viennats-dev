#ifndef MODELCl2_N2_ETCHING_H_
#define MODELCl2_N2_ETCHING_H_

#include <stack>
#include "../Statistics.h"
#include <cassert>
#include <string>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"


namespace model{
    template<class PPT>
    class Cl2_N2Etching {
        double StartDirection[3];
        double temperature;
        double pressure;
        double flux;

        double const stickingCoeff = 5e-7;     //just test different values
        double Flux_ev;            //evaporation flux of TiCl from surface
        double const rho_TiN = 5.0804e22;          //particles/cm^3
        double const rho_Cl = 4.2399e19;
        double const kB_over_m_Cl = 234.6938;

        std::vector<int> ActiveLayers;

    public:
        static const int CoverageStorageSize=1;
		static const int RatesStorageSize=2;

        static const unsigned int NumberOfParticleTypes=1;
        unsigned int NumberOfParticleClusters[1];

        static const bool OutputFluxes=false;
		static const bool SpatiallyEqualDistributedFlux=true;

		static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;

        static const bool CalculateNormalVectors=false;

        static const bool ReemissionIsMaterialDependent=true;

        class ParticleType {
        public:
            double Direction[3];
            double Probability;
            double Energy;
            double Flux;
            int Type;
        };

        Cl2_N2Etching(typename std::list<PPT>::const_iterator p){
            ActiveLayers = p->ActiveLayers;

            using namespace boost::spirit::classic;

            double Accuracy;

            bool b = parse(
                p->ModelParameters.begin(),
                p->ModelParameters.end(),
                *(
                    (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                    (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';') |
                    (str_p("flux") >> '=' >> real_p[assign_a(flux)] >> ';') |
                    (str_p("pressure") >> '=' >> real_p[assign_a(pressure)] >> ';') |
                    (str_p("temperature") >> '=' >> real_p[assign_a(temperature)] >> ';')
                ),
                space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

            NumberOfParticleClusters[0]=static_cast<unsigned int>(Accuracy);

            Flux_ev = flux*std::exp(-0.5/(8.617e-5*temperature));
            //std::cout << "Evaporation Flux: " << Flux_ev << std::endl;

        }

        template <class VecType>
		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
            if(Material<(int)ActiveLayers.size()){
                Velocity = -Rates[0]*Coverages[0]/(rho_TiN*100);        //m/s
            }
            //else if(Material==(int)ActiveLayers.size()) Velocity = -Rates[0]*Coverages[0]/(rho_TiN*100)/5;
            else{
                Velocity = 0;
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

        static void UpdateCoverage(double *Coverages, const double *Rates) {
            Coverages[0] = Rates[0]==0?0:Rates[0]/(Rates[0]+Rates[1]);
        }

        template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {
            p.Type = ParticleType;
            my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
            p.Energy=1.;
            p.Probability=1.;
            p.Flux=0.25*rho_Cl*std::sqrt(3*kB_over_m_Cl*temperature);        //flux just by thermal energy TODO: add pressure dependence rho_Cl
        }

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {
            Rates[0] += p.Flux*stickingCoeff;
            Rates[1] = Flux_ev;
        }

        template <class PT, class VecType> void ParticleReflexion(
        							const PT& p,
        							std::stack<PT>& particle_stack,
        							const VecType& NormalVector,
        							const double* Coverages,
        							int Material) const {
            double new_probability = p.Probability*(1.-stickingCoeff);
            if(new_probability >= (1.-stickingCoeff)){
                std::cout << p.Probability;
                particle_stack.push(p);
				PT& p_new = particle_stack.top();
				p_new.Probability=0.;
				my::stat::CosineNDistributedRandomDirection(1.,NormalVector,p_new.Direction);
            }
        }


    };
}

#endif
