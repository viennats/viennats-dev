#ifndef MODELN2_FLASH_H_
#define MODELN2_FLASH_H_


namespace model{

    class N2_FLASH {
        double StartDirection[3];
        double temperature;
        double pressure;
        double flux;

        double const stickingCoeff = 1e-11;     //just test different values
        double const Flux_ev = 1e16;            //evaporation flux of SiN from surface
        double const rho_CF = 1.27e22;          //particles/cm^3
        double const rho_N2 = 5.38e19;

    public:
        static const int CoverageStorageSize=1;
		static const int RatesStorageSize=1;

        static const unsigned int NumberOfParticleTypes=3;

        static const bool OutputFluxes=false;
		static const bool SpatiallyEqualDistributedFlux=true;

		static const bool CalculateVisibilities=false;
        static const bool CalculateConnectivities=false;

        static const bool CalculateNormalVectors=true;

        static const bool ReemissionIsMaterialDependent=false;

        class ParticleType {
        public:
            double Direction[3];
            double Probability;
            double Energy;
            double Flux;
            int Type;
        };

        N2_FLASH(const std::string & Parameters){
            using namespace boost::spirit::classic;

            double Accuracy;

            bool b = parse(
                Parameters.begin(),
                Parameters.end(),
                *(
                    (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                    (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';') |
                    (str_p("flux") >> '=' >> real_p[assign_a(flux)] >> ';') |
                    (str_p("pressure") >> '=' >> real_p[assign_a(pressure)] >> ';') |
                    (str_p("temperature") >> real_p[assign_a(temperature)] >> ';')
                ),
                space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

            unsigned int num_particles=static_cast<unsigned int>(Accuracy);

        }

        template <class VecType>
		void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {
            if(Material==0){
                Velocity = 1e7*Rates[0]*Coverages[0]/rho_CF;        //nm/s
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
            Coverages[0] = Rates[0]/(Rates[0]+Flux_ev);
        }

        template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {
            p.Type = ParticleType;
            my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
            p.Energy=1.;
            p.Probability=1.;
            p.Flux=rho_N2*std::sqrt(2*kB*T/m_N2);        //flux just by thermal energy
        }

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {
            Rates[0] += p.Flux*stickingCoeff;
        }

        template <class PT, class VecType> void ParticleReflexion(
        							const PT& p,
        							std::stack<PT>& particle_stack,
        							const VecType& NormalVector,
        							const double* Coverages,
        							int Material) const {}


    }
}
