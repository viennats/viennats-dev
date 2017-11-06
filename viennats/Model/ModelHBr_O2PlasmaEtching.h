#ifndef MODELHBr_O2PLASMAETCHING_H_
#define MODELHBr_O2PLASMAETCHING_H_

#include <stack>
#include "../Statistics.h"
#include <cassert>
#include <string>
#include <boost/spirit/include/classic.hpp>
#include "../message.h"

namespace model {

///HBr/O2 plasma etching model

class HBr_O2_PlasmaEtching {

	static const double kB;	// m2 kg s-2 K-1
	static const double rho_polySi; // 	//in atoms/cm^3 //2.6e22
	static const double rho_p;	//in atoms/cm^3

	static const double k_ei;
	static const double k_ev;

	static const double Gamma_ee;
	static const double Gamma_p;
	static const double Gamma_ep;

//		static const double A_p=0.0337;
	static const double Ae_ei;
	static const double Ae_sp;
	static const double Ap_ei;

	static const double B_sp;

	static const double Eth_e_e;
	static const double Eth_e_p;
	static const double Eth_e_sp;//

	double min_ion_energy;
	double delta_ion_energy;

	double temperature;

	static const double ion_exponent;

	static const double Phi_min;

	double StartDirection[3];

	double FluxIon;
	double FluxO;
	double FluxBr;

	double Flux_ev;



public:

	class ParticleType {
	public:
		double Direction[3];
		double Probability;
		double Energy;
		double Flux;
		int Type;
	};

	unsigned int NumberOfParticleClusters[3];

	static const int CoverageStorageSize=3;
	static const int RatesStorageSize=7;

	static const unsigned int NumberOfParticleTypes=3;

	static const bool OutputFluxes=false;
	static const bool SpatiallyEqualDistributedFlux=true;

	static const bool CalculateVisibilities=false;
	static const bool CalculateConnectivities=false;

	static const bool CalculateNormalVectors=true;

	static const bool ReemissionIsMaterialDependent=true;

//        int stacked_mat;


	HBr_O2_PlasmaEtching(const std::string & Parameters) {//, int stacked_material) {

		using namespace boost::spirit::classic;
//		    stacked_mat=stacked_material;

		double Accuracy;

		bool b = parse(
				Parameters.begin(),
				Parameters.end(),
				*(
						(str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
						(str_p("flux_ion")  >> '='  >> real_p[assign_a(FluxIon)]  >> ';') |
						(str_p("flux_oxygen")  >> '='  >> real_p[assign_a(FluxO)]  >> ';') |
						(str_p("flux_bromine")  >> '='  >> real_p[assign_a(FluxBr)]  >> ';') |
						(str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';') |
						(str_p("min_ion_energy")  >> '='  >>real_p[assign_a(min_ion_energy)]  >> ';') |
						(str_p("temperature")  >> '='  >>real_p[assign_a(temperature)]  >> ';') |
						(str_p("delta_ion_energy")  >> '='  >>real_p[assign_a(delta_ion_energy)]  >> ';')
				),
				space_p | comment_p("//") | comment_p("/*", "*/")).full;

		if (!b) msg::print_error("Failed interpreting process parameters!");

		unsigned int num_particles=static_cast<unsigned int>(Accuracy);

		NumberOfParticleClusters[0]=(FluxIon>0.)?num_particles:0;
		NumberOfParticleClusters[1]=(FluxO>0.)?num_particles:0;
		NumberOfParticleClusters[2]=(FluxBr>0.)?num_particles:0;

//		    temperature=300;
		Flux_ev=0.0027*FluxBr*std::exp(-0.168/(kB*temperature));
	}

	template <class VecType>
	void CalculateVelocity(double &Velocity, const VecType& NormalVector, const double *Coverages, const double *Rates, int Material, bool Connected, bool Visible) const {

		double ER_p = (Rates[2]*Coverages[0])/rho_p;
		double DR_p = (Rates[3])/rho_p;
		//std::cout << ER_p << " " << DR_p << std::endl;

		if (Coverages[1]>1) {	//Deposition
			if (Material==0) {	//Si
				Velocity = (DR_p - ER_p);
			} else {
				Velocity=0;
			}
		} else { 			//Etching
			if (Material<=1) {	//Si
				Velocity = -(Rates[1]*Coverages[2]+Rates[0]*(1.-Coverages[2])+Rates[6]*Coverages[2])/rho_polySi;
			} else {
				Velocity = 0;
			}
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

		//Active sites on polymer coverage e,P
		Coverages[0]=(Rates[5]==0)?0.:Rates[5]/(Rates[5]+Rates[2]);

		//Polymer coverage
		Coverages[1]=(Rates[3]==0.)?0.:((Coverages[0]==0.)||(Rates[2]==0.))?1.:Rates[3]/(Rates[2]*Coverages[0]);

		//Neutral coverage
		if (Coverages[1]<1.) {
			Coverages[2]=(Rates[4]==0)?0:Rates[4]*(1.-Coverages[1])/(Rates[4]+k_ei*Rates[1]+k_ev*Rates[6]);//*Flux_ev);
		} else {
			Coverages[2]=0.;//
		}
	}

	template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {

		p.Type=ParticleType;
		switch(ParticleType) {
			case 0: //ION
				my::stat::CosineNDistributedRandomDirection(ion_exponent,StartDirection,p.Direction);
				p.Probability=1;
				p.Energy=min_ion_energy+delta_ion_energy*my::stat::RandomNumber();
				p.Flux=FluxIon;
//					if (p.Flux==0) std::cout << "p.Flux = " << p.Flux << "\n";
				break;
			case 1:	//POLYMER				//OXYGEN -
				my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
				p.Energy=1.;
				p.Probability=1.;
				p.Flux=FluxO;
				break;
			case 2:	//NEUTRAL				//Bromine
				my::stat::Cosine1DistributedRandomDirection(StartDirection,p.Direction);
				p.Energy=1.;
				p.Probability=1.;
				p.Flux=FluxBr;	//5.5e18
				break;
			default:
				assert(0);
		}
	}

	template <class PT, class NormVecType> void ParticleCollision(
								const PT& p,
								const NormVecType& NormalVector,
								double* Rates,
								const double* Coverages,
								double RelTime) const {

		switch(p.Type) {
			case 0: {	//ION

				const double dot=NormalVector[0]*p.Direction[0]+NormalVector[1]*p.Direction[1]+NormalVector[2]*p.Direction[2];
				const double cos=std::min(std::max(-dot,0.),1.);

				const double sqrtE=std::sqrt(p.Energy);

				//------------NEW---------------------
				double f_e_ei=cos;
				double f_e_sp=std::max((1+B_sp*(1.-cos*cos))*cos,0.);

				const double Y_sp=Ae_sp *std::max(sqrtE-std::sqrt(Eth_e_sp),0.)*f_e_sp;
				const double Ye_ei=Ae_ei*std::max(sqrtE-std::sqrt(Eth_e_e),0.) *f_e_ei;
				const double Yp_ei=Ap_ei*std::max(sqrtE-std::sqrt(Eth_e_p),0.) *f_e_ei;

				Rates[0]+=p.Flux*Y_sp;
				Rates[1]+=p.Flux*Ye_ei;
				Rates[2]+=p.Flux*Yp_ei;

				Rates[6]=Flux_ev;//0.0027*p.Flux*std::exp(-0.168/(kB*temperature));

				break;
			}
			case 1: {		//POLYMER						OXYGEN
				const double S_p=std::max(Gamma_p*(1.-Coverages[2]-Coverages[0]*Coverages[1]),0.);
				Rates[3]+=p.Flux*S_p;

				break;

			}
			case 2: {		//NEUTRAL					Br
				const double S_ee=std::max(Gamma_ee*(1.-Coverages[2]),0.);//*(1-Coverages[0])),0.);//-Coverages[1]),0.);
				const double S_ep=std::max(Gamma_ep*((1.-Coverages[0])),0.);
				Rates[4]+=p.Flux*S_ee;
				Rates[5]+=p.Flux*S_ep;

				break;
			}
			default: {
				assert(0);
			}
		}
	}


	template <class PT, class VecType> void ParticleReflexion(
						const PT& p,
						std::stack<PT>& particle_stack,
						const VecType& NormalVector,
						const double* Coverages,
						int Material//,
//                            int D,
//                            double dot // dot product between the incoming particle direction and the normal vector
						) const {
	}

};

double const HBr_O2_PlasmaEtching::kB=0.000086173324;    // m2 kg s-2 K-1
double const HBr_O2_PlasmaEtching::rho_polySi=2.328e22;  //in atoms/cm^3 //2.6e22
double const HBr_O2_PlasmaEtching::rho_p=4.8e21;         //density of SiBr_xO_y in atoms/cm^3

double const HBr_O2_PlasmaEtching::k_ei=2;
double const HBr_O2_PlasmaEtching::k_ev=2;

double const HBr_O2_PlasmaEtching::Gamma_ee=0.04;		//Belen2006
double const HBr_O2_PlasmaEtching::Gamma_p=0.05;
double const HBr_O2_PlasmaEtching::Gamma_ep=0.05;

//		double const SiO2_PlasmaEtching::A_p=0.0337;
double const HBr_O2_PlasmaEtching::Ae_ei=0.015;//0.0361;
double const HBr_O2_PlasmaEtching::Ae_sp=0.005;//0.06;
double const HBr_O2_PlasmaEtching::Ap_ei=0.1;//0.1444;

double const HBr_O2_PlasmaEtching::B_sp=9.3;

double const HBr_O2_PlasmaEtching::Eth_e_e=4;
double const HBr_O2_PlasmaEtching::Eth_e_p=4;
double const HBr_O2_PlasmaEtching::Eth_e_sp=30;

double const HBr_O2_PlasmaEtching::ion_exponent=1000.;

double const HBr_O2_PlasmaEtching::Phi_min=1.3962634;

}

#endif /*MODELHBR_O2_PLASMAETCHING_H_*/
