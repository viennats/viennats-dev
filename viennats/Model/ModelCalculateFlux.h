#ifndef MODELCALCULATEFLUX_H_
#define MODELCALCULATEFLUX_H_

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <stack>
#include <set>
#include "../Statistics.h"
#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_insert_key_actor.hpp>
#include <boost/spirit/include/classic_insert_at_actor.hpp>
#include "../parser_actors.h"
#include "../message.h"


namespace model {
///Model to calculate flux along an interface depending on the reflection used

  class CalculateFlux {

    double TotalFlux;
    double StickingProbability;
    double Yield;
    double end_probability;
    double StartAngleDistribution;
    double ReemittedAngleDistribution;

    double StartDirection[3];

    enum ReflectionModelType {NONE, DIFFUSIVE_SINGLE, DIFFUSIVE_MULTIPLE, SPECULAR, DIFFUSIVE_COMBINED} ;
    ReflectionModelType ReflectionModel;

    std::vector<int> reflection_materials_temp;
    std::set<int> reflection_materials;

    int reflection_diffusive_upper_bound;
//    int reference_start;


  public:
    static const bool OutputFluxes=true;
    static const bool SpatiallyEqualDistributedFlux=true;
    static const bool ReemissionIsMaterialDependent=true;
    static const bool CalculateConnectivities=false;
    static const bool CalculateVisibilities=false;
    static const bool CalculateNormalVectors=false;
    static const int CoverageStorageSize=0;
    static const int RatesStorageSize=1;
    static const unsigned int NumberOfParticleTypes=1;
    unsigned int NumberOfParticleClusters[NumberOfParticleTypes];

    class ParticleType {
    public:
      double Probability;
      double Direction[3];
      double Flux;
    };


    CalculateFlux(const std::string & Parameters):StartAngleDistribution(1.),ReemittedAngleDistribution(1.) {

        double Accuracy;
            using namespace boost::spirit::classic;

            // Default reflection model
            ReflectionModel = DIFFUSIVE_SINGLE;

            // Default number of generated particles for diffusive reflection
            reflection_diffusive_upper_bound = 4;

            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                            (str_p("direction")  >> '='  >> '{' >> real_p[assign_a(StartDirection[0])]  >> "," >> real_p[assign_a(StartDirection[1])] >> "," >> real_p[assign_a(StartDirection[2])] >> '}' >> ';') |
                            (str_p("reflection_model")  >> '='  >> (
                                str_p("DIFFUSIVE_SINGLE")   [parser_actors::assign_enum<ReflectionModelType>(ReflectionModel,DIFFUSIVE_SINGLE)] |
                                str_p("DIFFUSIVE_MULTIPLE") [parser_actors::assign_enum<ReflectionModelType>(ReflectionModel,DIFFUSIVE_MULTIPLE)] |
                                str_p("SPECULAR")           [parser_actors::assign_enum<ReflectionModelType>(ReflectionModel,SPECULAR)] |
                                str_p("DIFFUSIVE_COMBINED") [parser_actors::assign_enum<ReflectionModelType>(ReflectionModel,DIFFUSIVE_COMBINED)] |
                                str_p("NONE")               [parser_actors::assign_enum<ReflectionModelType>(ReflectionModel,NONE)]
                                )  >> ';') |
                            (str_p("reflection_materials")  >> '='  >>  '{' >> (int_p[push_back_a(reflection_materials_temp)] % ',') >> '}'  >> ';') |
                            (str_p("reflection_diffusive_upperbound")  >> '='  >> int_p[assign_a(reflection_diffusive_upper_bound)]  >> ';') |
//                            (str_p("reference_start")  >> '='  >> real_p[assign_a(reference_start)]  >> ';') |
                            (str_p("flux")  >> '='  >> real_p[assign_a(TotalFlux)]  >> ';') |
                            (str_p("sticking_probability")  >> '='  >> real_p[assign_a(StickingProbability)]  >> ';') |
                            (str_p("start_angle_distribution")  >> '='  >> real_p[assign_a(StartAngleDistribution)]  >> ';') |
                            (str_p("reemitted_angle_distribution")  >> '='  >> real_p[assign_a(ReemittedAngleDistribution)]  >> ';') |
                            (str_p("stop_criterion")  >> '='  >> real_p[assign_a(end_probability)]  >> ';') |
                            (str_p("yield")  >> '='  >> real_p[assign_a(Yield)]  >> ';') |
                            (str_p("statistical_accuracy")  >> '='  >>real_p[assign_a(Accuracy)]  >> ';')

                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            // transfer reflection materials from vector to set for faster lookup
            std::copy(reflection_materials_temp.begin(), reflection_materials_temp.end(), std::inserter(reflection_materials, reflection_materials.begin()));
            reflection_materials_temp.clear();

            if (!b) msg::print_error("Failed interpreting process parameters!");
/*

             end_probability=0.01;
             IonFlux=3.125e15;          // atoms/(cm²s)
             CFxFlux=2e18;              // atoms/(cm²s)
             yDPolyI=10e-24;            // cm³/atom
             yDPolyN=0.5e-24;           // cm³/atom
             StickingProbabilityCFxonPolymer=0.1;
             StickingProbabilityCFxonSi=0.1;
             StickingProbabilityCFxonMask=0.1;
             IonAngleDistribution=410.;
             CFxAngleDistribution=1.;*/

             NumberOfParticleClusters[0]=(TotalFlux>0.)?static_cast<unsigned int>(Accuracy):0;
        }

    template<class VecType>
        void CalculateVelocity(
            double &Velocity,
            const VecType& NormalVector,
            const double *Coverages,
            const double *Rates,
            int Material, bool Connected, bool Visible) const {

        Velocity=0.;//Rates[0]*Yield;
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

        template <class DropletType, class ParameterType, class PartitionType>
        void DropletGeneration(DropletType& d, double* Position, const ParameterType& Parameter, const PartitionType& Partition) const {}

//        template <class PT, class ParameterType, class PartitionType>
//        void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position, const ParameterType& Parameter, const PartitionType& Partition) const {
    template <class PT> void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) const {
      my::stat::CosineNDistributedRandomDirection(StartAngleDistribution,StartDirection,p.Direction);
      p.Probability=1.;
      p.Flux=TotalFlux;
    }

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const //, int Material) const
    {
         Rates[0]+=p.Flux*p.Probability;
    }

    template <class PT, class VecType> void ParticleReflexion(
                            const PT& p,
                            std::stack<PT>& particle_stack,
                            const VecType& NormalVector,
                            const double* Coverages,
                            int Material,
                            int D,
                            double dot // dot product between the incoming particle direction and the normal vector
                            ) const {

          if(ReflectionModel == NONE) return;

     // perform the reflection if either the material ID is member of the materials which is suppsoed to reflect
     // or if the material list is empty, which has the implicit meaning that all materials should reflect.
     //
     if((reflection_materials.find(Material) != reflection_materials.end()) || reflection_materials.empty())
     {
       // this is the old approach, generate one new particle and assign a random direction to it
       //
       if(ReflectionModel == DIFFUSIVE_SINGLE)
       {
          double new_probability=p.Probability*(1.-StickingProbability);
          if(new_probability>=end_probability) {
            particle_stack.push(p);
            PT& p_new=particle_stack.top();
            p_new.Probability=new_probability;
            my::stat::CosineNDistributedRandomDirection(ReemittedAngleDistribution,NormalVector,p_new.Direction);
          }
       }
       // diffusive reflection: generate several particles and assign random directions to them
       //
       else if(ReflectionModel == DIFFUSIVE_MULTIPLE)
       {
//          double new_probability=(p.Probability*(1.-StickingProbability))/(1.*reflection_diffusive_upper_bound);
          double new_probability=p.Probability*(1.-StickingProbability);
          if(new_probability>=end_probability)
          {
//            std::cout << particle_stack.size() << " vs. ";
            // evenly distribute the probability among the newly generated particles
            new_probability /= reflection_diffusive_upper_bound;
            for(int i = 0; i < reflection_diffusive_upper_bound; i++)
            {
              particle_stack.push(p);
              PT& p_new=particle_stack.top();
              p_new.Probability=new_probability;
              my::stat::CosineNDistributedRandomDirection(ReemittedAngleDistribution,NormalVector,p_new.Direction);
//              std::cout << p_new.Direction[0] << " "<< p_new.Direction[0] << " " << std::endl;
            }
//            exit(0);
//            std::cout << particle_stack.size() << std::endl;
          }
       }
       // specular reflection: generate one new particle and assign the reflected angle to it
       //
       else if(ReflectionModel == SPECULAR)
       {
          double new_probability=p.Probability*(1.-StickingProbability);
          if(new_probability>=end_probability) {
            particle_stack.push(p);
            PT& p_new=particle_stack.top();
            p_new.Probability=new_probability;

            // assign the new particle the reflection vector as direction: R = V - 2N<V,N>
            // http://mathworld.wolfram.com/Reflection.html
            for (int d=0;d<D;++d) p_new.Direction[d] = p.Direction[d] - 2.0 * NormalVector[d] * dot;
          }
       }
       else if(ReflectionModel == DIFFUSIVE_COMBINED)
       {
          double new_probability=p.Probability*(1.-StickingProbability);
          if(new_probability>=end_probability)
          {
            int reflection_diffusive_upper_bound_temp = 1;
            if(p.Probability >= 0.1)
              reflection_diffusive_upper_bound_temp = reflection_diffusive_upper_bound;

            new_probability /= reflection_diffusive_upper_bound_temp;

            for(int i = 0; i < reflection_diffusive_upper_bound_temp; i++)
            {
              particle_stack.push(p);
              PT& p_new=particle_stack.top();
              p_new.Probability=new_probability;
              my::stat::CosineNDistributedRandomDirection(ReemittedAngleDistribution,NormalVector,p_new.Direction);
            }
          }
       }
       else
       {
          std::cerr << "Reflection model is not supported! - Aborting .." << std::endl;
          exit(-1);
       }

     }
    }
  };


}




#endif /*MODELCALCULATEFLUX_H_*/
