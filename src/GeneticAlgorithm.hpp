//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

namespace galgo {

//=================================================================================================

template <typename T>
class GeneticAlgorithm
{
   static_assert(std::is_same<float,T>::value || std::is_same<double,T>::value, "variable type can only be float or double, please amend.");

   template <typename K>
   friend class Population;
   template <typename K>
   friend class Chromosome;

   template <typename K>
   using Func = std::vector<K> (*)(const std::vector<K>&);

private:
   Population<T> pop;             // population of chromosomes
   std::vector<T> lowerBound;     // parameter(s) lower bound
   std::vector<T> upperBound;     // parameter(s) upper bound
   std::vector<T> initialSet;     // initial set of parameter(s)
   std::vector<int> idx;          // indexes for chromosome breakdown

public: 
   // objective function pointer
   Func<T> Objective; 
   // selection method initialized to roulette wheel selection                                   
   void (*Selection)(Population<T>&) = RWS;  
   // cross-over method initialized to 1-point cross-over                                
   void (*CrossOver)(const Population<T>&, CHR<T>&, CHR<T>&) = P1XO;
   // mutation method initialized to single-point mutation 
   void (*Mutation)(CHR<T>&) = SPM;  
   // adaptation to constraint(s) method                                        
   void (*Adaptation)(Population<T>&) = nullptr; 
   // constraint(s)                               
   std::vector<T> (*Constraint)(const std::vector<T>&) = nullptr; 

   T covrate = .50;   // cross-over rate
   T mutrate = .05;   // mutation rate   
   T SP = 1.5;        // selective pressure for RSP selection method 
   T tolerance = 0.0; // terminal condition (inactive if equal to zero)
                 
   int elitpop = 1;   // elit population size
   int matsize;       // mating pool size, set to popsize by default
   int tntsize = 10;  // tournament size
   int genstep = 10;  // generation step for outputting results
   int precision = 5; // precision for outputting results

   // constructor
   GeneticAlgorithm(Func<T> objective, int popsize, int nbgen, bool output, int nbparam);
   // run genetic algorithm
   void run();
   // return best chromosome 
   const CHR<T>& result() const;

private:
   int nbgen;     // number of generations
   int nogen = 0; // numero of generation
   int nbparam;   // number of parameters to be estimated
   int popsize;   // population size
   bool output;   // control if results must be outputted

   // print results for each new generation
   void print() const;
};

/*-------------------------------------------------------------------------------------------------*/
   
// constructor
template <typename T>
GeneticAlgorithm<T>::GeneticAlgorithm(Func<T> objective, int popsize, int nbgen, bool output, int nbparam)
{
   this->Objective = objective;
   // getting total number of bits per chromosome
   this->nbgen = nbgen;
   // getting number of parameters in the pack
   this->nbparam = nbparam;
   this->popsize = popsize;
   this->matsize = popsize;
   this->output = output;
}

/*-------------------------------------------------------------------------------------------------*/
   
// run genetic algorithm
template <typename T>
void GeneticAlgorithm<T>::run()
{
   // setting adaptation method to default if needed
   if (Constraint != nullptr && Adaptation == nullptr) {
      Adaptation = DAC;
   }

   // initializing population
   pop = Population<T>(*this);

   if (output) {
      std::cout << "\n Running Genetic Algorithm...\n";
      std::cout << " ----------------------------\n";
   }

   // creating population
   pop.creation();
   // initializing best result and previous best result
   T bestResult = pop(0)->getTotal();
   T prevBestResult = bestResult;
   // outputting results 
   if (output) print();
    
   // starting population evolution
   for (nogen = 1; nogen <= nbgen; ++nogen) {
      // evolving population
      pop.evolution();
      // getting best current result
      bestResult = pop(0)->getTotal();
      // outputting results
      if (output) print();
      // checking convergence
      if (tolerance != 0.0) {
         if (fabs(bestResult - prevBestResult) < fabs(tolerance)) {
            break;
         }
         prevBestResult = bestResult;
      }
   } 

   // outputting contraint value
   if (Constraint != nullptr) {
      // getting best parameter(s) constraint value(s)
      std::vector<T> cst = pop(0)->getConstraint(); 
      if (output) {
         std::cout << "\n Constraint(s)\n";
         std::cout << " -------------\n";
         for (unsigned i = 0; i < cst.size(); ++i) {
            std::cout << " C"; 
            if (nbparam > 1) {
               std::cout << std::to_string(i + 1);
            }
            std::cout << "(x) = " << std::setw(6) << std::fixed << std::setprecision(precision) << cst[i] << "\n"; 
         }
         std::cout << "\n"; 
      }
   }   
}

/*-------------------------------------------------------------------------------------------------*/

// return best chromosome
template <typename T>
inline const CHR<T>& GeneticAlgorithm<T>::result() const
{
   return pop(0);
}

/*-------------------------------------------------------------------------------------------------*/
   
// print results for each new generation
template <typename T>
void GeneticAlgorithm<T>::print() const
{
   // getting best parameter(s) from best chromosome
   std::vector<T> bestParam = pop(0)->getParam();
   std::vector<T> bestResult = pop(0)->getResult();

   if (nogen % genstep == 0) {
      std::cout << " Generation = " << std::setw(std::to_string(nbgen).size()) << nogen << " |";
      for (int i = 0; i < nbparam; ++i) {
	      std::cout << " X";
         if (nbparam > 1) {
            std::cout << std::to_string(i + 1);
         }
         std::cout << " = "  << std::setw(9) << std::fixed << std::setprecision(precision) << bestParam[i] << " |";
	   }
      for (unsigned i = 0; i < bestResult.size(); ++i) {
	      std::cout << " F";
         if (bestResult.size() > 1) {
            std::cout << std::to_string(i + 1);
         }
         std::cout << "(x) = " << std::setw(12) << std::fixed << std::setprecision(precision) << bestResult[i];
         if (i < bestResult.size() - 1) {
            std::cout << " |";
         } else {
            std::cout << "\n";
         }
	   }
 
   }
}
   
//=================================================================================================

}

#endif
