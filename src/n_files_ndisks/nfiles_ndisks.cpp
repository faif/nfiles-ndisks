/*
    nfiles_ndisks.cpp Try to fit a random number of files into a random 
                      number of hard disks
    Copyright (C) 2006 Athanasios Kasabalis

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <ga/GASimpleGA.h>	  // we're going to use the simple GA
#include <ga/GABin2DecGenome.h>   // our representation will be
                                  // binary to decimal
#include <ga/garandom.h>          // every GA needs random numbers...
#include <iostream>               // the standard C++ I/O stream
#include <fstream>                // the standard C++ file stream
#include <cmath>                  // the standard C math header
using namespace std;

#define OFILE "output.dat"        // the name of the output file
#define FFILE "fitness.dat"       // the name of the file with the fitness values

// the size of a disk (2-400 GB)
// it is global because we want to use it
// in the fitness function and it must
// remain unchanged for each gene, but 
// must be random for each execution :-)
unsigned int DISK_SIZE;
unsigned int N_DISKS;            // the number of disks
ofstream ofile;                 // the output file to write the results
unsigned int counter = 1;       // a global counter of individuals

// the max phenotype value of the individual
// one value for each disk
// it must be global so that it can be used
// in the fitness function, but at this point
// the number of the disks that will be used 
// is unknown, that's why we assign the maximum
// allowed number of disks
float max_values[15];                  

// the max fitness value of the best individual
float max_fitness;                     
float Objective(GAGenome &);           // the prototype of the fitness function
void total_size(GAGenome &, float[]);  // the prototype of the phenotype function

int
main()
{
  // create the output file
  ofile.open(OFILE);
  if(!ofile)
    {
      cerr << "Could not create output file, exiting..." 
	   << endl;
      return 1;
    }
  
  /* 
     By default GARandomSeed uses the current time multiplied by the 
     application's Process Identifier (on systems that support PIDS)
     as the seed for a pseudo-random number generator. On systems 
     which do not support PIDS, only the current time is used.
  */
  GARandomSeed(); 

  // -- initialize the GA parameters --

  // the maximum fitness of all executions
  // (and of all the valid individuals)
  max_fitness                   = 0.0f;

  // the size of a disk is 2-400 GB
  DISK_SIZE                     = GARandomInt(2, 400);

  // the number of files is 2-30 for each atom
  const unsigned int N_FILES    = GARandomInt(2, 30);

  // the number of disks is half of the number of files
  N_DISKS                       = GARandomInt(2, N_FILES/2);

  // the number of atoms (file teams) is 2-100
  const unsigned int POPULATION = GARandomInt(2, 100);

  // the possibility to use crossover (%)
  const float P_CROSS           = 0.7f;

  // the possibility to use mutation (%)
  const float P_MUT             = 0.1f;

  // the number of maximum genes
  //const unsigned int MAX_GENES  = 200;

  // threshhold for when we have converged (%)
  const float P_CONV            = 0.99f;
  // how many generations back to look
  const unsigned int N_CONV     = 50;
 
  // -- end of initialization --

  // this is not exactly the phenotype,
  // it just only shows in which disk
  // a file goes,
  // but we have to use it cause we
  // will use a real representation for
  // the genome and we want an integer
  // representation for the phenotype
  GABin2DecPhenotype phen;
  for (unsigned int i=0; i<N_FILES; i++)
    phen.add(N_FILES, 0, N_DISKS);
  
  // create the genome
  GABin2DecGenome genome(phen, Objective);

  cout << "executing the GA, please wait";

  // create the GA, set the parameters and execute it
  // this practically means call the fitness function for
  // each atom and evaluate it
  GASimpleGA ga(genome);
  ga.populationSize(POPULATION);
  //ga.nGenerations(MAX_GENES);
  ga.pMutation(P_MUT);
  ga.pCrossover(P_CROSS);
  // write the fitness of its individual in a file
  ga.scoreFilename(FFILE);
  // dump scores to disk every 20th generation
  ga.flushFrequency(20);
  // in case we want to use convergence instead of max genes
  ga.pConvergence(P_CONV);
  ga.nConvergence(N_CONV);
  ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);

  // execute the GA and initialize the genome (chromosome)
  ga.evolve();
  genome.initialize();

  // give some information about the execution
  ofile << "Size of each disk is " << DISK_SIZE << "Gb" << endl;
  ofile << "Number of disks is " << N_DISKS << endl;
  ofile << "Number of files for each individual: "
	<< N_FILES << endl;
  /*ofile << "Population: " << POPULATION
	<< " individuals" << endl;
  ofile << "Max generations: " << MAX_GENES
  << endl;*/

  // print the phenotype
  // this is in decimal format, we don't print it...
  //ofile << "The GA found:";
  //ofile << ga.statistics().bestIndividual() << " ";
  
  // print the best results 
  ofile << endl << "The amounts of filled disks are:" 
	<< endl;
  for (unsigned int i=0; i<N_DISKS; i++)
    {
      ofile << (i+1) << ": " 
	    << max_values[i] << "Gb " <<endl;
    }

  // close the output file
  ofile.close();

  // print some information to the user
  // after the execution
  cout << "the results have been written to " 
       << OFILE << endl;
  cout << "the fitness results have been written "
       << "to " << FFILE << endl;

  return 0;
}

// the fitness function (evaluates each individual/atom)
float
Objective(GAGenome& g) 
{
  float total[N_DISKS]; // the total sizes of files for each disk

  // save the total sizes of files for each disk
  // in the total matrix
  total_size(g, total);

  float fitness=0.0f;   // maximum fitness value: ?

  //get the fileteam's total size from the phenotype
  for (unsigned int i=0; i<N_DISKS; i++)
    {
      //ofile << "total[i] " << total[i] << endl;
      float remaining = abs((float)(DISK_SIZE - total[i])); //count in gigabytes
      // after this save the disk is full (best situation)
      if (total[i] == DISK_SIZE)
	{
	  fitness = 1;
	}
      // the individual(fileteam) can be saved in the hard disk
      else if (total[i] < DISK_SIZE)
	{
	  fitness = 1/(float)remaining;
	}
      // the individual cannot be saved in the hard disk
      // punish him by assigning a smaller fitness value
      else 
	{
	  fitness = (1/(float)remaining * 0.1);
	}
    }

  ofile << "fitness of individual is " 
	<< fitness << endl << endl;

  // although we found a better fitness, it is possible
  // that some of the total matrix values exceed the
  // DISK_SIZE limit, that's why we should check it
  if(fitness > max_fitness) 
    {
      bool better=true;
      for (unsigned int i=0; i<N_DISKS; i++)
	if (total[i]  > DISK_SIZE) // the individual is inappropriate
	  better=false;

      if (better) // change only if the individual is acceptable
	{
	  for (unsigned int i=0; i<N_DISKS; i++)
	    max_values[i]=total[i];
	  
	  max_fitness = fitness; // this is the best fitness until now
	  //ofile << "best fitness is " << max_fitness << endl;
	}
    }

  return fitness;
}


// creates the phenotype (total size of files of an atom)
void
total_size(GAGenome& g, float tot[])
{
  cout << ".";

  GABin2DecGenome & genome = (GABin2DecGenome &)g;
  unsigned int count[N_DISKS];   // the number of files which will go in each disk

  // fill the matrices with zeros to avoid garbage...
  for (unsigned int i=0; i<N_DISKS; i++)
    {
      count[i] = 0;
      tot[i] = 0;
    }

  // calculate the number of files for each disk
  for (int i=0; i<genome.nPhenotypes(); i++)
    {
      // keep just the integer value
      // the real representation is useless
      int tmp=round(genome.phenotype(i)); 

      count[tmp]++;
    }

  /*ofile << "count matrix" << endl;
  for (unsigned int i=0; i<N_DISKS; i++)
    ofile << count[i] << " ";
    ofile << endl;*/

  // add the file sizes to later check if they fit
  // in the disk (total[1] is the total size of files that 
  // will try to be saved in disk 1, total[2] in disk 1,
  // etc)
  for (unsigned int i=0; i<N_DISKS; i++)
  {
    for (unsigned int j=0; j<count[i]; j++)
      {
	//ofile << "count i " << count[i] << endl;
	float rnd = GARandomFloat(0.004, DISK_SIZE); // count in gigabytes
	tot[i] += rnd;
	/*ofile << "file " << (j+1) << " in disk "  
	  << (i+1) << " has size " << rnd << " Gb" << endl;*/
      }
  }

  ofile << "individual " << counter++ << endl;

  // write the total size of files for each disk in the output file 
  for(unsigned int i=0; i<N_DISKS; i++)
    {
      ofile << "total size of files for disk " << (i+1)
	    << " is " << tot[i] 
	    << "Gb" <<endl;
    }
}
