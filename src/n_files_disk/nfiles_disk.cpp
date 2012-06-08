/*
    nfiles_disk.cpp Try to fit a random number of files into one hard disk
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
#include <ga/GA1DBinStrGenome.h>  // our representation will be binary string
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
ofstream ofile;                 // the output file to write the results
unsigned int counter = 1;       // a global counter of individuals
float max_value = 0.0f;         // the max phenotype value of the individual (teamfile size)

float Objective(GAGenome &);    // the prototype of the fitness function
float total_size(GAGenome &);   // the prototype of the phenotype function

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

  // the size of a disk is 2-400 GB
  DISK_SIZE                     = GARandomInt(2, 400);
  
  // the number of files is 1-20 for each atom
  const unsigned int N_FILES    = GARandomInt(2, 30);

  // the number of atoms (file teams) is 2-100
  const unsigned int POPULATION = GARandomInt(2, 100);

  // the possibility to use crossover (%)
  const float P_CROSS           = 0.7f;

  // the possibility to use mutation (%)
  const float P_MUT             = 0.1f;

  // the number of maximum genes
  const unsigned int MAX_GENES    = 200;
 
  // -- end of initialization --

  // create a genome for the GA, it is necessary to be 
  // able to create a GA

  // GA1DBinaryStringGenome genome(N_FILES, Objective);
  GA1DBinaryStringGenome genome(N_FILES, Objective);

  cout << "executing the GA, please wait";

  // create the GA, set the parameters and execute it
  // this practically means call the fitness function for
  // each atom and evaluate it
  GASimpleGA ga(genome);
  ga.populationSize(POPULATION);
  ga.nGenerations(MAX_GENES);
  ga.pMutation(P_MUT);
  ga.pCrossover(P_CROSS);
  // write the fitness of its individual in a file
  ga.scoreFilename(FFILE);
  // dump scores to disk every 20th generation
  ga.flushFrequency(20);
  //ga.elitist(true);
  ga.evolve();

  // give some information about the execution
  ofile << "Number of files for each individual: "
	<< N_FILES << endl;
  /*ofile << "Population: " << POPULATION
	<< " individuals" << endl;
  ofile << "Max generations: " << MAX_GENES
  << endl;*/
  ofile << "Disk size: " << DISK_SIZE
	<< "Gb" << endl;
  // print out the best genome that the GA found
  ofile << "The GA found:" << ga.statistics().bestIndividual() << endl;
  ofile << "The size of the best individual is: " << max_value
	<< "Gb" << endl;

  // close the output file
  ofile.close();

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
  GA1DBinaryStringGenome & genome = (GA1DBinaryStringGenome &)g;
  float fitness=0.0f; // maximum fitness value: ?

  //get the fileteam's total size from the phenotype
  float total = total_size(genome);
  float remaining = abs((float)(DISK_SIZE - total)); //count in gigabytes

  // after this save the disk is full (best situation)
  if (total == DISK_SIZE)
    {
      fitness = 1;
      max_value = total;
    }
  // the individual(fileteam) can be saved in the hard disk
  else if (total < DISK_SIZE)
    {
      fitness = 1/(float)remaining;
      if (total > max_value) // we found a better result...
	max_value = total;
    }
  // the individual cannot be saved in the hard disk
  // punish him by assigning a smaller fitness value
  else 
    {
	  fitness=(1/(float)remaining * 0.1);
    }

  ofile << "fitness of individual is " 
	<< fitness << endl << endl;

  return fitness;
}


// create the phenotype (total size of files of an atom)
float
total_size(GAGenome& g)
{
  cout << ".";

  GA1DBinaryStringGenome & genome = (GA1DBinaryStringGenome &)g;

  float size[genome.length()];
  float total = 0.0f; // the total size of files

  ofile << "individual " << counter++ << endl;


  for(int i=0; i<genome.length(); i++)
    {
      // the size of a file is 4K up to DISK_SIZE
      // 1Kb is 0.001 Gb 
      if(genome.gene(i) == 1) // the file exists
	size[i] = GARandomFloat(0.004, DISK_SIZE); // count in gigabytes
      else size[i] = 0; // the file does not exist
      total +=size[i];
    }

  for(int i=0; i<genome.length(); i++)
    {
      ofile << "size of file " << i 
	    << " is " << size[i] 
	    << "Gb" <<endl;
    }

  ofile << "the total size of individual is: " << total 
	<< "Gb" << endl;

  return total; // the size is returned in gigabytes
}
