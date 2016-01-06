//
//  main.cpp
//  mirna_aligner
//
//  Created by Alex P on 12/30/15.
//  Copyright (c) 2015 Alex Pankov. All rights reserved.
//


#include <iostream>
#include <stdio.h>
#include "experiment_reads.h"



static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   mirna_align shingle_size dl_distance_threshold shingle_match_percentage hairpin_file.fa.gz mature_file.fa.gz <first_mirna_library.fastq.gz[,...,nth_mirna_library.fastq.gz]> \n");
  fprintf(stderr, "\n");
  return 1;
}


int main(int argc, const char * argv[]) {
  clock_t t1 = std::clock();
  time_t tm = std::time(NULL);
  std::srand(static_cast<unsigned int>( tm ));
  std::cerr << std::asctime(std::localtime(&tm));
  
  unsigned int k = 8;
  int thresh = -1;
  float cutoff = 0.9;
  
  std::string hairpin_file, matures_file, srna_libraries;
  
  if ( argc == 7 )
  {
    k = (unsigned int)atoi(argv[1]) ;
    thresh = atoi(argv[2]) ;
    cutoff = atof(argv[3]);
    
    hairpin_file = argv[4];
    matures_file = argv[5];
    
    srna_libraries = argv[6];
  }
  
  
  std::vector<std::string> libraries_vec = split(srna_libraries, ',');
  
  
  experiment_reads my_experiment(libraries_vec);
  
  
  my_experiment.mat_to_map(matures_file);
  my_experiment.hair_to_map(hairpin_file);
  my_experiment.make_shingles(k);
  
  my_experiment.relate_mature_to_hairpin();

  std::cerr << "\nIndexing took: " << (float) (std::clock() - t1)/CLOCKS_PER_SEC << " seconds\n\n";
  
  my_experiment.run_libraries(cutoff, thresh);
  

  std::cerr << "\nRunning Time in secs: " << (float) (std::clock() - t1)/CLOCKS_PER_SEC << '\n';
  return 0;
}
