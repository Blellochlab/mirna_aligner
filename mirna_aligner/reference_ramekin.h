//
//  reference_ramekin.h
//  mirna_aligner
//
//  Created by Alex P on 1/1/16.
//  Copyright (c) 2016 Alex Pankov. All rights reserved.
//

#ifndef mirna_aligner_reference_ramekin_h
#define mirna_aligner_reference_ramekin_h

#include <vector>
#include <unordered_map>

struct mature;

struct hairpin{
  std::vector< std::string > hairpin_names;
  std::vector<unsigned long> start_pos;
  std::vector< std::unordered_map<std::string, mature >::iterator > mat_its;

  std::vector<unsigned long > counts_per_file;
  
  
  hairpin& operator+=(const std::string & myname){
    hairpin_names.push_back(myname);
    return *this;
  }

  void add_mat(std::unordered_map<std::string, mature >::iterator & it, unsigned long pos ){
    mat_its.push_back(it);
    start_pos.push_back(pos);
  }
  void initialize_size(size_t number_of_libraries){
    counts_per_file = std::vector<unsigned long>(number_of_libraries, 0);
  }

  void update_mat(std::unordered_map<std::string, mature >::iterator & mat_it, size_t start){
    mat_its.push_back(mat_it);
    start_pos.push_back(start);
  }
};


struct mature{
  std::vector<std::string> mature_names;
  std::vector< std::unordered_map<std::string, hairpin>::iterator > hair_its;
    
  std::vector<unsigned long > counts_per_file;
  
  
  void add(std::string & mature_name){
    mature_names.push_back(mature_name );
  }
  void initialize_size(size_t number_of_libraries){
    counts_per_file = std::vector<unsigned long>(number_of_libraries, 0);
  }

    
};


#endif
