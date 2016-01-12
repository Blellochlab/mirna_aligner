//
//  some_functions.h
//  mirna_aligner
//
//  Created by Alex P on 1/1/16.
//  Copyright (c) 2016 Alex Pankov. All rights reserved.
//

#ifndef mirna_aligner_some_functions_h
#define mirna_aligner_some_functions_h

#include <vector>
#include <sstream>
#include <unordered_set>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}



std::unordered_set<std::string> shingle_string(std::string my_seq, int shing_len){
  std::unordered_set<std::string> shingle_set;
  for (unsigned int i = 0; i < my_seq.length() - shing_len + 1; ++i) {
    shingle_set.insert(my_seq.substr(i, shing_len));
  }
  return shingle_set;
}


//Need to check that this is correct for 2 right-open, left-closed intervals
int overlap_amount(int r1_min, int r1_max, int r2_min, int r2_max){
  if(r1_min >= r2_max || r2_min >= r1_max ) return 0;
  else if (r1_min < r2_max) return r1_max - r2_min;
  else return r2_max - r1_min;
}



#endif
