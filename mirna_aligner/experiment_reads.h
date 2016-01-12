//
//  experiment_reads.h
//  mirna_aligner
//
//  Created by Alex P on 1/1/16.
//  Copyright (c) 2016 Alex Pankov. All rights reserved.
//

#ifndef mirna_aligner_experiment_reads_h
#define mirna_aligner_experiment_reads_h

#include <cmath>

#include <ctime>
#include <unistd.h>
#include <sys/types.h>


#include <zlib.h>
#include "kseq.h"
#include "reference_ramekin.h"

#include "some_functions.h"
#include "edit_dist.h"
#include "huff_codes.h"

KSEQ_INIT(gzFile, gzread)



struct experiment_reads{
  experiment_reads(std::vector< std::string > & library_names_vec){
    library_names = library_names_vec;
    number_of_libraries = library_names.size();
    library_index = 0;
  }
  size_t number_of_libraries, library_index;
  
  std::vector< std::string > library_names ;
  std::string mature_file_name, hairpin_file_name;
  
  std::unordered_map<std::string, mature > mat_map;
  std::unordered_map<std::string, hairpin > hair_map;
  
  std::unordered_map<std::string, std::vector< std::unordered_map<std::string, mature >::iterator > > mat_shingles;
  std::unordered_map<std::string, std::vector< std::unordered_map<std::string, hairpin >::iterator > > hair_shingles;

  unsigned int shingle_size;
  unsigned long new_reads;
  
  float cutoff;
  
  int threshold;

  struct readInfo{
    std::vector<unsigned long> counts;
    
    float hair_match_prob = 0., mat_match_prob = 0.;
    
    bool in_mature = false, in_hairpin = false, in_multi_mature = false, in_multi_hairpin = false, overhang = false, mismatch = false;
    std::unordered_map<std::string, hairpin>::iterator hair_it;
    std::unordered_map<std::string, mature>::iterator mat_it;
    
    void increment(experiment_reads & exp_reads){
      ++counts[exp_reads.library_index];
      if(in_mature) ++((mat_it->second).counts_per_file[exp_reads.library_index]);
      else if (in_hairpin) ++((hair_it->second).counts_per_file[exp_reads.library_index]);
    }
    
    void search_approximate(std::string & seq, std::unordered_set< std::string > & matched_hair_s, std::unordered_set< std::string > & matched_mat_s, experiment_reads & exp_reads){

      int min_mat_dist = INT_MAX, min_hair_dist = INT_MAX, min_dist = INT_MAX;
      std::vector<std::string> min_mat_v, min_hair_v;
      if(mat_match_prob > exp_reads.cutoff){
        
        for(auto i : matched_mat_s){
          int dist = dldist(seq, i);
          int excess_dist = std::abs( int(seq.length()) - int(i.length()) );
          int read_dist = dist - excess_dist;
          
          if(read_dist < min_mat_dist){
            min_mat_dist = read_dist;
            min_mat_v.clear();
            min_mat_v.push_back(i);
          } else if(read_dist == min_mat_dist){
            min_mat_v.push_back(i);
          }
        }
      }
      if(hair_match_prob > exp_reads.cutoff){
        for(auto i : matched_hair_s){
          int dist = dldist(seq, i);
          int excess_dist = std::abs( int(seq.length()) - int(i.length()) );
          int read_dist = dist - excess_dist;
          
          if(read_dist < min_hair_dist){
            min_hair_dist = read_dist;
            min_hair_v.clear();
            min_hair_v.push_back(i);
          } else if(read_dist == min_hair_dist){
            min_hair_v.push_back(i);
          }
        }
      }

      if(min_mat_dist < min_hair_dist){
        if(min_mat_dist <= exp_reads.threshold){
          min_dist = min_mat_dist;
          mismatch = true;
          in_mature = true;
          in_hairpin = true;
          if( min_mat_v.size() > 1){
            in_multi_mature = true;
          }
          mat_it = exp_reads.mat_map.find(min_mat_v[0]);
        }else if (min_hair_dist < min_mat_dist){
          if(min_hair_dist <= exp_reads.threshold){
            min_dist = min_hair_dist;
            mismatch = true;
            in_hairpin = true;
            if(min_hair_v.size() > 1){
              in_multi_hairpin = true;
            }
            hair_it = exp_reads.hair_map.find(min_hair_v[0]);
          }
        }else if (min_hair_dist == min_mat_dist && min_mat_dist < INT_MAX){
          if(min_mat_dist <= exp_reads.threshold){
            min_dist = min_mat_dist;
            mismatch = true;
            if(min_mat_v.size() <= min_hair_v.size()){
              in_mature = true;
              in_hairpin = true;
              if(min_mat_v.size() > 1) in_multi_mature = true;
              mat_it = exp_reads.mat_map.find(min_mat_v[0]);
            }else{
              in_hairpin = true;
              if(min_hair_v.size() > 1) in_multi_hairpin = true;
              hair_it =exp_reads.hair_map.find(min_hair_v[0]);
//  Havent checked if overlapping with a mature (thus relating it to a single mature read)
            }
          }
        }
      }
//      if( min_dist == 0 && in_mature && !in_multi_mature)
//        std::cout << seq << '\n' << mat_it->first << '\n';
//        std::cout << min_dist << '\n';
    }
    
    void initialize(std::string & seq, experiment_reads & exp_reads){
      counts = std::vector<unsigned long>(exp_reads.number_of_libraries, 0);
      unsigned int hair_match_int(0), mat_match_int(0);
      std::unordered_set< std::string > matched_hair_s, matched_mat_s;
      
      auto temp_set = shingle_string(seq, exp_reads.shingle_size);
      for(auto i = temp_set.begin(); i != temp_set.end(); ++i){
        auto search = exp_reads.hair_shingles.find( (*i) );
        auto mat_search = exp_reads.mat_shingles.find( (*i));
        if(search != exp_reads.hair_shingles.end()){
          ++hair_match_int;
          for(auto hair: (search->second)){
            matched_hair_s.insert(hair->first);
          }
        }
        if(mat_search != exp_reads.mat_shingles.end()){
          ++mat_match_int;
          for(auto mat: (mat_search->second)){
            matched_mat_s.insert(mat->first);
          }
        }
      }
      hair_match_prob = float(hair_match_int)/float(temp_set.size());
      mat_match_prob = float(mat_match_int)/float(temp_set.size());
      for(auto k: matched_mat_s){
        auto j = exp_reads.mat_map.find(k);
        auto seq_in_mat = (j->first).find( seq );
        if(seq_in_mat != std::string::npos ){
          if(!in_mature){
            in_mature = true;
            in_hairpin = true;
            mat_it = j;
          }
          else{
            in_multi_mature = true;
            in_multi_hairpin = true;
          }
        }
      }
      
      if(!in_hairpin){
        for(auto k: matched_hair_s){
          auto j = exp_reads.hair_map.find(k);
          auto seq_in_hair = (j->first).find( seq );
          if(seq_in_hair != std::string::npos ){
            
            if(!in_hairpin){
              in_hairpin = true;
              hair_it = j;
            }
            else{
              in_multi_hairpin = true;
            }
          }
        }
        if(in_hairpin && !in_multi_hairpin){
          unsigned times_in_mat = 0;
          auto seq_loc = (hair_it->first).find(seq);
          for(auto j = 0; j < (hair_it->second).start_pos.size(); ++j){
            auto this_start = (hair_it->second).start_pos[j];
            auto ol = overlap_amount(seq_loc, seq_loc + seq.length(), this_start, ((hair_it->second).mat_its[j]->first).length());
            if (ol > 0){
              in_mature = true;
              overhang = true;
              ++times_in_mat;
            }
          }
          if(times_in_mat > 1) in_mature = false;
        }
        if(!in_hairpin && !in_mature && exp_reads.threshold >= 0){
          search_approximate(seq, matched_hair_s, matched_mat_s, exp_reads);
        }
      }
    }
  };
  
  std::unordered_map<std::vector<bool>, readInfo> fq_counts;
  
  
  
  void mat_to_map(std::string mature_file){
    mature_file_name = mature_file;
    
    gzFile fp;
    kseq_t *seq;
    int l;
    
    fp = gzopen(mature_file_name.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      
      if ( std::string(seq->name.s).find("mmu") == std::string::npos) continue;
      std::string myseq(seq->seq.s);
      std::replace( myseq.begin(), myseq.end(), 'U', 'T');
      
      auto seq_name = std::string(seq->name.s) + " " + std::string(seq->comment.s);
      
      auto mat_map_it = mat_map.find(myseq);
      if ( mat_map_it == mat_map.end() ) {
        mat_map[myseq].initialize_size(number_of_libraries);
      }
      mat_map[myseq].add(seq_name);
    }
    
    kseq_destroy(seq);
    gzclose(fp);
  }
  
  
  
  void hair_to_map(std::string hairpin_file){
    hairpin_file_name = hairpin_file;
    
    gzFile fp;
    kseq_t *seq;
    int l;
    
    fp = gzopen(hairpin_file_name.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) {
      if ( std::string(seq->name.s).find("mmu") == std::string::npos) continue;
      std::string myseq(seq->seq.s);
      std::replace( myseq.begin(), myseq.end(), 'U', 'T');
      
      auto hair_map_it = hair_map.find(myseq);
      if ( hair_map_it == hair_map.end() ) {
        hair_map[myseq].initialize_size(number_of_libraries);
      }
      hair_map[myseq] += std::string(seq->name.s) + ' ' + std::string(seq->comment.s);
    }
    
    kseq_destroy(seq);
    gzclose(fp);
  }

  bool make_shingles(unsigned int k){
    if( k < 5) return false;
    if( k > 15) return false;
    shingle_size = k;
    
    std::unordered_map< std::string, std::unordered_set<std::string> > mat_shing_s, hair_shing_s;
    
    for(auto j = mat_map.begin(); j != mat_map.end(); ++j){
      auto temp_set = shingle_string(j->first, shingle_size);
      for (auto i: temp_set){
        (mat_shing_s[i]).insert( j->first );
      }
    }
    
    for(auto j = mat_shing_s.begin(); j != mat_shing_s.end(); ++j){
      for(auto k: j->second){
        mat_shingles[j->first].push_back( mat_map.find(k) );
      }
    }
    
    for(auto j = hair_map.begin(); j != hair_map.end(); ++j){
      auto temp_set = shingle_string(j->first, shingle_size);
      for (auto i: temp_set){
        (hair_shing_s[i]).insert( j->first );
      }
    }

    for(auto j = hair_shing_s.begin(); j != hair_shing_s.end(); ++j){
      for(auto k: j->second){
        hair_shingles[j->first].push_back( hair_map.find(k) );
      }
    }

    return true;
  }
  
  void run_seq(std::string & myseq){
    std::vector<bool> fq_bit;
    convert_seq_bool(myseq, fq_bit);
    auto fq_search = fq_counts.find(fq_bit);
    if (fq_search == fq_counts.end()) {
      fq_counts[fq_bit].initialize(myseq, *this);
      ++new_reads;
    }
    fq_counts[fq_bit].increment(*this);
  }
  
  void run_libraries(float exp_cutoff, int exp_thresh){
    cutoff = exp_cutoff;
    threshold = exp_thresh;
    for(auto library: library_names){
      std::cerr << "Mapping library: " << library << '\n';
      clock_t t1 = std::clock();
      new_reads = 0;
      
      gzFile fp;
      kseq_t *seq;
      int l;
      fp = gzopen(library.c_str(), "r");
      seq = kseq_init(fp);
      while ((l = kseq_read(seq)) >= 0) {
        auto my_seq = std::string(seq->seq.s);
        run_seq(my_seq);
      }
      kseq_destroy(seq);
      gzclose(fp);
      ++library_index;
      std::cerr << "Finished Mapping in " << (float) (std::clock() - t1)/CLOCKS_PER_SEC  << ". Observed " <<  new_reads << " from new reads.\n\n";
    }
  }
  
  void relate_mature_to_hairpin(){

    for(auto j = mat_map.begin(); j != mat_map.end(); ++j){
      auto temp_set = shingle_string(j->first, shingle_size);
      std::unordered_set< std::string > hairpin_set;
      
      for (auto i: temp_set){
        auto hair_it_vec = hair_shingles[i];
        for(auto k : hair_it_vec) hairpin_set.insert(k->first);
      }
      
      for(auto i: hairpin_set){
        auto mat_in_hair = i.find(j->first);
        if( mat_in_hair != std::string::npos){
          hair_map[i].add_mat(j, mat_in_hair);
        }
      }
    }
  }
  
};

#endif
