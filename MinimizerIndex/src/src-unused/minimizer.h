#ifndef _MM_SKETCH_
#define _MM_SKETCH_

#include "loader.h"
#include "util_func.h"
#include "murmur3.hpp"
#include "xxhash64.h"
#include "ksort.h"

#include <boost/filesystem.hpp>
#include <algorithm>
#include <functional>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <tuple>
#include <string>
#include <list>
#include <time.h>
#include <sys/stat.h>
#include<set>
#include <cstring>
#include <sstream>
#include <locale>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <limits>
#include <omp.h>
#include <bitset>


#define MAX_BLOCK 128

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;


struct MMIndex{
    uint64_t info_;

    MMIndex(const uint64_t& info ){
      info_ = info;
    }
 

    //bool operator<(const MMIndex &rhs) const { return targetSeq < rhs.targetSeq; };
   

};
 
//NEW MMHASH impl -- same as minimap2's mm128_t struct
struct MMHash{  
    uint64_t x, y; 
    
   /* MMHash(const uint64_t& hashval_, const uint64_t& mm_info_ ){
        x= hashval_;
        y= mm_info_;
    }*/
      
};

//OLD MMHASH impl
/*struct MMHash{
    uint64_t hashValue;
    int position;
    int strand;

    MMHash(const int& hashvalue_, const int& position_, const int& strand_ ){
        hashValue = hashvalue_;
        position = position_;
        strand = strand_ ;
    }

};*/

struct MMMatch{
    uint32_t id;
    int strand;
    int start;
    int end;

    MMMatch(const uint32_t& id_, const int& strand_, const int& start_, const int& end_){
        id = id_;
        strand = strand_;
        start = start_;
        end = end_;
    }
};

/*struct comparer{
    inline bool operator()(const MMMatch& one, const MMMatch& two){

            return (one.targetSeq < two.targetSeq || one.targetSeq==two.targetSeq && one.strand < two.strand ||
                one.targetSeq==two.targetSeq && one.strand == two.strand && one.start < two.start ||
               one.targetSeq==two.targetSeq && one.strand == two.strand && one.start == two.start && one.end < two.end);

    }
};*/

struct hash_comparer{
    inline bool operator()(const MMHash& one, const MMHash& two){

            return (one.x < two.x || one.x==two.x && one.y < two.y);

    }
};

struct comparer{
    inline bool operator()(const MMIndex& one, const MMIndex& two){

            return (one.info_ <= two.info_);

    }
};

struct qcomparer{
    inline bool operator()(const MMMatch& one, const MMMatch& two){

            return (one.id < two.id || one.id==two.id && one.strand < two.strand ||
                one.id==two.id && one.strand == two.strand && one.start < two.start ||
               one.id==two.id && one.strand == two.strand && one.start == two.start && one.end < two.end);

    }
};

struct Cell {
    int score;
    int i;
    int j;
};


class MMSketch {
 public:
    MMSketch(std::vector<FastaInfo>& seqinfo, int& kmer, int& window);
    ~MMSketch(void);
    void PrintMMHashTable(void);
    void PrintIndexTable(void);
    void PrintMinimizers(void);
    std::string Map(FastaInfo& seq_entry, int& kmer, int& window, int epsilon);
    void checkQuery(std::vector<MMHash> & qhash, std::map<uint64_t, std::vector<MMIndex>>& mhash_table);
    
 private:
    void init(std::vector<FastaInfo>& seqinfo, int& kmer, int& window);
    void buildSketch();
    std::vector<std::pair<char, int>> x_drop_align(const std::string& s1, const std::string& s2);
    std::vector<std::pair<char, int>> get_cigar(const std::string& s1, const std::string& s2, const std::vector<Cell>& traceback);
    std::vector<int> longest_increasing_subset(std::vector<MMMatch>& matches, std::string orientation);
    std::vector<MMMatch> slice(std::vector<MMMatch>const & matches, int start, int end);
    std::map<uint64_t, std::vector<MMIndex>> loadIndex(std::map<uint64_t, std::vector<MMIndex>>& mhash_table);
    void computeMinimizerSketch(FastaInfo& sequence, int& kmer, int &window,std::vector<MMHash>& hashes);   
     std::vector<MMHash> computeMinimizerSketch(FastaInfo& seq_entry, int& k, int& w);
     void xDropAlignmentExtend(const std::string& seq1, const std::string& seq2,
                          int q_start, int q_end, int t_start, int t_end,
                          std::string& aligned_seq1, std::string& aligned_seq2);
     uint64_t hashSequence(const std::string& s);
     inline uint64_t invertableHash(uint64_t hash_val, uint64_t mask);
     std::string kmerSeq(const std::string& seq, int i, int k, int r);
     inline std::string reverseComplement(const std::string& seq);
     inline int hashValue(const char& n);
     std::vector<MMIndex> kmerTuples(std::string &sequence, int k);
     std::vector<std::string> kmerize(std::string & seq, int k);
    template<typename T>
    std::vector<T> vec2set(std::vector<T> kmers);
    uint64_t hash(const std::string & kmer);

 protected:
    int k, w;
    std::vector<FastaInfo> seqs;
    std::vector<MMHash> mm_hashes;
    std::unordered_map<uint64_t, std::vector<MMIndex>> MMHashTable;
    //std::vector<std::vector<MMIndex>> MMHashTable;
    std::unordered_map<std::string, std::string> TargetSeqMap;
    std::unordered_map<std::string, uint64_t> KHashMap;
    std::unordered_map<uint32_t, std::string> id_map;
    bool is_kmer_set_;
    bool is_mm_built_;
    bool is_window_set_;
 


};

#endif
