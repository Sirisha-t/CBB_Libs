#include "minimizer.h"
#include "kvec.h"


using namespace std;

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif


// Scoring parameters
const int match_score = 1;
const int mismatch_score = -1;
const int gap_score = -1;

// X-drop parameter
const int X_DROP = 50; // Adjust this value based on your requirements

const int GAP_PENALTY = 1; // gap penalty
const int X_DROP_THRESHOLD = 50; // X-drop threshold

//Protein overlap parameters
//const int MAPLEN = 7;
//const int MINHANG = 3;

//DNA overlap parameters
const int MAPLEN = 30;
const int MINHANG = 15;


unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static FORCE_INLINE uint64_t rotl64 ( uint64_t x, int8_t r )
{
  return (x << r) | (x >> (64 - r));
}

#define ROTL64(x,y)	rotl64(x,y)

#define BIG_CONSTANT(x) (x##LLU)
//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

#define getblock(p, i) (p[i])

static FORCE_INLINE uint64_t fmix64 ( uint64_t k )
{
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}


MMSketch::MMSketch(std::vector<FastaInfo>& seqinfo, int& kmer, int& window)  {
 init(seqinfo, kmer, window);
 //std::cout<<"MMSketch initialized successfully. Building sketch now for k="<<kmer<<",w="<<window<<"\n";
 MMSketch::buildSketch();
 //MMSketch::PrintIndexTable();
 return;
 
}

void MMSketch::init( std::vector<FastaInfo>& seqinfo, int& kmer, int& window){
  seqs = seqinfo;
  k = kmer;
  w = window;
  is_mm_built_ = false;
  is_kmer_set_ = true;
  is_window_set_ = true;

  for(auto seqs : seqinfo){
    id_map[seqs.rid_] = seqs.header_;
    TargetSeqMap[seqs.header_] = seqs.sequence_;
  }
}

MMSketch::~MMSketch(void) {
  return;
}

void MMSketch::PrintMMHashTable(void)  {
  /*(for(int i = 0; i < num_seqs_; ++ i) {
    cout << i << "  " << seq_entry[i].header_<<"  "<<seq_entry[i].sequence_<<"  "<<seq_entry[i].seq_len_<< endl;
  }*/
  return;
}

/** start of murmur hash implementation **/

inline string str_canonical(const std::string& kmer) {
    auto kmer_rev = kmer;
    std::reverse(kmer_rev.begin(), kmer_rev.end());
    for (size_t j = 0; j < kmer_rev.length(); ++j) {
        if (kmer_rev[j] == 'A') kmer_rev[j] = 'T';
        else if (kmer_rev[j] == 'T') kmer_rev[j] = 'A';
        else if (kmer_rev[j] == 'C') kmer_rev[j] = 'G';
        else if (kmer_rev[j] == 'G') kmer_rev[j] = 'C';
    }
    return kmer < kmer_rev ? kmer : kmer_rev;
}

uint64_t MurmurHash3_x64_128 ( const void * key, int len,
                           unsigned int seed )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 4;
  int i;

  uint64_t h1 = seed;

  uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f); 


  //----------
  // body

  const uint64_t * blocks = (const uint64_t *)(data);

  for(i = 0; i < nblocks; i++)
  {
    uint64_t k1 = getblock(blocks,i);

    k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;

    h1 = ROTL64(h1,27); h1 = h1*5+0x52dce729;

  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  uint64_t k1 = 0;

  switch(len & 8)
  {

  case  8: k1 ^= (uint64_t)(tail[ 7]) << 56;
  case  7: k1 ^= (uint64_t)(tail[ 6]) << 48;
  case  6: k1 ^= (uint64_t)(tail[ 5]) << 40;
  case  5: k1 ^= (uint64_t)(tail[ 4]) << 32;
  case  4: k1 ^= (uint64_t)(tail[ 3]) << 24;
  case  3: k1 ^= (uint64_t)(tail[ 2]) << 16;
  case  2: k1 ^= (uint64_t)(tail[ 1]) << 8;
  case  1: k1 ^= (uint64_t)(tail[ 0]) << 0;
           k1 *= c1; k1  = ROTL64(k1,31); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization
    h1 ^= len;
    h1 = fmix64(h1);

    return h1;
}

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed )
{
	const uint64_t m = 0xc6a4a7935bd1e995;
	const int r = 47;

	uint64_t h = seed ^ (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end)
	{
		uint64_t k = *data++;

		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h ^= k;
		h *= m; 
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7)
	{
	case 7: h ^= uint64_t(data2[6]) << 48;
	case 6: h ^= uint64_t(data2[5]) << 40;
	case 5: h ^= uint64_t(data2[4]) << 32;
	case 4: h ^= uint64_t(data2[3]) << 24;
	case 3: h ^= uint64_t(data2[2]) << 16;
	case 2: h ^= uint64_t(data2[1]) << 8;
	case 1: h ^= uint64_t(data2[0]);
	        h *= m;
	};
 
	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
} 

uint64_t MMSketch::hash(const std::string & kmer) {
    string canonical_kmer = str_canonical(kmer);
    const char *c = canonical_kmer.c_str();
    return MurmurHash3_x64_128(c, canonical_kmer.size(),42);
}

/** end of murmurhash implementation **/

inline uint64_t MMSketch::invertableHash(uint64_t key, uint64_t mask){
  key = (~key) + (key << 21) & mask; // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8) & mask; // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31) & mask;
  return key;
}


 inline int MMSketch::hashValue(const char& n){
    switch(n)
    {
        case 'A':
            return '0'; break;
        case 'C':
            return '1' ; break;
        case 'G':
            return '2'; break;
        case 'T':
            return '3'; break;
        default:
            break;

    }
    assert(true) ;return (' ');

}

 inline std::string MMSketch::reverseComplement(const std::string& seq){
    std::string rev_comp = seq;
    int seqlen = seq.length();
    std::reverse(rev_comp.begin(), rev_comp.end());
    for (std::size_t i = 0; i < seqlen; ++i){
        switch (seq[i]){
        case 'A':
            rev_comp[i] = 'T';
            break;    
        case 'C':
            rev_comp[i] = 'G';
            break;
        case 'G':
            rev_comp[i] = 'C';
            break;
        case 'T':
            rev_comp[i] = 'A';
            break;
        }
    }
    return rev_comp;
}


/**
 *
 * @param std::string seq: sequence
 * @param int i : position
 * @param int k : k-mer size
 * @param int r : reverse( if 1 return reverse complement, 0 for forward)
 * @return std::string: k-mer of size k at position i
 */

 std::string MMSketch::kmerSeq(const std::string& seq, int i, int k, int r){
    std::string kmer = seq.substr(i,k);

    if (r){
        return reverseComplement(kmer);
    }
    else
    {   
        return kmer;
    }
}

 uint64_t MMSketch:: hashSequence(const std::string& s ){
    uint64_t hash_value = 0;
    uint64_t final_hash = 0;
    int slen= s.length();

    /*for (int i = 0; i < slen; i++)
    {  
        uint64_t pow_val = pow(4, slen - (i + 1));
        hash_value += ((hashValue(s[i])) * pow(4, slen - (i + 1)));
    }*/

   /* ST: bitwise implementation of hash function */
    for(int i=0; i<slen; i++){
        size_t pow = slen-i-1;
        size_t val = 4;
        int bit_hash =0;
        int hash_char = hashValue(s[i]);
        uint64_t pow_res = 1;

        while(pow > 0){
            if( pow & 1){
                pow_res = pow_res * val;
            }
            val = val * val;
            pow = pow >> 1;
        }
        while(pow_res >0){
            if(pow_res & 1){
                bit_hash += hash_char;
            }
            hash_char = hash_char << 1;
            pow_res = pow_res >> 1;
        }
        final_hash += bit_hash; 
    }
    
       
    //uint64_t inverted_hash = invertableHash(final_hash);
    uint64_t inverted_hash=0;
    //std::cout<<s<<"\t"<<final_hash<<"\t"<<inverted_hash<<"\n";

    return inverted_hash;


}

template<typename T>
std::vector<T> MMSketch::vec2set(std::vector<T> kmers)
{
    std::set<MMIndex> s = std::set<MMIndex>(kmers.begin(), kmers.end());
    std::vector<MMIndex> ret = vector<MMIndex>(s.begin(), s.end());
    return ret;
}

std::vector<std::string> MMSketch:: kmerize(std::string & seq, int k){
    std::vector<std::string> kmer(seq.length()-k , "");
    for(int i=0; i < seq.length()-k; i++){
        std::string s = seq.substr(i, k);
        kmer[i] = s;
    }
    
    return kmer;
}

// Old function -- unused in current impl
std::vector<MMIndex> MMSketch:: kmerTuples(std::string &sequence, int k){
    std::vector<std::string> kmers = kmerize(sequence, k);
    std::vector<MMIndex> tuples;
    for(int i=0; i< kmers.size(); i++){
        //tuples.emplace_back(MMIndex(kmers[i],i,k));
        //MMIndex mm(kmers[i], i, k);
        //tuples = mm;
    }
    return tuples;

}

//// Modified sketch functions using Minimap's operations -- old implementation -- unused
/**std::vector<MMHash> MMSketch::computeMinimizerSketch(FastaInfo& seq_entry, int& k, int& w){
    // store hashes for each sequence into a vector
    std::vector<MMHash> hashes {};

    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1;
    uint64_t kmer[2] = {0,0}; //to store the forward and reverse kmers
	int i, j, l, buf_pos, min_pos, kmer_span = 0; // i is used to store current position in input seq, l stores length of current kmer
	MMHash buf[256], min = { UINT64_MAX, UINT64_MAX };

    assert(seq_entry.seq_len_ > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); 
	memset(buf, 0xff, w * 16);

    for (i = l = buf_pos = min_pos = 0; i < seq_entry.seq_len_; ++i) {
        
		int c = seq_nt4_table[(uint8_t)seq_entry.sequence_[i]];
		MMHash mm_info = { UINT64_MAX, UINT64_MAX };
       // std::cout<<i<<"  "<<l<<"\n";
        if (c < 4) { // not an ambiguous base
			int z;
			kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer -- here mask is used for clearing out any bits outside the given bit range
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
     
            if (l >= k && kmer_span < 256) {
				mm_info.x = invertableHash(kmer[z], mask) << 8 | kmer_span; //use upto 56 bits for hash value
				mm_info.y = (uint64_t)seq_entry.rid_<<32 | (uint32_t)i<<1 | z; //use only right most 32 bits of rid
			}
        }   
        buf[buf_pos] = mm_info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) 
                    {   
                        hashes.push_back(buf[j]);
                    }
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) 
                    hashes.push_back(buf[j]);
		}
        if (mm_info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) 
            {
                hashes.push_back(min);
            }
			min = mm_info, min_pos = buf_pos;
		}
        else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) 
            {
                hashes.push_back(min);
			}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y)
                        hashes.push_back(buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) 
                        hashes.push_back(buf[j]);
			}
		}
        if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
        hashes.push_back(min);


    return hashes;

}**/


//// Modified sketch functions using Minimap's operations 
void MMSketch::computeMinimizerSketch(FastaInfo& seq_entry, int& k, int& w, std::vector<MMHash> & hashes ){
    // store hashes for each sequence into a vector
    std::vector<MMHash> tmp_hashes;

    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1;
    uint64_t kmer[2] = {0,0}; //to store the forward and reverse kmers
	int i, j, l, buf_pos, min_pos, kmer_span = 0; // i is used to store current position in input seq, l stores length of current kmer
	MMHash buf[256], min = { UINT64_MAX, UINT64_MAX };

    assert(seq_entry.seq_len_ > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); 
	memset(buf, 0xff, w * 16);

//if(seq_entry.header_ == "Read_9311" || seq_entry.header_ == "Read_1591")
//{   
    //std::string sk = "";
    //std::cout<<"Target/Query : "<<seq_entry.header_<<"\n";
    for (i = l = buf_pos = min_pos = 0; i < seq_entry.seq_len_; ++i) {
		
        int c = seq_nt4_table[(uint8_t)seq_entry.sequence_[i]];
        //std::cout<<seq_entry.sequence_[i]<<"\t"<<c<<"\n";
		MMHash mm_info = { UINT64_MAX, UINT64_MAX };
       // std::cout<<i<<"  "<<l<<"\n";
        if (c < 4) { // not an ambiguous base
			int z;
			kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer -- here mask is used for clearing out any bits outside the given bit range
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
            //std::cout<<kmer[0]<<"\t"<<kmer[1]<<"\n";
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;

            //sk += seq_entry.sequence_[i];
            //std::cout<<"kmer seq: "<<sk<<"\n";

            if (l >= k && kmer_span < 256) {
                //sk = "";
				mm_info.x = invertableHash(kmer[z], mask) << 8 | kmer_span; //use upto 56 bits for hash value
				mm_info.y = (uint64_t)seq_entry.rid_<<32 | (uint32_t)i<<1 | z; //use only right most 32 bits of rid
			}
        }   
        //std::cout<<"MM info x andf y :"<<mm_info.x<<"\t"<<mm_info.y<<"\n";
        buf[buf_pos] = mm_info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
        if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) 
                    {   
                        hashes.push_back(buf[j]);
                        tmp_hashes.push_back(buf[j]);
                        //std::cout<<"First window 1st loop "<<buf[j].x<<"\t"<<buf[j].y<<"\n";
                    }
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y){
                        hashes.push_back(buf[j]);
                        tmp_hashes.push_back(buf[j]);
                        //std::cout<<"First window 2nd loop "<<buf[j].x<<"\t"<<buf[j].y<<"\n";
                }
                    
		}
        if (mm_info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) 
            {
                hashes.push_back(min);
                tmp_hashes.push_back(min);
            }
			min = mm_info, min_pos = buf_pos;
		}
        else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) 
            {
                hashes.push_back(min);
                tmp_hashes.push_back(min);
			}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y)
                        {
                            hashes.push_back(buf[j]);
                            tmp_hashes.push_back(buf[j]);
                        }
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) 
                        {
                            hashes.push_back(buf[j]);
                            tmp_hashes.push_back(buf[j]);
                        }
			}
		}
        if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX){
        hashes.push_back(min);
        tmp_hashes.push_back(min);
    }
        
    
  /* if(seq_entry.header_ == "Read_1816"){\
          std::cout<<"Read_1816 minimizers\n";
            
            for(auto &h : tmp_hashes){
                std::cout<<h.x<<"\t"<<h.y<<"\n";
           }
        }
    if(seq_entry.header_ == "Read_3181"){\
            std::cout<<"Read_3181 minimizers\n";
            
            for(auto &h : tmp_hashes){
                std::cout<<h.x<<"\t"<<h.y<<"\n";
           }
        }*/
        
    
//}
    
    return ;

}


//Old implementation-- needs to be removed

/*std::vector<MMHash> MMSketch::computeMinimizerSketch( std::string& sequence, int& k, int& w){
    // store hashes for each sequence into a vector
    std::vector<MMHash> hashes;
    int seqlen = sequence.length()-w-k+1;
    for(int i = 1; i <= seqlen; i++ ){
        uint64_t m = std::numeric_limits<uint64_t>::max();

        for( int j = 0; j < w; j++) {
            uint64_t u, v;
            std::string fwd, rev;
            //u = hashSequence(kmerSeq(sequence, i + j, k, 0));
            //v = hashSequence(kmerSeq(sequence, i + j, k, 1));

            fwd = kmerSeq(sequence, i + j, k, 0);
            rev = kmerSeq(sequence, i + j, k, 1);
            //const char* fwd_kseq = forward.c_str();
            //const char* rev_kseq = reverse.c_str();

           // uint64_t u = hash(forward);
            //uint64_t v = hash(reverse);


            //MurmurHash3_x64_128(fwd_kseq, strlen, 42, fhash);
            //MurmurHash3_x64_128(rev_kseq, strlen, 42, rhash);
            //uint64_t* u = fhash;
            //uint64_t* v = rhash;

            if(KHashMap.find(fwd) != KHashMap.end()){
                 u = KHashMap[fwd];
                 v = KHashMap[rev];
                
            }else{ 
                 u = hashSequence(fwd);
                 v = hashSequence(rev);

                 KHashMap[fwd] = u;
                 KHashMap[rev] = v;

            }

            //uint64_t u = hashSequence(kmerSeq(sequence, i + j, k, 0));
            //uint64_t v = hashSequence(kmerSeq(sequence, i + j, k, 1));

            if (u != v){
                m = std::min(m, std::min(u, v));
            }
                
        }
        for ( int j = 0; j < w; j++){
            //std::string forward = kmerSeq(sequence, i + j, k, 0);
            //std::string reverse = kmerSeq(sequence, i + j, k, 1);
            //const char* fwd_kseq = forward.c_str();
            //const char* rev_kseq = reverse.c_str();


           // uint64_t u = hash(forward);
           // uint64_t v = hash(reverse);

            //MurmurHash3_x64_128(fwd_kseq, strlen, 42, fhash);
            //MurmurHash3_x64_128(rev_kseq, strlen, 42, rhash);

            //uint64_t* u = fhash;
            //uint64_t* v = rhash;


            //uint64_t u,v;
            //std::string kseq = kmerSeq(sequence, i + j, k, 0);
            //std::string rev_kseq = kmerSeq(sequence, i + j, k, 1);
            /*if(KHashMap.find(forward) != KHashMap.end()){
                 u = KHashMap[forward];
                 v = KHashMap[reverse];
                
            }else{
                 u = KHashMap[forward];
                 v = KHashMap[reverse];
                
            }*/
        /*   uint64_t u = hashSequence(kmerSeq(sequence, i + j, k, 0));
            uint64_t v = hashSequence(kmerSeq(sequence, i + j, k, 1));
        
            if (u < v && u == m)
               hashes.emplace_back(MMHash(u,i+j,0));
            else if ( v < u && v ==m)
               hashes.emplace_back(MMHash(v,i+j,1));
        }
    }
    return hashes ;


}*/


/* Build  minimizer sketch for each target sequence
  param target_seq: target sequence 
  param int kmer: kmer size
  param int window : window size
  return a HashMap, with key as hash values and values as {targetSeq, position, strand}
 */

void MMSketch::buildSketch(){
    assert(is_kmer_set_);
    assert(is_window_set_);
    int count = 0;
    std::vector<MMIndex> mm_tuples;
    std::vector<MMIndex> ret;
    clock_t t;
    t = clock();
    for(FastaInfo target_entry: seqs){
        //for (auto &c : target_entry.sequence_) c = std::toupper(c);
        //mm_tuples = kmerTuples(target_entry.sequence_, k);
        //std::cout<<target_entry.header_<<"\n";
        //std::vector<MMHash> mm_hashes = computeMinimizerSketch(target_entry, k, w);
        computeMinimizerSketch(target_entry, k, w, mm_hashes);     
        count++;
  }
  t = clock() - t;
  std::cout<<"Collected minimizers in "<<((double)t)/CLOCKS_PER_SEC << " s\n";
  
  t=clock();
  std::sort(mm_hashes.begin(), mm_hashes.end(), hash_comparer());
  t = clock() - t;
  std::cout<<"Sorted minimizers in "<<((double)t)/CLOCKS_PER_SEC << " s\n";
  //PrintMinimizers();

  
  t=clock();
  // Inserting into hash table
 /* for( MMHash mm_entry : mm_hashes ){
        //    std::cout<<mm_entry.x<<"\t"<<mm_entry.y<<"\n";
            if( MMHashTable.find(mm_entry.x) != MMHashTable.end()){
                if(MMHashTable[mm_entry.x].back().info_ != mm_entry.y ){
                    MMHashTable[mm_entry.x].emplace_back(MMIndex(mm_entry.y));
                    std::sort(MMHashTable[mm_entry.x].begin(),MMHashTable[mm_entry.x].end(), comparer() );
                }   

            }
            else{
                MMHashTable[mm_entry.x].emplace_back(MMIndex(mm_entry.y));
            }

    }
  */  
    // Reserve memory to reduce reallocations
   MMHashTable.reserve(mm_hashes.size());

    for (const MMHash& mm_entry : mm_hashes) {
        auto it = MMHashTable.find(mm_entry.x);
        if (it != MMHashTable.end()) {
            if (it->second.empty() || it->second.back().info_ != mm_entry.y) {
                it->second.emplace_back(MMIndex(mm_entry.y));
            }
            if (!std::is_sorted(it->second.begin(), it->second.end(), comparer())) {
                std::sort(it->second.begin(), it->second.end(), comparer());
            }
        } else {
            MMHashTable[mm_entry.x] = { MMIndex(mm_entry.y) };
        }
    }
  t = clock() - t;
  std::cout<<"Created hash table index in "<<((double)t)/CLOCKS_PER_SEC << " s\n";

  //std::cout<<"Completed building Minimizer index for "<<count<<" target sequences"<<"\n";


}

void MMSketch::PrintIndexTable(void){
    //std::cout<<"Inside index table print\n";
    /*for(auto &elem: MMHashTable ){
        std::cout<<elem.first<<"\n";
        for(MMIndex &item : elem.second){
            std::cout<<"("<<item<<"), ";
        }
        std::cout<<"\n";
    }*/
    return;
}

void MMSketch::PrintMinimizers(void){
    std::cout<<"Inside mm_hashes print\n";
    for(auto &elem: mm_hashes ){
        std::cout<<elem.x<<"\t"<<elem.y<<"\n";
    }
    return;
}

void MMSketch::checkQuery(std::vector<MMHash> & qhash, std::map<uint64_t, std::vector<MMIndex>>& mhash_table)
{
    //for(MMHash &q : qhash){
        //std::cout<<q.hashValue<<"\n";
        /*for(auto it = MMHashTable.begin(); it != MMHashTable.end(); ++it ){
            std::cout<<(*it).first<<"\n";
            for(MMIndex &item : (*it).second){
                std::cout<<item.targetSeq<<"\t"<<item.position<<"\t"<<item.strand<<"\n";
            }
        }*/
        /*if(mhash_table.find(q.hashValue)!= mhash_table.end()){
            std::cout<<"Hash index exists\n";
        }
        else    std::cout<<"Something is wrong: Hash does not exist\n";*/

   // }

}


std::vector<MMMatch> MMSketch::slice(std::vector<MMMatch>const & matches, int start, int end)
{  
    auto first = matches.cbegin() +start;
    auto last = matches.cbegin() + end + 1;

    std::vector<MMMatch> sub_matches(first, last);

    return sub_matches;
}

std::vector<int> MMSketch::longest_increasing_subset(std::vector<MMMatch>& matches, std::string orientation)
{
    std::vector<int> map_len;
    for(MMMatch m : matches)    map_len.push_back(m.end);
    if(orientation == "-")  std::reverse(map_len.begin(), map_len.end());

    std::vector<int> P;
    std::vector<int> P_reversed;

    int n = map_len.size();
    for (int i = 0; i < n; i++) {
        auto it = lower_bound(P.begin(), P.end(), map_len[i]);
        if (it == P.end()) {
            P.push_back(map_len[i]);
        }
        else {
            *it = map_len[i];
        }
    }
    //for(int i=1; i<map_len.size(); i++){
    //Begin binary search
    /*for(MMMatch x : matches){
        if(P.empty() || P[P.size()-1].end < x.end){
            P.emplace_back(MMMatch(x.targetSeq, x.strand, x.start, x.end));
        }
        else{
            auto it = std::lower_bound(P.cbegin(), P.cend(), x.end);
            *it = x;
        }
    }*/

   /* M[0] = 0;
    int j, mid;
    int low = 0;
    int up = M.size();

    if( map_len[M[up-1]] < map_len[i]){
        j = up;
        M.insert(i);
    }
    else{
        while((up - low) > 1){
             mid = (low + up)/2;
              if(map_len[M[mid-1]]< map_len[i]){
                low = mid;
              }
              else  up = mid;
        }
        j = low;

        if(map_len[i]< map_len[M[j]]){
            M[j]= i;
        }
    }

    if(j==0) continue;
    else P.append(M[j-1]);

// NEEDS TO BE MODIFIED
    /*# Trace back using the predecessor array (P).
    def trace(i):
        if i is not None:
            yield from trace(P[i])
            yield i

    indices = list(trace(M[-1]))*/
    if(orientation == "+" || orientation == "-")  return P;
  /*  else{
        for(auto x: P)  P_reversed.push_back(map_len.size() - x - 1);
        std::reverse(P_reversed.begin(), P_reversed.end());
        return P_reversed;      
    } */
    

    
}

std::vector<std::pair<char, int>> MMSketch :: get_cigar(const std::string& s1, const std::string& s2, const std::vector<Cell>& traceback) {
    std::vector<std::pair<char, int>> cigar;
    int i = s1.size();
    int j = s2.size();
    for (int k = traceback.size() - 1; k >= 0; k--) {
        int i_prev = traceback[k].i;
        int j_prev = traceback[k].j;
        //std::cout<<i_prev<<"\t"<<j_prev<<"\n";
        if (i_prev == i) {
            cigar.emplace_back('D', j_prev - j);
        } else if (j_prev == j) {
            cigar.emplace_back('I', i_prev - i);
        } else {
            if (s1[i_prev-1] == s2[j_prev-1]) {
                cigar.emplace_back('M', 1);
            } else {
                cigar.emplace_back('X', 1);
            }
        }
        i = i_prev;
        j = j_prev;
    }
    reverse(cigar.begin(), cigar.end());
    
    return cigar;
}



std::vector<std::pair<char, int>>MMSketch :: x_drop_align(const std::string& s1, const std::string& s2) {
    int n = s1.size();
    int m = s2.size();
    int max_score = 0;
    int curr_score = 0;
    vector<Cell> traceback;
    vector<int> curr_row(n+1), prev_row(n+1);

    for (int i = 0; i <= n; i++) {
        curr_row[i] = -i * GAP_PENALTY;
    }

    for (int j = 1; j <= m; j++) {
        int max_score_col = -1;
        prev_row.swap(curr_row);
        curr_row[0] = -j * GAP_PENALTY;
        for (int i = 1; i <= n; i++) {
            int match_score = (s1[i-1] == s2[j-1]) ? 1 : -1;
            curr_score = std::max({prev_row[i-1] + match_score, curr_row[i-1] - GAP_PENALTY, prev_row[i] - GAP_PENALTY});
            curr_row[i] = curr_score;
            if (curr_score > max_score) {
                max_score = curr_score;
                max_score_col = i;
            }
        }
        if (max_score - curr_score > X_DROP_THRESHOLD) {
            break;
        }
        for (int i = 1; i <= n; i++) {
            if (i < max_score_col - X_DROP_THRESHOLD || i > max_score_col + X_DROP_THRESHOLD) {
                curr_row[i] = -1000000000;
            }
        }
        
        traceback.push_back({max_score, max_score_col, j});
    }
    //std::cout<<max_score<<"\n";
    
    //for(auto &x: traceback){
    //    std::cout<<x.score<<"\t"<<x.i<<"\t"<<x.j<<"\n";
    //}

    auto cigar = get_cigar(s1, s2, traceback);
    return cigar;
}


// Function to perform X-drop alignment for sequence extension
void MMSketch::xDropAlignmentExtend(const std::string& seq1, const std::string& seq2,
                          int q_start, int q_end, int t_start, int t_end,
                          std::string& aligned_seq1, std::string& aligned_seq2) {

    // Extend the alignment from the known overlap positions
    int q_len = q_end - q_start + 1;
    int t_len = t_end - t_start + 1;

    // Initialize the scoring matrix
    std::vector<std::vector<int>> score(q_len + 1, std::vector<int>(t_len + 1, 0));
    // Fill in the scoring matrix using X-drop algorithm
    for (int i = 1; i <= q_len; ++i) {
        for (int j = 1; j <= t_len; ++j) {
            int q_pos = q_start + i - 1;
            int t_pos = t_start + j - 1;

            int match = score[i - 1][j - 1] + (seq1[q_pos] == seq2[t_pos] ? match_score : mismatch_score);
            int gap1 = score[i - 1][j] + gap_score;
            int gap2 = score[i][j - 1] + gap_score;

            score[i][j] = std::max({0, match, gap1, gap2}); // X-drop threshold check
        }
    }

    int max_score = 0;
    int best_i = 0, best_j = 0;
    for (int i = 1; i <= q_len; ++i) {
        if (score[i][t_len] > max_score) {
            max_score = score[i][t_len];
            best_i = i;
            best_j = t_len;
        }
    }
    for (int j = 1; j <= t_len; ++j) {
        if (score[q_len][j] > max_score) {
            max_score = score[q_len][j];
            best_i = q_len;
            best_j = j;
        }
    }
    // Traceback 
    aligned_seq1 = "";
    aligned_seq2 = "";
    int i = best_i, j = best_j;
    while (i > 0 && j > 0 && score[i][j] > 0) {
        int q_pos = q_start + i - 1;
        int t_pos = t_start + j - 1;

        if (score[i][j] == score[i - 1][j - 1] + (seq1[q_pos] == seq2[t_pos] ? match_score : mismatch_score)) {
            aligned_seq1 = seq1[q_pos] + aligned_seq1;
            aligned_seq2 = seq2[t_pos] + aligned_seq2;
            i--;
            j--;
        } else if (score[i][j] == score[i - 1][j] + gap_score) {
            aligned_seq1 = seq1[q_pos] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
            i--;
        } else {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = seq2[t_pos] + aligned_seq2;
            j--;
        }
    }
}

/*
Mapping each query sequence to minimizer index
*/

std::string MMSketch::Map(FastaInfo& seq_entry, int& kmer, int& window, int epsilon)
{
    std::vector <MMMatch> query_matches;
    std::vector<MMMatch> potential_olap;
    std::vector<int> indices;
    std::string orientation = "";

    std::string asqg_out = "";
    std::string qname = seq_entry.header_;

    int tstart, tend, qstart, qend, orient, mapped_bp;
    std::string tname, first_target, last_target;
    int first_start = 0;
    int first_end = 0;
    int last_end = 0;
    int last_start = 0;

    

    // Get mm sketch for each query sequence
    std::vector<MMHash> qhash ;
    computeMinimizerSketch(seq_entry, kmer, window, qhash);
    for(MMHash q : qhash){
        try{
            for(MMIndex t: MMHashTable[q.x]){
                
                uint32_t qid = q.y>>32;
                int qstrand = q.y & 1;
                uint32_t qpos = (q.y >> 1) & 0x7fffffff;

                uint32_t tid = t.info_>>32;
                int tstrand = t.info_ & 1;
                uint32_t tpos = (t.info_ >> 1) & 0x7fffffff;
                
                //if(seq_entry.header_ == "Read_1591")
                   // if( tid == 29309)
                       // std::cout<<t.info_<<"\t"<<tid<<"  "<<qid<<"  "<<tstrand<<"\n";
               
                if(tid == qid)  
                    continue;
                if( tstrand == qstrand){
                    query_matches.emplace_back(MMMatch(tid, 0, qpos - tpos, tpos));
                } 
                else{
                    query_matches.emplace_back(MMMatch(tid, 1, qpos + tpos, tpos));
                }

            }          
        }
        catch(const std::out_of_range& oor) {};
    }
    std::sort(query_matches.begin(), query_matches.end(), qcomparer());
     /*if(seq_entry.header_ == "Read_3181"){
                for( auto &q: query_matches){
                    if(q.id == 1815){
                        std::cout<<q.id<<"\t"<<q.strand<<"\t"<<q.start<<"\t"<<q.end<<"\n";
                    }
                }   
        }*/

    int b=0;

    if(query_matches.size()>0){
        for(int e = 0; e < query_matches.size()-1; ++e){
            if(e == query_matches.size()-1 | 
                query_matches[e+1].id != query_matches[e].id | 
                query_matches[e+1].strand != query_matches[e].strand | 
                query_matches[e+1].start - query_matches[e].start >= epsilon)
                {   
                  
                    potential_olap = slice(query_matches, b, e);

                    /*if(seq_entry.header_ == "Read_34"){
                        //for( auto &q: query_matches){
                        //if(q.id == 1590){
                        //    std::cout<<q.id<<"\t"<<q.strand<<"\t"<<q.start<<"\t"<<q.end<<"\n";
                        //}
                        //}          
                        for(auto &e: potential_olap){
                            if(e.id == 1476)
                            std::cout<<e.id<<"\t"<<e.start<<"\t"<<e.end<<"\n";
                        }}*/
                 
                    if(query_matches[e].strand == 0){
                        orientation = "+"; 
                        /*tstart = query_matches[b].start+query_matches[b].end;
                        tend = query_matches[e].start+query_matches[e].end+k;
                        qstart = query_matches[b].end;
                        qend = query_matches[e].end+k;
                        tname = query_matches[e].targetSeq;
                        orient = 0;
                        target_indices.emplace_back(MMMatch(query_matches[e].targetSeq, 0, tstart, tend));
                        query_indices.emplace_back(MMMatch(qname, 0, qstart, qend));
                        */
                    } 
                    
                    else if (query_matches[e].strand == 1){
                        orientation = "-";
                        /*tstart = query_matches[b].start-query_matches[b].end;
                        tend = query_matches[e].start-query_matches[e].end+k;
                        qstart = query_matches[b].end;
                        qend = query_matches[e].end+k;
                        orient = 1;
                        target_indices.emplace_back(MMMatch(query_matches[e].targetSeq, 1, tstart, tend));
                        query_indices.emplace_back(MMMatch(qname, 1, qstart, qend ));
                        std::cout<<tstart<<"\t"<<tend<<"\t"<<orient<<"\n";
                        */
                        
                    }
                    //std::cout<<qname<<"\t"<<orientation<<"\n";
                
                    b = e+1;
                
                    if(orientation == "+"){
                        indices = longest_increasing_subset(potential_olap, orientation);   
                       // std::cout<<orientation<<"\n";
                     /*if(seq_entry.header_ == "Read_3181"){
                       // std::cout<<"inside or\n";
                        for(auto &p : potential_olap)
                            {
                                if(p.id == 1815){
                                    for(auto x: indices)
                                        std::cout<<x<<"\t";
                                    std::cout<<"\n";
                                }   
                            }
                       }*/
                    }
                    else{
                        std::reverse(potential_olap.begin(), potential_olap.end());
                        /*if(seq_entry.header_ == "NC_015678.1_1743835_1744360_2:0:0_1:0:0_2/1"){
                        for(auto &p : potential_olap)
                            {
                                if(p.id == 29309){
                                        std::cout<<p.id<<"\t"<<p.strand<<"\t"<<p.start<<"\t"<<p.end<<"\n";
                                    
                                }   
                            }
                       }*/
                        indices = longest_increasing_subset(potential_olap, orientation); 
                        /*if(seq_entry.header_ == "NC_015678.1_1743835_1744360_2:0:0_1:0:0_2/1"){
                        for(auto &p : potential_olap)
                            {
                                if(p.id == 29309){
                                    for(auto x: indices)
                                        std::cout<<x<<"\t";
                                    std::cout<<"\n";
                                }   
                            }
                       } */
                       // std::cout<<orientation<<"\n";
                        //for(auto x: indices)
                          //  std::cout<<x<<"\t";
                        //std::cout<<"\n";   
                    }
                    if(indices.size() <= 3) continue;

                    int first_mm = indices[0];
                    int last_mm = indices[indices.size()-1];

                   // std::cout<<first_mm<<"\t"<<last_mm<<"\n";

                    // testing
                    //MMMatch first_mm = potential_olap[indices[0]];
                    //MMMatch last_mm = potential_olap[indices[indices.size() - 1]];

                    /*if(seq_entry.header_ == "Read_34"){
                        for(auto &p : potential_olap)
                            {
                                if(p.id == 1476){
                                    std::cout<<indices[0]<<"\t"<<indices[indices.size()-1]<<"\n";
                                    std::cout<<"First mm "<<first_mm<<"\t"<<"Last mm "<<last_mm<<"\n";
                                }   
                            }
                            
                    }*/
                    //if(seq_entry.header_ == "Read_9311"){
                        for(MMMatch x : potential_olap){
                           // if(x.id == 1590){
                                if(x.end == first_mm){
                                //std::cout<<x.end<<"\t"<<first_mm<<"\n";
                                tname = id_map[x.id];
                                first_start = x.start;
                                first_end = x.end;
                            }
                            else if( x.end == last_mm){
                                //std::cout<<x.end<<"\t"<<last_mm<<"\n";
                                last_start = x.start;
                                last_end = x.end;
                            }
                            //std::cout<<tname<<"\t"<<first_start<<"\t"<<first_end<<"\t"<<last_start<<"\t"<<last_end<<"\n";
                            }

                            
                            //}
                        
                   // }
                    
                    //}
                        
                    int mapped_bp = indices.size() * kmer;
                    int qlen = seq_entry.seq_len_;
                    int tlen = seq_entry.seq_len_;
                    //tname = id_map[last_mm.id];
                   

                   /*if(orientation == "+"){
                    tstart = first_mm.end;
                    tend = last_mm.end + kmer;
                    qstart = first_mm.start + first_mm.end;
                    qend = last_mm.start + last_mm.end + kmer;

                   }
                   else{
                    tstart = tlen - first_mm.end - kmer;
                    tend = tlen - last_mm.end;
                    qstart = first_mm.start - first_mm.end;
                    qend = last_mm.start - last_mm.end + kmer;
                   }*/

                   //latest implementation
                /* if(orientation == "+"){
                    tstart = first_start;
                    tend = last_end;
                    qstart = first_end + first_start - kmer;
                    qend = last_end + last_start;

                   }
                   else{
                    tstart = first_end - kmer;
                    tend = last_end;
                    qstart = first_end - kmer;
                    qend = last_end ;
                   }
                */
                  //latest implementation
                 if(orientation == "+"){
                    tstart = first_start;
                    //tend = last_end + kmer;
                    tend = last_end;
                    qstart = first_end + first_start - kmer; 
                    //qstart = first_end + first_start - kmer;
                    qend = last_end + last_start;

                   }
                   else{
                    tstart = first_end - kmer;
                    tend = last_end;
                    qstart = first_end - kmer;
                    qend = last_end ;
                   }


                 
                       

                   /* if(orientation == "+" && first_start >=0 && first_end >=0 && last_start >= 0 && last_end >=0){
                        std::cout<<orientation<<"\t"<<first_start<<"\t"<<first_end<<"\t"<<last_start<<"\t"<<last_end<<"\t"<<kmer<<"\n";
                        tstart = tlen - first_end - kmer;
                        tend = tlen - last_end;
                        qstart = first_start + first_end;
                        qend = last_start - last_end + kmer;
                    }
                    else{
                        if(first_start >=0 && first_end >=0 && last_start >= 0 && last_end >=0)
                        {
                            std::cout<<orientation<<"\t"<<first_start<<"\t"<<first_end<<"\t"<<last_start<<"\t"<<last_end<<"\t"<<kmer<<"\n";
                            tstart = first_end;
                            tend = last_end + kmer;
                            qstart = first_start - first_end;
                            qend = last_start + last_end +  kmer;
                        }
 
                    }

                    /*if( orientation == "-"){
                        tstart = tlen - tend;
                        tend = tlen - tstart;
                    }*/

                    

                    int overhang = std::min(std::abs(qstart),std::abs(tstart)) + std::min(qlen-qend, tlen-tend);
                    int maplen = std::max(qend-std::abs(qstart), tend - std::abs(tstart));
                    if(maplen < MAPLEN) continue;
                    
                   /* if( seq_entry.header_ == "Read_34"){
                         std::cout<<mapped_bp<<"\t"<<qlen<<"\t"<<last_start<<"\t"<<last_end<<"\t"<<first_start<<"\t"<<first_end<<"\t"<<overhang<<"\t"<<maplen<<"\n";
                   }

                    if(seq_entry.header_ == "Read_34"){
                    std::cout<<qname<<"\t"<<qlen<<"\t"<<qstart<<"\t"<<qend<<"\t"<<orientation<<"\t"<<tname<<"\t"<<tlen<<"\t"<<std::abs(tstart)<<"\t"<<tend<<"\t"<<maplen<<"\n";
                    }*/
                    
                    
                    if(overhang > std::min(MINHANG, int(0.8*maplen))){
                     //   if(seq_entry.header_ == "Read_34")
                       // std::cout<<"overhang: "<<overhang<<"\t"<<"maplen :"<<maplen<<"\n";
                        continue;
                    }
                        
                    //if(std::abs(qstart) <= std::abs(tstart) && (qlen - qend) <= (tlen-tend))
                    //    continue;
                    //else if(std::abs(qstart) >= std::abs(tstart) && (qlen - qend) >= (tlen-tend))
                    //    continue;
                    //if(seq_entry.header_ == "Read_34")
                    //std::cout<<qname<<"\t"<<qlen<<"\t"<<std::abs(qstart)<<"\t"<<qend<<"\t"<<orientation<<"\t"<<tname<<"\t"<<tlen<<"\t"<<std::abs(tstart)<<"\t"<<tend<<"\n";
                    
                    //if(qstart >= tstart){
                        //std::cout<<"QSUB: "<<qname<<"\t"<<seq_entry.sequence_<<"\t"<<abs(qstart)<<"\t"<<abs(qend)<<"\n";
                        //std::cout<<"TSUB: "<<tname<<"\t"<<TargetSeqMap[tname]<<"\t"<<abs(tstart)<<"\t"<<abs(tend)<<"\n";
                        //std::cout<<qstart<<"\t"<<qlen<<"\t"<<"0"<<"\t"<<tend<<"\n";
                        //std::string q_sub = seq_entry.sequence_.substr(abs(qstart), qlen);
                        //std::string t_sub = TargetSeqMap[tname].substr(0, tend);
                        //std::cout<<q_sub<<"\t"<<t_sub<<"\n";
                        //auto cigar = x_drop_align(q_sub, t_sub);
                        //for (auto& p : cigar) {
                        // cout << p.second << p.first << " ";
                        //}
                        //std::cout<<"\n";
                        //std::cout<<q_sub<<"\t"<<t_sub<<"\n";
                        std::string aligned_seq1, aligned_seq2;
                        xDropAlignmentExtend(seq_entry.sequence_, TargetSeqMap[tname], abs(qstart), qend, abs(tstart), tend, aligned_seq1, aligned_seq2);

                        // Find the start and end positions of the overlap in the aligned sequences
                        int overlap_start = -1;
                        int overlap_end = -1;

                        for (int i = 0; i < aligned_seq1.length(); ++i) {
                            if (aligned_seq1[i] != '-' && aligned_seq2[i] != '-') {
                                if (overlap_start == -1) {
                                    overlap_start = i;
                                }
                                overlap_end = i;
                            }
                        }
                        // Print the extended aligned sequences
                        //std::cout << "Aligned sequence 1: " << aligned_seq1 << std::endl;
                        //std::cout << "Aligned sequence 2: " << aligned_seq2 << std::endl;
                       /* if (overlap_start != -1 && overlap_end != -1) {
                            std::cout << "Overlap region in query: " << abs(qstart) + overlap_start << " to " << abs(qstart) + overlap_end << std::endl;
                            std::cout << "Overlap region in target: " << abs(tstart) + overlap_start << " to " << abs(tstart) + overlap_end << std::endl;
                        } else {
                            std::cout << "No overlap found in the alignment." << std::endl;
                        }*/

                        //asqg_out += "ED\t"+qname+" "+tname+" "+qstart+" "+qend+" "+qlen+" "+tstart+" "+tend+" "+tlen+" "+orient+" "+0+"\n";
                        //if(qstart<=10 && (tlen-tend)<=10)
                            //std::cout<<"ED\t"<<qname<<"\t"<<tname<<" "<<qstart<<" "<<qend<<" "<<qlen<<" "<<tstart<<" "<<tend<<" "<<tlen<<" "<<"+"<<" 0"<<"\n";
                        //std::cout<<qname<<"\t"<<qlen<<"\t"<<qstart<<"\t"<<qend<<"\t"<<orientation<<"\t"<<tname<<"\t"<<tlen<<"\t"<<tstart<<"\t"<<tend<<"\n"; 
                    //}
                   /* else{
                        //std::cout<<"QSUB: "<<qname<<"\t"<<seq_entry.sequence_<<"\t"<<abs(qstart)<<"\t"<<qend<<"\t"<<"\n";
                       // std::cout<<"TSUB: "<<TargetSeqMap[tname]<<"\t"<<abs(tstart)<<"\t"<<tend<<"\n";
                        //std::cout<<qstart<<"\t"<<qlen<<"\t"<<"0"<<"\t"<<tend<<"\n";
                        //std::string q_sub = seq_entry.sequence_.substr(abs(qstart), qlen);
                        //std::string t_sub = TargetSeqMap[tname].substr(0, tend);
                        
                        std::string aligned_seq1, aligned_seq2;
                        xDropAlignmentExtend(seq_entry.sequence_, TargetSeqMap[tname], abs(qstart), qend, abs(tstart), tend, aligned_seq1, aligned_seq2);
                        //auto cigar = x_drop_align(q_sub, t_sub);
                        //for (auto& p : cigar) {
                        // cout << p.second << p.first << " ";
                       //}
                        //std::cout<<"\n";
                        //std::cout<<q_sub<<"\t"<<t_sub<<"\n";
                        // Find the start and end positions of the overlap in the aligned sequences
                        int overlap_start = -1;
                        int overlap_end = -1;

                        for (int i = 0; i < aligned_seq1.length(); ++i) {
                            if (aligned_seq1[i] != '-' && aligned_seq2[i] != '-') {
                                if (overlap_start == -1) {
                                    overlap_start = i;
                                }
                                overlap_end = i;
                            }
                        }
                        // Print the extended aligned sequences
                        //std::cout << "Aligned sequence 1: " << aligned_seq1 << std::endl;
                        //std::cout << "Aligned sequence 2: " << aligned_seq2 << std::endl;
                        /*if (overlap_start != -1 && overlap_end != -1) {
                            std::cout << "Overlap region in query: " << abs(qstart) + overlap_start << " to " << abs(qstart) + overlap_end << std::endl;
                            std::cout << "Overlap region in target: " << abs(tstart) + overlap_start << " to " << abs(tstart) + overlap_end << std::endl;
                        } else {
                            std::cout << "No overlap found in the alignment." << std::endl;
                        }*/
                        //continue;
                        //asqg_out += "ED\t"+qname+" "+tname+" "+qstart+" "+qend+" "+qlen+" "+tstart+" "+tend+" "+tlen+" "+orient+" "+1+"\n";
                        //if(tstart<=10 && (qlen-qend)<=10)
                            //std::cout<<"ED\t"<<qname<<" "<<tname<<" "<<qstart<<" "<<qend<<" "<<qlen<<" "<<tstart<<" "<<tend<<" "<<tlen<<" "<<"-"<<" 0"<<"\n"; 
                        //std::cout<<qname<<"\t"<<qlen<<"\t"<<qstart<<"\t"<<qend<<"\t"<<orientation<<"\t"<<tname<<"\t"<<tlen<<"\t"<<tstart<<"\t"<<tend<<"\n";
                    //}

                    
                    
                    //asqg_out += "ED\t"+qname+" "+tname+" "+qstart+" "+qend+" "+qlen+" "+tstart+" "+tend+" "+tlen+" "+orient+" "+0+"\n";               
                    //std::cout<<"ED\t"<<qname<<" "<<tname<<" "<<qstart<<" "<<qend<<" "<<qlen<<" "<<tstart<<" "<<tend<<" "<<tlen<<" "<<orientation<<" 0"<<"\n"; 
                
            }
        }
    } 
    return asqg_out;
}  


