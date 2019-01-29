#include "bitmorev5.h"
#include "constexpr_for.h"
#include <algorithm>
#include <sys/time.h>
#include <iostream>
#include <string>
#include <chrono>

double getTimeElapsedMilliSec(struct timeval end, struct timeval start)
{
  //tv_usec is microseconds
  //tv_sec is in sec
  return ((end.tv_sec - start.tv_sec)*1000.00) + ((end.tv_usec - start.tv_usec)/1000.00);//milliseconds
}

double getTimeElapsedMicroSec(struct timeval end, struct timeval start)
{
  //tv_usec is microseconds
  //tv_sec is in sec
  return ((end.tv_sec - start.tv_sec)*1000000.00) + ((end.tv_usec - start.tv_usec));//microseconds
}

template <size_t nbytes_per_word>
struct word
{ 
  static constexpr size_t len256 = std::ceil(nbytes_per_word / static_cast<double>(sizeof(__m256i)));
  __m256i m256[len256];
  inline word<nbytes_per_word> & operator^=(const word<nbytes_per_word> & other)
  {
    for (size_t i = 0; i < len256; ++i) m256[i] = _mm256_xor_si256(m256[i], other.m256[i]);
    return *this;
  }
};

template <size_t nbytes_per_word, size_t nwords_per_row>
struct record
{
private:
  word<nbytes_per_word> words[nwords_per_row];
public:
  record(word<nbytes_per_word> & word) { std::fill_n(words, nwords_per_row, word); }
  inline constexpr word<nbytes_per_word> & operator[](const size_t j) { return words[j]; }
};

int main(int argc, char * argv[])
{
  using namespace dpf;
  int ntrials = 8, ntrials_pruned = 5, nmetrics = 4, confidence = 2;
  char* metrics[] = {"All query generation (us):\t", "Query expansion per server (us):", "Response generation per server (ms):", "Record reconstruction (ns):\t"};
  struct timeval tvalBefore, tvalAfter;
  FILE *f = fopen("result_chor_Jan27.txt", "w");
  // FILE *f1 = fopen("result_chor_details_Jan27.txt", "w");
  int err;

  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  constexpr std::array<uint8_t, 7> Lognitems = {22, 20, 18, 16, 14, 12, 10};//22, 20, 18, 16, 14, 12, 
  constexpr std::array<uint8_t, 5> Nwords_per_row = {31, 15, 7, 3, 1};
  for_constexpr<for_bounds<0, 5>, for_bounds<0, 7>>(
      [&](auto iw, auto ir)
  {
    constexpr uint8_t lognitems = std::get<ir.value>(Lognitems);
    constexpr size_t nitems = (1ULL << lognitems);
    constexpr size_t nbytes_per_row = 512;
    constexpr size_t item = 11;

    constexpr uint8_t nservers = 2;
    constexpr size_t nwords_per_row = std::get<iw.value>(Nwords_per_row);
    constexpr size_t nbytes_per_word = std::ceil(nbytes_per_row / static_cast<double>(nwords_per_row));

    typedef record<nbytes_per_word, nwords_per_row> record;

    size_t alloc_size = nitems * sizeof(record);

    void * database_;
    err = posix_memalign((void**)&database_, sizeof(__m256i), alloc_size);
    if(err) perror("Error in memalign for database_");
    arc4random_buf(database_, alloc_size);
    
    //printf("***********************************\nEach row: %f KiB, Database: %f GiB.\n", (float)nbytes_per_row/(1024), (float)(nitems*nbytes_per_row)/(1024*1024*1024));
    fprintf(f,"***********************************\nEach row: %.3f KiB, Database: %.5f GiB, In memory: %.5f GiB.\n", (float)nbytes_per_row/(1024), (float)(nitems*nbytes_per_row)/(1024*1024*1024),(float)(alloc_size)/(1024*1024*1024));    
    //fprintf(f1,"***********************************\nEach row: %.3f KiB, Database: %.5f GiB, In memory: %.5f GiB.\n", (float)nbytes_per_row/(1024), (float)(nitems*nbytes_per_row)/(1024*1024*1024),(float)(alloc_size)/(1024*1024*1024));    

    double sum[nmetrics];
    for (int i = 0; i < nmetrics; ++i) { sum[i] = 0.0;  }
      
    double time_elapses[nmetrics][ntrials];
    //printf("\t\trows: 2^%u, ell: %u\n", lognitems, nservers);
    fprintf(f,"\t\trows: 2^%u, nwords/row: %zu, ell: %u\n", lognitems, nwords_per_row, nservers);
    //fprintf(f1,"\t\trows: 2^%u, nwords/row: %zu, ell: %u\n", lognitems, nwords_per_row, nservers);

    for(int trials=0;trials<ntrials;++trials)
    {
      __m128i ** s = (__m128i**)malloc(nservers * sizeof(__m128i *));
      uint8_t ** t = (uint8_t**)malloc(nservers * sizeof(uint8_t *));
      for (size_t i = 0; i < nservers; ++i)
      {
        err = posix_memalign((void**)&s[i], sizeof(__m256i), dpf::dpf_key<nitems>::output_length * sizeof(__m128i));
        if(err) perror("Error in memalign for s");
        t[i] = (uint8_t *)malloc(dpf::dpf_key<nitems>::output_length * sizeof(uint8_t));
      }

      record * database = reinterpret_cast<record *>(database_);
    
      int query_length = std::ceil(nitems / 64);
      dpf_key<nitems> dpfkey[2];

      //Client: query generation
      gettimeofday(&tvalBefore, NULL);
      dpf::gen(aeskey, item, dpfkey);
      gettimeofday(&tvalAfter, NULL);
      time_elapses[0][trials] = getTimeElapsedMicroSec(tvalAfter, tvalBefore);
      //fprintf(f1,"(0) %f\n", time_elapses[0][trials]);

      //Query expansion at s1 th server
      gettimeofday(&tvalBefore, NULL);
      dpf::evalfull(aeskey, dpfkey[0], s[0], t[0]);
      uint64_t query[query_length];
      for(int i=0;i<query_length/2;i++) {
        query[i*2]  = s[0][i][0];
        query[i*2+1] = s[0][i][1];
      }
      gettimeofday(&tvalAfter, NULL);
      time_elapses[1][trials] = getTimeElapsedMicroSec(tvalAfter, tvalBefore);

      //response record from s1 th server
      gettimeofday(&tvalBefore, NULL);
      word<nbytes_per_word> zerow = { 0 };
      record result = { zerow };
      int record_count = 0;
      for (size_t k = 0; k < query_length; ++k) 
      {
        uint64_t bitset = query[k];
        // printf("%d\n", record_count);
        while (bitset != 0) {

          const int nextbit = __builtin_ctzll(bitset);//trailing zero, i.e. the row number where is a 1
          if( k * 64 + nextbit >= nitems ){break;}
          for(int p=0;p<nwords_per_row;p++) 
          {
            result[p] ^= database[k*64+nextbit][p];
          }
          // std::cout << "bitset: " << std::bitset<64>(bitset) << "\n";
          // printf("nextbit: %d\n", nextbit );
          
          bitset ^= bitset & -bitset;//vanishes LSB 1
          record_count++;
        }
      }
      gettimeofday(&tvalAfter, NULL);
      time_elapses[2][trials] = getTimeElapsedMilliSec(tvalAfter, tvalBefore);

      //Query expansion at s2 th server
      gettimeofday(&tvalBefore, NULL);
      dpf::evalfull(aeskey, dpfkey[1], s[1], t[1]);
      uint64_t query1[query_length];
      for(int i=0;i<query_length/2;i++) 
      {
        query1[i*2]  = s[1][i][0];
        query1[i*2+1] = s[1][i][1];
      }
      gettimeofday(&tvalAfter, NULL);
      time_elapses[1][trials] += getTimeElapsedMicroSec(tvalAfter, tvalBefore);
      
      //response record from s2 th server
      gettimeofday(&tvalBefore, NULL);
      record result1 = { zerow };
      record_count = 0;
      for (size_t k = 0; k < query_length; ++k) 
      {
        uint64_t bitset = query1[k];
        // printf("%d\n", record_count);
        while (bitset != 0)
        {
          const int nextbit = __builtin_ctzll(bitset);//trailing zero, i.e. the row number where is a 1
          if(k*64+nextbit>=nitems){break;}
          for(int p=0;p<nwords_per_row;p++) 
          {
            result1[p] ^= database[k*64+nextbit][p];
          }
          // std::cout << "bitset: " << std::bitset<64>(bitset) << "\n";
          // printf("nextbit: %d\n", nextbit );
          
          bitset ^= bitset & -bitset;//vanishes LSB 1
          record_count++;
        }
      }
      gettimeofday(&tvalAfter, NULL);
      time_elapses[2][trials] += getTimeElapsedMilliSec(tvalAfter, tvalBefore);
      //fprintf(f, "%lf\n",time_elapses[trials]);
      //sum[2]+=time_elapses[2][trials];

      time_elapses[1][trials]/=nservers;
      time_elapses[2][trials]/=nservers;
      //fprintf(f1,"(1) %f\n", time_elapses[1][trials]);
      //fprintf(f1,"(2) %f\n", time_elapses[2][trials]);
      sum[0]+=time_elapses[0][trials];
      sum[1]+=time_elapses[1][trials];
      sum[2]+=time_elapses[2][trials];


      //Client: record reconsteuction
      auto start = std::chrono::steady_clock::now();
      for(int p=0;p<nwords_per_row;p++) 
      {
        result[p] ^= result1[p];
      }
      auto end = std::chrono::steady_clock::now();
      time_elapses[3][trials] = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
      //fprintf(f1,"(3) %f\n", time_elapses[3][trials]);
      sum[3]+=time_elapses[3][trials];

      // for (int i = 0; i < nwords_per_row; ++i)
      // {
      //  printf("%llu=?=%llu\n", result[i].m256[0][0], database[item][i].m256[0][0]);
      // }

      bool correct_ness = true;
      for (size_t i = 0; i < nwords_per_row; ++i)
      {
        result[i] ^= database[item][i];
      }
      for (size_t i = 0; i < nwords_per_row; ++i)
      {
        if(result[i].m256[0][0]!=0x00) {correct_ness = false; printf("%u[%d]->Incorrect!\n",nservers,trials);break;}
      }
      if (correct_ness) printf("%zu, %u [%d]-> c\n", nwords_per_row, lognitems, trials);

      for (size_t i = 0; i < nservers; ++i)
      {
        free(s[i]);
        free(t[i]);
      }
      free(s);
      free(t);
    }
    double means[nmetrics], sum_for_sd[nmetrics], s_ds[nmetrics];
    for (int i = 0; i < nmetrics; ++i) { means[i]  = sum[i]/ntrials;  sum_for_sd[i]=0;}
    for (int i = 0; i < nmetrics; ++i) 
    {
      for(int tt=0;tt<ntrials;++tt){
        sum_for_sd[i]+=(time_elapses[i][tt]-means[i])*(time_elapses[i][tt]-means[i]);
      }
      s_ds[i] = std::sqrt(sum_for_sd[i]/(ntrials-1));
    }

    //pruned outliers
    double prune_upper[nmetrics], prune_lower[nmetrics];
    for(int i = 0; i < nmetrics; ++i)
    {
      prune_upper[i] = means[i]+confidence*s_ds[i];
      prune_lower[i] = means[i]-confidence*s_ds[i];
    }

    double time_elapses_pruned[nmetrics][ntrials_pruned];
    int support[nmetrics];
    memset(support,0,sizeof(support)*nmetrics);
    for(int i = 0; i < nmetrics; ++i)
    {
      for(int j=0;j<ntrials;++j)
      {
        if(time_elapses[i][j]<=prune_upper[i] && time_elapses[i][j]>=prune_lower[i] && support[i] < ntrials_pruned)
          time_elapses_pruned[i][support[i]++]=time_elapses[i][j];
      }
    }

    for (int i = 0; i < nmetrics; ++i) { sum[i] = 0.0;  }
    for (int i = 0; i < nmetrics; ++i) 
    {
      for(int trials=0;trials<support[i];++trials)
      {
        sum[i]+=time_elapses_pruned[i][trials];
      }
    }

    for (int i = 0; i < nmetrics; ++i) 
    { 
        means[i]  = sum[i]/support[i];
        sum_for_sd[i]=0;
    }
    //printf("m: %f\t",mean_time_elapse);
    //double sum_for_sd[nmetrics], s_ds[nmetrics];
    for (int i = 0; i < nmetrics; ++i) 
    {
      for(int tt=0;tt<support[i];++tt)
      {
        sum_for_sd[i]+=(time_elapses_pruned[i][tt]-means[i])*(time_elapses_pruned[i][tt]-means[i]);
      }
      s_ds[i] = std::sqrt(sum_for_sd[i]/(support[i]-1));
    }



    for(int ii=0;ii<nmetrics;++ii)
    {
      //printf("%s\tmean: %f\tsd: %f\tsupport: %i\n", metrics[ii], means[ii], s_ds[ii], support[ii]);
      fprintf(f,"%s\tmean: %.5f\tsd: %.5f\tsupport: %i\n", metrics[ii], means[ii], s_ds[ii], support[ii]);
      //fprintf(f1,"%s\tmean: %.5f\tsd: %.5f\tsupport: %i\n", metrics[ii], means[ii], s_ds[ii], support[ii]);
    }
  
  free(database_);
  }
);
  fclose(f);
  //fclose(f1);
  return 0;
  
}