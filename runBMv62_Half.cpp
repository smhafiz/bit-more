#include "bitmorev5.h"
#include "constexpr_for.h"
#include <algorithm>
#include <sys/time.h>
#include <iostream>
#include <string>

template <int begin, int end, typename lambda>
inline void static_for(lambda const & f)
{
  if constexpr (begin < end)
  {
    f(std::integral_constant<int, begin>{});
    static_for<begin + 1, end>(f);
  }
}

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

double getTimeElapsedNanoSec(struct timeval end, struct timeval start)
{
  //tv_usec is microseconds
  //tv_sec is in sec
  return ((end.tv_sec - start.tv_sec)*1000000000.00) + ((end.tv_usec - start.tv_usec) * 1000.00);//nanoseconds
}

int main(int argc, char * argv[])
{
  int ntrials = 150, ntrials_pruned = 100, nmetrics = 2, confidence = 2;
  std::string metrics[] = {"All query generation (ns):\t", "Query expansion per server (ms):"};
  struct timeval tvalBefore, tvalAfter;
  FILE *f = fopen("power2_half_Jan28.txt", "w");
  int err;

  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  constexpr std::array<uint8_t, 12> Lognitems = {32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10};
  
  for_constexpr<for_bounds<0, 2>, for_bounds<0, 12>>(
      [&](auto il, auto ir) 
      
  {
    constexpr uint8_t lognitems = std::get<ir.value>(Lognitems);
    constexpr size_t nitems = (1ULL << lognitems);
    constexpr size_t item = 101;

    uint8_t * expanded_query;
    err = posix_memalign((void**)&expanded_query, sizeof(__m256i), nitems * sizeof(uint8_t));
    if(err) 
      {
        perror("Error in memalign for expanded_query");
        //break;
      }

    constexpr std::array<uint8_t, 2> Nservers = {4, 2};

    constexpr uint8_t nservers = std::get<il.value>(Nservers);
    constexpr size_t nwords_per_row = nservers-1;
    constexpr size_t soundness = 0;
    constexpr size_t nkeys = std::ceil(std::log2(nservers)) + soundness;
    
    double sum[nmetrics];
    for (int i = 0; i < nmetrics; ++i) { sum[i] = 0.0;  }
      
    double time_elapses[nmetrics][ntrials];
    //printf("\t\trows: 2^%u, ell: %u\n", lognitems, nservers);
    fprintf(f,"\t\trows: 2^%u, ell: %u\n", lognitems, nservers);
    for(int trials=0;trials<ntrials;++trials)
    {
      __m128i ** s = (__m128i**)malloc(nkeys * sizeof(__m128i *));
      uint8_t ** t = (uint8_t**)malloc(nkeys * sizeof(uint8_t *));
      for (size_t i = 0; i < nkeys; ++i)
      {
        err = posix_memalign((void**)&s[i], sizeof(__m256i), dpf::dpf_key<nitems>::output_length * sizeof(__m128i));
        if(err) perror("Error in memalign for s");
        t[i] = (uint8_t *)malloc(dpf::dpf_key<nitems>::output_length * sizeof(uint8_t));
      }

//      typedef bitmore::record<nbytes_per_word, nwords_per_row> record;
      typedef bitmore::query<nitems, nkeys, nservers> query;

//      record * database = reinterpret_cast<record *>(database_);
    
      //Client: query generation
      gettimeofday(&tvalBefore, NULL);
      query q(aeskey, item);
      auto Q = q.get_request(nwords_per_row);
      gettimeofday(&tvalAfter, NULL);
      time_elapses[0][trials] = getTimeElapsedNanoSec(tvalAfter, tvalBefore);
      //fprintf(f, "%lf\n",time_elapses[trials]);
      sum[0]+=time_elapses[0][trials];


      //Query expansion at s+1 th server
      gettimeofday(&tvalBefore, NULL);
      bitmore::evalfull<nkeys,nitems,nservers>(aeskey, Q, s, t, expanded_query);
      gettimeofday(&tvalAfter, NULL);
      time_elapses[1][trials] = getTimeElapsedMilliSec(tvalAfter, tvalBefore);
      //fprintf(f, "%lf\n",time_elapses[trials]);
      //sum[1]+=time_elapses[1][trials];

      sum[1]+=time_elapses[1][trials];
      // sum[2]+=time_elapses[2][trials];

      printf("%u, %u [%d]->!\n",nservers, lognitems, trials);
      for (size_t i = 0; i < nkeys; ++i)
      {
        free(s[i]);
        free(t[i]);
      }
      free(s);
      free(t);
    }
    double means[nmetrics], sum_for_sd[nmetrics], s_ds[nmetrics];
    for (int i = 0; i < nmetrics; ++i)
    {
      means[i]  = sum[i]/ntrials; 
      sum_for_sd[i] = 0.0;
    }
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
      for(int tt=0;tt<support[i];++tt){
        sum_for_sd[i]+=(time_elapses_pruned[i][tt]-means[i])*(time_elapses_pruned[i][tt]-means[i]);
      }
      s_ds[i] = std::sqrt(sum_for_sd[i]/(support[i]-1));
    }



    for(int ii=0;ii<nmetrics;++ii)
    {
      //printf("%s\tmean: %f\tsd: %f\tsupport: %i\n", metrics[ii], means[ii], s_ds[ii], support[ii]);
      fprintf(f,"%s\tmean: %.5f\tsd: %.5f\tsupport: %i\n", metrics[ii].c_str(), means[ii], s_ds[ii], support[ii]);
    }
  
  //free(database_);
  free(expanded_query);
  }
);
  fclose(f);
  return 0;
  
}
