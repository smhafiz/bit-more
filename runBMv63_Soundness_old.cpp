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
  int ntrials = 4, ntrials_pruned = 2, nmetrics = 4, confidence = 2;
  std::string metrics[] = {"All query generation (us):\t\t", "Query expansion per server (ms):\t", "Response generation per server (us):", "Record reconstruction (ns):\t\t"};
  struct timeval tvalBefore, tvalAfter;
  FILE *f = fopen("soundness_pir_Jan21.txt", "w");
  int err;

  AES_KEY aeskey;
  AES_set_encrypt_key(_mm_set_epi64x(597349, 121379), &aeskey);

  constexpr std::array<uint8_t, 1> Lognitems = {22};
  
  for_constexpr<for_bounds<0, 1>, for_bounds<0, 8>>(
      [&](auto ir, auto is) 
      
  {
    constexpr uint8_t lognitems = std::get<ir.value>(Lognitems);
    constexpr size_t nitems = (1ULL << lognitems);
    constexpr size_t nbytes_per_row = 512;
    //printf("***********************************\nEach row: %f KiB, Database: %f GiB.\n", (float)nbytes_per_row/(1024), (float)(nitems*nbytes_per_row)/(1024*1024*1024));
    fprintf(f,"***********************************\nEach row: %.3f KiB, Database: %.5f GiB, In memory: %.5f GiB.\n", (float)nbytes_per_row/(1024), (float)(nitems*nbytes_per_row)/(1024*1024*1024),(float)(alloc_size)/(1024*1024*1024));    
    constexpr size_t item = 101;

    void * database_;
    err = posix_memalign((void**)&database_, sizeof(__m256i), nitems * nbytes_per_row * 1.5);
    if(err) perror("Error in memalign for database_");
    arc4random_buf(database_, nitems * nbytes_per_row * 1.5);

    uint8_t * expanded_query;
    err = posix_memalign((void**)&expanded_query, sizeof(__m256i), nitems * sizeof(uint8_t));
    if(err) perror("Error in memalign for expanded_query");

    constexpr std::array<uint16_t, 8> Soundness = {32, 48, 64, 96, 112, 128, 192, 256};
    //constexpr std::array<uint8_t, 4> Nservers = {3, 6, 10, 15};

    constexpr uint8_t nservers = 5;
    constexpr size_t nwords_per_row = nservers-1;
    constexpr size_t nbytes_per_word = std::ceil(nbytes_per_row / static_cast<double>(nwords_per_row));
    constexpr size_t soundness = std::get<is.value>(Soundness);
    constexpr size_t nkeys = std::ceil(std::log2(nservers)) + soundness;
    
    double sum[nmetrics];
    for (int i = 0; i < nmetrics; ++i) { sum[i] = 0.0;  }
      
    double time_elapses[nmetrics][ntrials];
    //printf("\t\trows: 2^%u, ell: %u\n", lognitems, nservers);
    fprintf(f,"\t\trows: 2^%u, ell: %u, soundness: %zu\n", lognitems, nservers, soundness);
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

      typedef bitmore::record<nbytes_per_word, nwords_per_row> record;
      typedef bitmore::query<nitems, nkeys, nservers> query;

      record * database = reinterpret_cast<record *>(database_);
    
      //Client: query generation
      gettimeofday(&tvalBefore, NULL);
      query q(aeskey, item);
      auto Q = q.get_request(nwords_per_row);
      gettimeofday(&tvalAfter, NULL);
      time_elapses[0][trials] = getTimeElapsedMicroSec(tvalAfter, tvalBefore);
      //fprintf(f, "%lf\n",time_elapses[trials]);
      //sum[0]+=time_elapses[0][trials];


      //Query expansion at s+1 th server
      gettimeofday(&tvalBefore, NULL);
      bitmore::evalfull<nkeys,nitems,nservers>(aeskey, Q, s, t, expanded_query);
      gettimeofday(&tvalAfter, NULL);
      time_elapses[1][trials] = getTimeElapsedMilliSec(tvalAfter, tvalBefore);
      //fprintf(f, "%lf\n",time_elapses[trials]);
      //sum[1]+=time_elapses[1][trials];

      //response word from s+1 th server
      gettimeofday(&tvalBefore, NULL);
      bitmore::word<nbytes_per_word> qtemplate = { 0 };
      for (size_t i = 0; i < nitems; ++i)
      {
        if (expanded_query[i] != nwords_per_row) qtemplate ^= database[i][expanded_query[i]];
      }
      gettimeofday(&tvalAfter, NULL);
      time_elapses[2][trials] = getTimeElapsedMicroSec(tvalAfter, tvalBefore);
      //fprintf(f, "%lf\n",time_elapses[trials]);
      //sum[2]+=time_elapses[2][trials];


      //Query expansion at the rest s number of servers, 0...s-1
      // gettimeofday(&tvalBefore, NULL);
      // for (size_t k = 0; k < nservers-1; ++k)
      // {
      //   Q = q.get_request(k);
      //   bitmore::evalfull<nkeys,nitems,nservers>(aeskey, Q, s, t, expanded_query);
      // }
      // gettimeofday(&tvalAfter, NULL);
      // time_elapses[1][trials] = getTimeElapsed(tvalAfter, tvalBefore);
      // //fprintf(f, "%lf\n",time_elapses[trials]);
      // sum[1]+=time_elapses[1][trials];


      //response word generation at the rest s number of servers, 0...s-1
      
      record result = { qtemplate };
      for (size_t k = 0; k < nservers-1; ++k)
      {
        gettimeofday(&tvalBefore, NULL);
        Q = q.get_request(k);
        gettimeofday(&tvalAfter, NULL);
        time_elapses[0][trials] += getTimeElapsedMicroSec(tvalAfter, tvalBefore);

        gettimeofday(&tvalBefore, NULL);
        bitmore::evalfull<nkeys,nitems,nservers>(aeskey, Q, s, t, expanded_query);
        gettimeofday(&tvalAfter, NULL);
        time_elapses[1][trials] += getTimeElapsedMilliSec(tvalAfter, tvalBefore);
        //fprintf(f, "%lf\n",time_elapses[trials]);
        //sum[1]+=time_elapses[1][trials];

        gettimeofday(&tvalBefore, NULL);
        for (size_t i = 0; i < nitems; ++i)
        {
          if (expanded_query[i] != nwords_per_row) result[k] ^= database[i][expanded_query[i]];
        }
        gettimeofday(&tvalAfter, NULL);
        time_elapses[2][trials] += getTimeElapsedMicroSec(tvalAfter, tvalBefore);
        //fprintf(f, "%lf\n",time_elapses[trials]);
        //sum[2]+=time_elapses[2][trials];
      }
      time_elapses[1][trials]/=nservers;
      time_elapses[2][trials]/=nservers;
      sum[0]+=time_elapses[0][trials];
      sum[1]+=time_elapses[1][trials];
      sum[2]+=time_elapses[2][trials];

      //Client: record reconsteuction
      gettimeofday(&tvalBefore, NULL);
      for (size_t i = 0; i < nwords_per_row; ++i)
      {
        result[i] ^= database[item][i];
      }
      gettimeofday(&tvalAfter, NULL);
      time_elapses[3][trials] = getTimeElapsedNanoSec(tvalAfter, tvalBefore);
      //fprintf(f, "%lf\n",time_elapses[trials]);
      sum[3]+=time_elapses[3][trials];
      bool correct_ness = true;
      for (size_t i = 0; i < nwords_per_row; ++i)
      {
        if(result[i].m256[0][0]!=0x00) {correct_ness = false; printf("%d->Incorrect!\n",nservers);break;}
      }
      if (correct_ness) printf("%zu[%d]->Correct protocol!\n",soundness, trials);

      for (size_t i = 0; i < nkeys; ++i)
      {
        free(s[i]);
        free(t[i]);
      }
      free(s);
      free(t);
    }
    double means[nmetrics];
    for (int i = 0; i < nmetrics; ++i) { means[i]  = sum[i]/ntrials;  }
    //printf("m: %f\t",mean_time_elapse);
    double sum_for_sd[nmetrics], s_ds[nmetrics];
    for (int i = 0; i < nmetrics; ++i) 
    {
      for(int tt=0;tt<ntrials;++tt){
        sum_for_sd[i]+=(time_elapses[i][tt]-means[i])*(time_elapses[i][tt]-means[i]);
      }
      s_ds[i] = std::sqrt(sum_for_sd[i]/(ntrials-1));
    }

    // for(int ii=0;ii<nmetrics;++ii)
    // {
    //   std::cout << "sum[" << ii << "] " << sum[ii] << std::endl;
    //   for (int jj = 0; jj < ntrials; ++jj)
    //   {
    //     std::cout << time_elapses[ii][jj] << std::endl;
    //   }
    //   std::cout << "mean: " << means[ii] << "\tsd: " << s_ds[ii] << std::endl;
    // }

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
      // std::cout << "sum[" << ii << "] " << sum[ii] << std::endl;
      // for (int jj = 0; jj < support[ii]; ++jj)
      // {
      //   std::cout << time_elapses_pruned[ii][jj] << std::endl;
      // }
      //std::cout << metrics[ii] << "\tmean: " << means[ii] << "\tsd: " << s_ds[ii] << "\tsupport: " << support[ii] << std::endl;
      //printf("%s\tmean: %f\tsd: %f\tsupport: %i\n", metrics[ii], means[ii], s_ds[ii], support[ii]);
      fprintf(f,"%s\tmean: %.5f\tsd: %.5f\tsupport: %i\n", metrics[ii].c_str(), means[ii], s_ds[ii], support[ii]);
    }
  
  free(database_);
  free(expanded_query);
  }
);
  fclose(f);
  return 0;
  
}