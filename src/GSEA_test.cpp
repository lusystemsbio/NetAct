#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix GSEA_test(Rcpp::List GSDB, Rcpp::NumericVector stats_vector, int nperm) {
  try
  {
    
    // total number of transcription factors in the database
    int nTF = GSDB.size();
    //number of genes
    int nGenes = stats_vector.size();
    
    // Reformat the database from gene names to arbitrary gene numbers
    std::vector< std::vector<int> > database;
    
    std::vector<int> targetCount(nTF); //Number of targets of each TF
    for (int i=0; i<nTF; ++i) {
      database.push_back(GSDB[i]);
      targetCount[i] = database[i].size();
    }
    
    // Gene index of each gene--used to permute the genes
    std::vector<int> geneIndex(nGenes);
    for(int i =0; i<nGenes; ++i)
    {     geneIndex[i]=i;}
    
    //Output matrix
    Rcpp::NumericMatrix ES_cal_mat(nTF, nperm);
    
    // weighted sum vector
    std::vector<double> stat_vector_sum(nGenes);
    
    for(int permCount=0;permCount<nperm; ++permCount)
    {
      if(permCount % (nperm/10)==0) Rcout<<"=====";
      
      
      // Reorder the gene index
  //    std::random_shuffle ( geneIndex.begin(), geneIndex.end() );
      
      for ( int i=0; i<database.size(); ++i )
      {
        for ( int j=0;j<database[i].size(); ++j)
        {
          database[i][j] = geneIndex[database[i][j]];
        }
      }
      
      double miss[nTF];
      for ( int i=0; i<nTF; ++i )
      {
        miss[i] = -1.0/(double(nGenes - targetCount[i]));
        
      }
      
      for ( int i=0; i<nTF; ++i )
      {
        
        for(int temp1=0;temp1<nGenes;++temp1) {stat_vector_sum[temp1] = miss[i];}
        
        double hit = 1.0/double(targetCount[i]);
        double Nr = 0;
        for ( int j=0;j<targetCount[i]; ++j)
        {
          Nr += stats_vector[database[i][j]]*hit;
        }
        
        for ( int j=0;j<targetCount[i]; ++j)
        {
          stat_vector_sum[database[i][j]] = stats_vector[database[i][j]]*hit/Nr;
        }
        
        for(int temp1=1;temp1<nGenes;++temp1) {stat_vector_sum[temp1] += stat_vector_sum[(temp1 - 1)];}
        double max = *std::max_element(stat_vector_sum.begin(), stat_vector_sum.end());
        double min = *std::min_element(stat_vector_sum.begin(), stat_vector_sum.end());
        
        ES_cal_mat(i, permCount) = (max< (-min))?min:max;
        
      }
      
    }
    
    return ES_cal_mat;  
    
  }
  catch (Rcpp::internal::InterruptedException& e)
  {
    Rcout << "Program interrupted!" << std::endl;
  }
  return 0;
  
}

