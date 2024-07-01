# Unit tests for test_sharp_null function to see if it runs under reasonable configurations of parameters


test_that("If test_sharp_null runs under reasonable configurations of the parameters 
          using our application (for each possible type of test).", {
            
  testdf = MedBounds::kerwin_data;
  
  methods = c("ARP", "CS", "FSST", "CR");
  # bootsize = c(100, 1000, 5000);
  cluster = c(NULL, "Age", "Gender");
  weightmats = c("diag", "avar", "identity"); # Weight matrix used to implement FSST.
  numbins = as.integer(seq.int(from = 2, to = length(unique(testdf[["EL_EGRA_PCA_Index"]])), length.out = 3));
  # alphas = c(0.01, 0.05, 0.1);
  rearranges = c(T, F);
  fix_n1s = c(T, F);
  lambdas = c("dd", "ndd");
  use_ncs = c(T, F);
  analytical_vars = c(T, F); #only available if method is ARP or CS
  defiers_shares = c(0, 0.1, 0.4);
  new_dof_CSs = c(T, F)
  
  # See if the function works with every method
  
  # See if the function works for specific methods with specifications
  for(method in methods){
    
    test_sharp_null_result <- 
      MedBounds::test_sharp_null(df = testdf,
                                 d = "treated",
                                 m = "primarily_leblango",
                                 y = "EL_EGRA_PCA_Index",
                                 method = method,
                                 cluster = cluster)
    
    expect_error(test_sharp_null_result, regexp = NA)
    
    if(method == "FSST"){
      for(weightmat in weightmats){
        test_sharp_null_result <- 
          MedBounds::test_sharp_null(df = testdf,
                                     d = "treated",
                                     m = "primarily_leblango",
                                     y = "EL_EGRA_PCA_Index",
                                     method = method,
                                     # test for different weighting matrix
                                     weight.matrix = weightmat) 
        
        expect_error(test_sharp_null_result, regexp = NA)
        
      }
    }else if(method %in% c("ARP", "CS")){
      test_sharp_null_result <- 
        MedBounds::test_sharp_null(df = testdf,
                                   d = "treated",
                                   m = "primarily_leblango",
                                   y = "EL_EGRA_PCA_Index",
                                   method = method,
                                   # test for analytic variance set to TRUE
                                   analytic_variance = TRUE) 
      
      expect_error(test_sharp_null_result, regexp = NA)
    }else{
      test_sharp_null_result <- 
        MedBounds::test_sharp_null(df = testdf,
                                   d = "treated",
                                   m = "primarily_leblango",
                                   y = "EL_EGRA_PCA_Index",
                                   method = method)
      
      expect_error(test_sharp_null_result, regexp = NA)
    }
  }
  
  # See if the cluster argument works
  
  for(method in methods){
    for(clust in cluster){
      test_sharp_null_result <- 
        MedBounds::test_sharp_null(df = testdf,
                                   d = "treated",
                                   m = "primarily_leblango",
                                   y = "EL_EGRA_PCA_Index",
                                   method = method,
                                   cluster = clust)
      
      expect_error(test_sharp_null_result, regexp = NA)
    }
  }
  
  # Try a number of (64) argument combinations: rearrange, fix_n1, lambda, and use_nc
  for(method in methods){
    for(rearrange in rearranges){
      for(fix_n1 in fix_n1s){
        for(lambda in lambdas){
          for(use_nc in use_ncs){
            
            test_sharp_null_result <- 
              MedBounds::test_sharp_null(df = testdf,
                                         d = "treated",
                                         m = "primarily_leblango",
                                         y = "EL_EGRA_PCA_Index",
                                         method = method,
                                         rearrange = rearrange,
                                         fix_n1 = fix_n1,
                                         lambda = lambda,
                                         use_nc = use_nc)
            
            expect_error(test_sharp_null_result, regexp = NA)
          }
        }
      }
    }
  }
  
  # Try another number of (64) argument combinations: num_Ybins, defiers_shares, new_dof_CSs
  for(method in methods){
    for(numbin in numbins){
      for(def in defiers_shares){
        for(new in new_dof_CSs){
          
            test_sharp_null_result <- 
              MedBounds::test_sharp_null(df = testdf,
                                         d = "treated",
                                         m = "primarily_leblango",
                                         y = "EL_EGRA_PCA_Index",
                                         method = method,
                                         num_Ybins = numbin,
                                         defiers_share = def,
                                         new_dof_CS = new)
            
            expect_error(test_sharp_null_result, regexp = NA)
          
        }
      }
    }
  }
})

# n=1
# for(method in methods){
#   for(boot in bootsize){
#     for(clust in cluster){
#       for(alpha in alphas){
#         for(numbin in numbins){
#           for(rearrange in rearranges){
#             for(fix_n1 in fix_n1s){
#               for(lambda in lambdas){
#                 for(use_nc in use_ncs){
#                   for(def in defiers_shares){
#                     for(dof in new_dof_CSs){
# 
#                       print(n)
#                       n=n+1
# 
#                     }
#                   }
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#   }
# }

# total of 46080 combos 




