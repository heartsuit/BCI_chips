BCI_chips
=========

**Optimization of a space representation for a better SVM classification**

The main file is main_all_sub_test_multi_chan.m

All these optimization methods are undertaken under the **leave-one-out** meaning that for each optimization method, the training set is every signal unless one and we test the signal which is left. When you make a loop over every possible training set, you can then make statistics on the accuracy of the process.  

For more information about the topic please read the PhD thesis of Xavier Artusi.

Once it is opened and launched you will trigger this sucession of function:

wavelet transform with classic db4 wavelet (without any optimization)
---------------------------------------------------------------------
if dwt_without_optim = 1 in the main file, you will run:

* programme_test_descripteur_xval.m  
      * calc_feature (used to have a first set of coefficient computed on db4)
      * svm_learning (used to train the set of coefficient computed on db4)
      * test_class_ovr (used to compute class labeling on the test set)
      * function_acp (let us see in a better way of understanding the axes and sample influence on the results)
                                                              
                                     
wavelet transform with error probability criterion estimation
-------------------------------------------------------------
if optim_pce = 1 in the main file, you will run:                                    
                                     
* programme_test_descripteur_xval.m  
     * create_dic (computes every set of classical wavelet coefficients on a parameterized wavelet over a theta)
     * calc_feature (used to have a first set of coefficient computed on db4, in order to compare the effectiveness of our parameterization)
     * optim_wavelet_multi_cv (gives the best theta and its coefficients according to the pce criterion)
     * svm_learning (used to train the set of coefficient computed on a theta optimal)
     * test_class_ovr (used to compute class labeling on the test set)
     * function_acp (let us see in a better way of understanding the axes and sample influence on the results)

wavelet transform with Fisher criterion estimation
--------------------------------------------------
if optim_fisher = 1 in the main file, you will run:                                    
                                     
* programme_test_descripteur_xval.m  
     * create_dic (computes every set of classical wavelet coefficients on a parameterized wavelet over a theta)
     * calc_feature (used to have a first set of coefficient computed on db4, in order to compare the effectiveness of our parameterization)
     * optim_wavelet_multi_cv (gives the best theta and its coefficients according to the fisher criterion)
     * svm_learning (used to train the set of coefficient computed on a certain theta)
     * test_class_ovr (used to compute class labeling on the test set)
     * function_acp (let us see in a better way of understanding the axes and sample influence on the results)

wavelet packet transform: the best basis adventure
--------------------------------------------------
if optim_basis = 1 in the main file, you will run:                                    
                                     
* programme_test_descripteur_xval.m  
     * set_wp_marg(computes all the mavelet packet marginals)
     * best_basis search (search what will be the best basis among all possible basis from the total set of marginals)
          * calc_fisher_wp_features
          * find_best_basis_fisher
     * svm_learning (used to train the set of marginals chosen thanks to the previous function)
     * test_class_ovr (used to compute class labeling on the test set)
     * function_acp (let us see in a better way of understanding the axes and sample influence on the results)


