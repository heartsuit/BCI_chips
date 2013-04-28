BCI_chips
=========

optimization of a space representation for a better SVM classification

the main file is main_all_sub_test_multi_chan.m

all these optimization are undertaken under the **leave-one-out** meaning that for each optimization method, the training set is every signal unless one we test the signal which is left. when you make a loop on every possible training set, you can then make statistic on the accuracy of the process.  

once it is opened and launched you will trigger this sucession of function:

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
     * create_dic (computes every set of classical wavelet coefficients on a parameterized wavelet)
     * calc_feature (used to have a first set of coefficient computed on db4, in order to compare the effectiveness of our parameterization)
     * optim_wavelet_multi_cv (gives the best theta and its coefficients)
     * svm_learning (used to train the set of coefficient computed on a theta optimal)
     * test_class_ovr (used to compute class labeling on the test set)
     * function_acp (let us see in a better way of understanding the axes and sample influence on the results)

wavelet transform with Fisher criterion estimation
--------------------------------------------------
if optim_fisher = 1 in the main file, you will run:                                    
                                     
* programme_test_descripteur_xval.m  
     * create_dic (computes every set of classical wavelet coefficients on a parameterized wavelet)
     * calc_feature (used to have a first set of coefficient computed on db4, in order to compare the effectiveness of our parameterization)
     * optim_wavelet_multi_cv (gives the best theta and its coefficients)
     * svm_learning (used to train the set of coefficient computed on a certain theta)
     * test_class_ovr (used to compute class labeling on the test set)
     * function_acp (let us see in a better way of understanding the axes and sample influence on the results)

