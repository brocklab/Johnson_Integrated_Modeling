# Dose_optimization

For the classifier model from scRNAseq data (Python files)
Prior to running scripts, download the data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154932. Local path will need to be indicated in the following scripts.
1. KJ_choose_classifer.py- Loads in the scRNAseq data, assigns pretreatment cells as sensitive or resistant resistant based on change in lineage abundance, uses labeled data set to test out the accuracy of PCA+KNN classifiers and linear SVM classifier using 5-fold CV. Determines a reasonable threshold probability for each classifier using the entire pre-treatment data set. 
2. KJ_PCAKNN.py- classifies the unlabeled cells in the scRNAseq data from all time points using PCA. Makes corresponding figures of projections in PC space, and outputs the phenotypic composition (phi_t), and the variance in PCs- all intended for the supplement. Uses the hyperparameters in KJ_findhyperps.py
3. KJ_SVM.py - classifies the unlabeled cells in the scRNAseq data from all time points using linear SVM. Makes the UMAPS colored by class estimate for each time point and for the whole data set, and outputs the phenotypic composition (phi_t).    


For the dynamic model calibration (MATLAB files)
1. Load_raw_data.m – reads in the excel files from each experiment ( for all cell types and treatment regimens)- generates trajraw.mat
2. Filter_data_231.m loads in trajraw.mat and selects only the 231 cells from this structure, and then selects only those that have only received one treatment. Assigns a color and fits a bi-exponential model to each individual cell well. Records time to reach 2* baseline (critical time).  Generates trajfit231.mat
3. Fit_N_231.m loads in trajfit231.mat and combines wells from the same treatment scenario to get a mean and standard deviation vector for each treatment. Also truncates the data just below the carrying capacity, and adds an effective dose U(t) based on the concentration of dox for each treatment. Then fits the N_t data to the model using the fit_fxn_Greene.m which has 3 ways of calibrating- using normpdf fxn, using weighted least-squares, and using just least-squares. Generated pfitN and trajsumfit231.mat
4. Fit_N_phi_231.m loads in trajsumfit231.mat to get the N(t) data and phi_t_est_pyth_SVM.csv to get the phi(t) data. Fits the joint data using fit_fxn_N_phi.m which does a fit using both the weighted relative error in N and the weighted relative error in phi. It also does the fitting using just the N data and just the weighted relative error in N. Bother outputs of parameter estimates and models are saved.
5. Test_param_predictions loads in trajsumfit231 and the fit parameters from integrated fit and fit on N alone. It then runs everything forward with those parameters- comparing the fit to N(t) for the integrated fit and the N(t) alone fit, and comparing the predictions for each. 
6. Green_sim_mult_doses.m produces the example model predicted figures from Fig.1 comparing pulsed and constant dosing regimens for resistance-preserving (alpha =0) and resistance-inducing (alpha >0) scenarios. 


