###########################################################################
###################		                            		 ####################
###################  		LOAD FILES AND FUNCTIONS  		 ####################
###################		 																 ####################
###########################################################################


# set data
devtools::load_all()
devtools::document()
library(magrittr)
experiment_fname <- "Data/center-out-experiment3D-1all.Rdata"
load(experiment_fname)
opt_lag   <- 5
dead_channels <- which(colSums(UnsortedMat.EncLAG[[opt_lag]])==0)

cat("*** DATA LOADED *** \n")

###########################################################################
###################		                            		 ####################
###################  			BUILD SC AND WF MATRIX  		 ####################
###################		 																 ####################
###########################################################################

lag_used <- 1:13
UnsortedMat.EncLAG <- UnsortedMat.EncLAG[lag_used]
UnsortedMat.DecLAG <- UnsortedMat.DecLAG[lag_used]
Z_transf <- get_transformed_scdata(UnsortedMat.EncLAG, UnsortedMat.DecLAG,
																dead_channels=dead_channels)
Z        <- Z_transf$Z
Ztest    <- Z_transf$Ztest

cat("*** SPIKE COUNT TRANSFORMATION COMPLETED *** \n")

build_wf_mat <- function(wf_mat, dead_channels)
	wf_mat %>%
	lapply(function(x) lapply(x, function(y) y[,-c(dead_channels)]))

WF_Enc    <- build_wf_mat(WF_Enc, dead_channels)
WF_Enc_sq <- build_wf_mat(WF_Enc_sq, dead_channels)
WF_Dec    <- build_wf_mat(WF_Dec, dead_channels)
WF_Dec_sq <- build_wf_mat(WF_Dec_sq, dead_channels)

WF_Enc_tot <- cbind(do.call(cbind, do.call(cbind, WF_Enc)),
									 do.call(cbind, do.call(cbind, WF_Enc_sq)))
WF_Dec_tot <- cbind(do.call(cbind, do.call(cbind, WF_Dec)),
									 do.call(cbind, do.call(cbind, WF_Dec_sq)))

Z_tot <- cbind(Z, abs(WF_Enc_tot))
Ztest_tot <- cbind(Ztest, abs(WF_Dec_tot))

cat("*** SC'S AND WF'S TRANSFORMATION COMPLETED *** \n")





###########################################################################
###################		                            		 ####################
###################  				LOAD EQUATIONS INFO  			 ####################
###################		 																 ####################
###########################################################################

# compute the number of transformations
n_units   <- ncol(UnsortedMat.DecLAG[[opt_lag]]) - length(dead_channels)
n_lags    <- length(UnsortedMat.DecLAG)
n_transf  <- ncol(Z)/(n_units*n_lags)
n_wf      <- 2 #ampl + ampl2
n_feat    <- 4


# obtain electrodes information for only spike counts model
equationinfo_fname <- "Data/tot_model_info.RData"
# get_equations_info(experiment_fname, Z_tot,
#                    n_units, n_lags, n_wf, n_feat, n_transf,
#   								 save_name = equationinfo_fname, vel = velocity.Train)
load(equationinfo_fname)
eqs_sorted_eq <- eqs_info
cat("*** SPIKE COUNTS INFORMATION EXTRACTED *** \n")
