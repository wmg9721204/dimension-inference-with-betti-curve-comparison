include("src/unified_real_ctrl_fvecs_bettis_PIs_new.jl") ## load the data analysis utils

## Usage demonstration
data_path = "your_data_path.jld" ## data source path
d_up = 5 ## dimension upper limit for persistent homology computation
cutoff_size = 1 ## row removal cutoff threshold
cutoff_ind_ratio = 0.50 ## filtration index cutoff

real_data_fvec_bettis_PIs(data_path, d_up, cutoff_size, cutoff_ind_ratio)


