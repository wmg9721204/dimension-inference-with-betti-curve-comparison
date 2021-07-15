using Distributed 
addprocs(30);
#####################################################
## data source path
@eval @everywhere data_path = "your_data_path.jld"

#############################################
## set the dimension upper limit for skeleton computation
@eval @everywhere d_up = 5
###############################################################################
## set the cutoff index ratio (i.e. the upper limit for the normalized indices)
@eval @everywhere cutoff_ind_ratio = 0.50
###############################################################################
## set the cutoff_size; i.e. ## if a row has less than "cutoff_size" nonzero entries, kill it
@eval @everywhere cutoff_size = 1
###############################################################################

@eval @everywhere include("src/unified_real_ctrl_fvecs_bettis_PIs_new.jl") ## load the data analysis utils

@everywhere function ctrl(dir_type::String, pt_type::String, d::Int, seed::Int, 
        theta = pi/4, r=0.2)
    global data_path
    global d_up
    global cutoff_ind_ratio
    global cutoff_size
    ctrl_data_fvec_bettis_PIs(data_path::String, d_up::Int,
        cutoff_ind_ratio::Float64, cutoff_size::Int,
        dir_type::String, pt_type::String, seed::Int, d::Int, theta::Float64, r::Float64)
end

@everywhere function to_parallelize_pos_uniform((d, seednum))
    ctrl("pos", "uniform", d, seednum)
    return
end

@everywhere function to_parallelize_pos_sphere((d, seednum))
    ctrl("pos", "sphere", d, seednum)
    return
end

@everywhere function to_parallelize_unif_uniform((d, seednum))
    ctrl("unif", "uniform", d, seednum)
    return
end

@everywhere function to_parallelize_unif_sphere((d, seednum))
    ctrl("unif", "sphere", d, seednum)
    return
end

@everywhere function to_parallelize_unif_de_sitter((d, seednum))
    ctrl("unif", "de-sitter", d, seednum, pi/4, 0.2)
    return
end

@everywhere function to_parallelize_angle_uniform((theta, d, seednum))
    ctrl("angle", "uniform", d, seednum, theta, 0.2)
    return
end

# ## hypers for angle controls
num_angles = 10
num_seeds = 100

Angles = collect(range(pi/4,pi,length=num_angles))
dims = collect(2:10)

hypers_all = reshape([(theta, d, seednum) for seednum = 1:num_seeds, 
            d in dims, theta in Angles],
            length(dims)*num_angles*num_seeds)

# #################################################################################
# ## (new). uniform points with angle-ctrl directions
# ################################################################################
# # run the processes in parallel
pmap(to_parallelize_angle_uniform, hypers_all)

num_seeds = 10
dims = collect(2:10)
seeds = collect(1:num_seeds)
## all possible hyperparameters
hypers_all = reshape([(d, seed) for seed in seeds, d in dims], length(dims)*num_seeds)

# #################################################################################
# ## 1. uniform points with positive directions
# ################################################################################
# # run the processes in parallel
pmap(to_parallelize_pos_uniform, hypers_all)

# #################################################################################
# ## 2. spherical points with positive directions
# #################################################################################
# # run the processes in parallel
# pmap(to_parallelize_pos_sphere, hypers_all)

# #################################################################################
# ## 3. uniform points with uniform directions
# #################################################################################
# #run the processes in parallel
# pmap(to_parallelize_unif_uniform, hypers_all)

# #################################################################################
# ## 4. shperical points with uniform directions
# #################################################################################
# #run the processes in parallel
# pmap(to_parallelize_unif_sphere, hypers_all)

# #################################################################################
# ## 5. de-sitter points with uniform directions
# #################################################################################
# #run the processes in parallel
# pmap(to_parallelize_unif_de_sitter, hypers_all)

## remove the workers
rmprocs(workers())
