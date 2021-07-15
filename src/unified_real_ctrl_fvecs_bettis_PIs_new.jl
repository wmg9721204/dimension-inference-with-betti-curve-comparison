###############################################
## load the required function: FiltSparseMatrix
###############################################
#include("FiltSparseMatrix.jl")
#include("utils-de-error.jl");
include("FiltSparseMatrix.jl")
include("utils-de-error.jl");

using JLD
using DelimitedFiles ## for using "writedlm"

###########################################
## For Real-Data
###########################################
"""
"real_data_fvec_bettis_PIs" computes 
the fvector, betti curves, and persistence intervals of a truncated filtration

Inputs:
(1) data_path: data source path
(2) d_up: dimension upper limit for persistent homology computation
(3) cutoff_size: row removal cutoff threshold
(4) cutoff_ind_ratio: filtration index cutoff

Output: 
Nothing; but saves the results (and time peeks) in folders
"""
function real_data_fvec_bettis_PIs(data_path::String, d_up::Int, 
        cutoff_size::Int, cutoff_ind_ratio::Float64)
    ## load the source data as a matrix M
    Data = load(data_path)
    M = Data["response"]

    ## the directory name of data source
    dir_data = join(split(data_path, "/")[1:end-1], "/")
    ## the directory name to save results; make directory if necessary
    dir_results = dir_data*"_results_real_0$(Int(cutoff_ind_ratio*100))"
    if !(isdir(dir_results))
        mkdir(dir_results)
    end

    ## the directory name to save progress peek; make directory if necessary
    dir_prog_peek = dir_data*"_prog_peek"
    if !(isdir(dir_prog_peek))
        mkdir(dir_prog_peek)
    end

    base_name = split(data_path, "/")[end][1:end-4] ## get the base name from data_path
    ctrl_name = "dup_$(d_up)_cutoff_$(cutoff_size)_0$(Int(cutoff_ind_ratio*100))"
    Progress_peek_name = dir_prog_peek*"/prog_peek_$(base_name)_$(ctrl_name)"

    ## apply the doubling trick
    M = Array{Float64,2}([M;-M])

    Time_record = [] ## for recording progress peek
    #########################
    ## compute the filtration
    #########################
    println("Working on filtration computation")
    push!(Time_record, "Working on filtration computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    t_filt = @elapsed FS, ind_rem, true_birth_dic, denom = FiltSparseMatrix2(M, cutoff_size);    
    push!(Time_record, t_filt)
    writedlm(Progress_peek_name*".txt", Time_record)

    #############################
    ## compute normalized indices
    #############################
    ind_norm = [true_birth_dic[k]/denom for k = 1:FS.birth[end]]
    ## compute the cutoff (un-normalized) index
    cutoff_ind = findlast(ind_norm.<=cutoff_ind_ratio)
    println("LOOK HERE (FOR DEBUG)")
    println(ind_norm[1:cutoff_ind])
    ## compute the index-truncated filtration
    birth_cut = FS.birth[FS.birth.<=cutoff_ind]
    faces_cut = FS.faces[FS.birth.<=cutoff_ind]
    FS_cut = FiltrationOfSimplicialComplexes(faces_cut, birth_cut);
    
    ####################################
    ## compute the skeleton and fvectors
    ####################################
    println("Working on Skeleton/f_vector computation")
    push!(Time_record, "Working on Skeleton/f_vector computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    t_sk_fvec = @elapsed FS_Skeleton, f_vector = Skeleton_and_fvector(FS_cut, d_up)
    push!(Time_record, t_sk_fvec)
    writedlm(Progress_peek_name*".txt", Time_record)
    ## set the file name to save fvectors
    fvec_name = dir_results* "/$(base_name)_f_vectors_$(ctrl_name)"
    save(fvec_name*".jld", 
        "f_vector", f_vector, "true_birth_dic", true_birth_dic, "denominator", denom)

    ####################################
    ## compute the persistence intervals
    ####################################
    println("Working on PI computation")
    push!(Time_record, "Working on PI computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    t_PI = @elapsed PI = PersistenceIntervals(FS_Skeleton, Inf, "$(base_name)_$(ctrl_name)")
    push!(Time_record, t_PI)
    writedlm(Progress_peek_name*".txt", Time_record)
    ## save computed persistence intervals
    println("Saving computed PI")
    push!(Time_record, "Saving computed PI")
    ## set the file name to save PIs
    PIs_name = dir_results*"/$(base_name)_PIs_$(ctrl_name)"
    t_save_PI = @elapsed save(PIs_name*".jld", 
        "PI", PI, "true_birth_dic", true_birth_dic, "denominator", denom)
    push!(Time_record, t_save_PI)
    writedlm(Progress_peek_name*".txt", Time_record)

    ###########################
    ## compute the betti curves
    ###########################
    println("Working on betti curves computation")
    push!(Time_record, "Working on betti curves computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    PI_aug = [[zeros(2,2)]; PI]
    t_bettis = @elapsed bettis = Intervals2Bettis(PI_aug, FS_cut.birth[end])
    push!(Time_record, t_bettis)
    writedlm(Progress_peek_name*".txt", Time_record)
    ## set the file name to save betti curves
    bettis_name = dir_results*"/$(base_name)_bettis_$(ctrl_name)"
    save(bettis_name*".jld", 
        "bettis", bettis, "true_birth_dic", true_birth_dic, "denominator", denom)
end

# ## Usage demonstration
# data_path = pwd()*"/flyData/fly_020_data.jld" ## data source path

# d_up = 5 ## dimension upper limit for persistent homology computation
# cutoff_size = 60 ## row removal cutoff threshold
# cutoff_ind_ratio = 0.65 ## filtration index cutoff

# real_data_fvec_bettis_PIs(data_path, d_up, cutoff_size, cutoff_ind_ratio)


###########################################
## For Controlled-Data
###########################################

include("unified_dir_pt_sample.jl") ## loading functions for generating random directions and points
"""
"ctrl_data_fvec_bettis_PIs" computes 
the fvector, betti curves, and persistence intervals of a truncated filtration

Inputs:
(1) data_path: data source path
(2) d_up: dimension upper limit for persistent homology computation
(3) cutoff_size: row removal cutoff threshold
(4) cutoff_ind_ratio: filtration index cutoff
(5) dir_type: a string, the direction type (in ["pos", "unif"])
(6) pt_type: a string, the point type (in ["uniform", "sphere"])
(7) seed: an integer, for random number generation
(8) d: an integer, the dimension d where the controlled points and directions are in R^d

Output: 
Nothing; but saves the results (and time peeks) in folders
"""
function ctrl_data_fvec_bettis_PIs(data_path::String, d_up::Int,
        cutoff_ind_ratio::Float64, cutoff_size::Int,
        dir_type::String, pt_type::String, seed::Int, d::Int, 
        theta=pi/4::Float64, r=0.2::Float64)
    ## load the source data as a matrix M
    Data = load(data_path)
    M = Data["response"]

    ## the directory name of data source
    dir_data = join(split(data_path, "/")[1:end-1], "/")
    ## the directory name to save results; make directory if necessary
    dir_results = dir_data*"_results_ctrl_0$(Int(cutoff_ind_ratio*100))"
    if !(isdir(dir_results))
        mkdir(dir_results)
    end

    ## the directory name to save progress peek; make directory if necessary
    dir_prog_peek = dir_data*"_prog_peek"
    if !(isdir(dir_prog_peek))
        mkdir(dir_prog_peek)
    end

    base_name = split(data_path, "/")[end][1:end-4] ## get the base name from data_path
    ctrl_name = "dup_$(d_up)_dir_$(dir_type)_pt_$(pt_type)_d_$(d)_seed_$(seed)_cutoff_$(cutoff_size)_0$(Int(cutoff_ind_ratio*100))"
    if dir_type=="angle"
        name_angle = Int(round(theta*100))
        ctrl_name = "dup_$(d_up)_dir_$(dir_type)_$(name_angle)_pt_$(pt_type)_d_$(d)_seed_$(seed)_cutoff_$(cutoff_size)_0$(Int(cutoff_ind_ratio*100))"
    end
    Progress_peek_name = dir_prog_peek*"/prog_peek_$(base_name)_$(ctrl_name)"

    m, n = size(M)

    ## compute the number of nonzero entries in each row
    Rho = [length(findall(M[i,:].!=0)) for i = 1:m]
    directions = zeros(m,d)
    if dir_type=="angle"
        directions = dir_sample(dir_type::String, d::Int, m::Int, seed::Int, theta::Float64)
    else
        directions = dir_sample(dir_type::String, d::Int, m::Int, seed::Int)
    end
    pts = pt_sample(pt_type::String, d::Int, n::Int, seed::Int, r)

    ## generate the control matrix
    A_ctrl = directions*pts
    println(size(A_ctrl))
    for i = 1:m
        ## compute the column indices in row i with the top Rho_rem[i] values
        ind_keep_i = Array{Int,1}(sortslices([A_ctrl[i,:] collect(1:n)], dims = 1)[end-Rho[i]+1:end, 2])
        ## make the to-keep values all positive
        if ind_keep_i!=[]
            A_ctrl[i, ind_keep_i].-=(minimum(A_ctrl[i,ind_keep_i])-1)
        end
        ## make the to-remove values all zero
        ind_rm_i = setdiff(collect(1:n), ind_keep_i)
        #println("remove_$(i): ", ind_rm_i)
        A_ctrl[i, ind_rm_i] = zeros(n-Rho[i])
    end
    A_ctrl

    ## apply the doubling trick
    M_ctrl = [A_ctrl; -A_ctrl];

    Time_record = [] ## for progress peeking
    ## compute the filtration
    println("Working on filtration computation")
    push!(Time_record, "Working on filtration computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    t_filt = @elapsed FS, ind_rem, true_birth_dic, denom = FiltSparseMatrix2(M_ctrl, cutoff_size);
    push!(Time_record, t_filt)
    writedlm(Progress_peek_name*".txt", Time_record)

    ## compute normalized indices
    ind_norm = [true_birth_dic[k]/denom for k = 1:FS.birth[end]]
    ## compute the cutoff (un-normalized) index
    cutoff_ind = findlast(ind_norm.<=cutoff_ind_ratio)

    birth_cut = FS.birth[FS.birth.<=cutoff_ind]
    faces_cut = FS.faces[FS.birth.<=cutoff_ind]
    FS_cut = FiltrationOfSimplicialComplexes(faces_cut, birth_cut);

    ## compute the skeleton and fvectors
    println("Working on Skeleton/f_vector computation")
    push!(Time_record, "Working on Skeleton/f_vector computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    t_sk_fvec = @elapsed FS_Skeleton, f_vector = Skeleton_and_fvector(FS_cut, d_up)
    push!(Time_record, t_sk_fvec)
    writedlm(Progress_peek_name*".txt", Time_record)
    ## set the file name to save fvectors
    fvec_name = dir_results*"/$(base_name)_f_vectors_$(ctrl_name)"
    save(fvec_name*".jld", 
        "f_vector", f_vector, "true_birth_dic", true_birth_dic, "denominator", denom)

    ## compute the persistence intervals
    println("Working on PI computation")
    push!(Time_record, "Working on PI computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    t_PI = @elapsed PI = PersistenceIntervals(FS_Skeleton, Inf, "$(base_name)_$(ctrl_name)")
    push!(Time_record, t_PI)
    writedlm(Progress_peek_name*".txt", Time_record)
    ## save computed persistence intervals
    println("Saving computed PI")
    push!(Time_record, "Saving computed PI")
    ## set the file name to save PIs
    PIs_name = dir_results*"/$(base_name)_PIs_$(ctrl_name)"
    t_save_PI = @elapsed save(PIs_name*".jld", 
        "PI", PI, "true_birth_dic", true_birth_dic, "denominator", denom)
    push!(Time_record, t_save_PI)
    writedlm(Progress_peek_name*".txt", Time_record)

    ## compute the betti curves
    println("Working on betti curves computation")
    push!(Time_record, "Working on betti curves computation")
    writedlm(Progress_peek_name*".txt", Time_record)
    PI_aug = [[zeros(2,2)]; PI]
    t_bettis = @elapsed bettis = Intervals2Bettis(PI_aug, FS_cut.birth[end])
    push!(Time_record, t_bettis)
    writedlm(Progress_peek_name*".txt", Time_record)
    ## set the file name to save betti curves
    bettis_name = dir_results*"/$(base_name)_bettis_$(ctrl_name)"
    save(bettis_name*".jld", 
        "bettis", bettis, "true_birth_dic", true_birth_dic, "denominator", denom)
end


# ## Usage sample code for controlled data

# include("unified_real_ctrl_fvecs_bettis_PIs.jl");
# ## parameters
# data_path = pwd()*"/flyData/fly_020_data.jld" ## data source path

# d_up = 5 ## dimension upper limit for persistent homology computation
# cutoff_ind_ratio = 0.3 ## filtration index cutoff
# cutoff_size = 15 ## row removal cutoff threshold
# dir_type = "pos" ## direction control, choices in ["pos", "unif"]
# pt_type = "uniform" ## point control, choices in ["uniform", "sphere", "hyperbolic"]
# seed = 1
# d = 2

# ctrl_data_fvec_bettis_PIs(data_path::String, d_up::Int,
#         cutoff_ind_ratio::Float64, cutoff_size::Int,
#         dir_type::String, pt_type::String, seed::Int, d::Int)