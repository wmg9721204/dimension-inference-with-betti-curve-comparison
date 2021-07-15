function remove_zero_cols(M::Array{Float64,2})
    ## find out the indices of the zero columns of Mi
    ind_zero_column = []
    for j = 1:size(M,2)
        if M[:,j]==zeros(size(M,1))
            push!(ind_zero_column, j)
        end
    end
    ind_zero_column
    println("Number of zero columns: $(length(ind_zero_column))")
    return M[:, setdiff(1:size(M,2), ind_zero_column)]
end

function remove_zero_rows(M::Array{Float64,2})
    ## find out the indices of the zero columns of Mi
    ind_zero_rows = []
    for i = 1:size(M,1)
        if M[i,:]==zeros(size(M,2))
            push!(ind_zero_rows, i)
        end
    end
    ind_zero_rows
    println("Number of zero rows: $(length(ind_zero_rows))")
    return M[setdiff(1:size(M,1), ind_zero_rows), :]
end


using Simplicial
#############################################################
## Idea: use normalized speed to construct a filtration

function FiltSparseMatrix(M::Array{Float64,2}, cutoff_size::Int)
    Ord_dic = Dict()
    for i = 1:size(M,1)
        ind_nz_i = findall(M[i,:].!=0)
        Ord_i = Array{Int,1}(sortslices([M[i,ind_nz_i] ind_nz_i], dims = 1)[:,2])
        Ord_dic[i] = Ord_i
    end
    Ord_dic

    Rho = zeros(Int,size(M,1))
    for i = 1:size(M,1)
        Rho[i] = length(Ord_dic[i])
    end

    ind_rem = [i for i = 1:size(M,1) if Rho[i]>=cutoff_size] ## rem = remain
    println("number of rows remained: ", length(ind_rem))
    M_rem = M[ind_rem, :];

    ## Construct the corresponding order matrices. 
    ## Note: the indices are reset to 1:length(ind_rem)
    Ord_dic_rem = Dict()
    for i = 1:length(ind_rem)
        Ord_dic_rem[i] = Ord_dic[ind_rem[i]]
    end
    Ord_dic_rem

    Rho_rem = BigInt[length(Ord_dic_rem[i]) for i = 1:length(ind_rem)];

    ## find out the filtration indices
    ind_filt = BigInt[]
    ind_top = lcm(Rho_rem)

    for i = 1:length(ind_rem)
        temp = Array{BigInt, 1}(Array{BigInt, 1}(collect(1:Rho_rem[i])).*ind_top./Rho_rem[i])
        ind_filt = BigInt[ind_filt; temp]
    end
    ind_filt = sort(collect(Set(ind_filt)));

    ## To efficiently compute the filtration, 
    ## compute the dictionary 
    ## that maps each filtration index to the matrix position
    ind_filt_2_pos = Dict{BigInt,Array{Tuple{Int64,Int64}, 1}}(ind=>[] for ind in ind_filt)
    for i = 1:length(ind_rem)
        diff_i = BigInt(ind_top/Rho_rem[i])
        ##println(diff_i)
        ##println(Rho_rem[i])
        for l = 1:Rho_rem[i]
            ##println("l:", l)
            ##println("diff_i*l: ", diff_i*l)
            ##println("diff_i*l: ", BigInt(diff_i)*l)
            push!(ind_filt_2_pos[diff_i*l], (i,Ord_dic_rem[i][l]))
        end
    end
    ind_filt_2_pos;

    ## Construct the filtration
    Birth = BigInt[]
    Faces = CodeWord[]
    column_2_face = Dict(j=>CodeWord() for j = 1:size(M_rem,2))
    for t = 1:length(ind_filt)
        ind_t = ind_filt[t]
        for pair in ind_filt_2_pos[ind_t]
            push!(Birth, ind_filt[t])
            push!(column_2_face[pair[2]], pair[1])
            push!(Faces, copy(column_2_face[pair[2]]))
        end
    end
    
    ## This block is run in case 
    ## the indices in Birth is too big to be handled by FiltrationOfSimplicialComplexes
    true_birth_dic = Dict{Int, BigInt}(1=>Birth[1])
    new_birth = [1]
    i = 1
    for j = 2:length(Birth)
        if Birth[j]!=Birth[j-1]
            i+=1
            true_birth_dic[i] = Birth[j]
        end
        push!(new_birth, i)
    end
    
    FS = FiltrationOfSimplicialComplexes(Faces, new_birth);
    return (FS, ind_rem, true_birth_dic)
end

function FiltSparseMatrix2(M::Array{Float64,2}, cutoff_size::Int)
    Ord_dic = Dict()
    for i = 1:size(M,1)
        ind_nz_i = findall(M[i,:].!=0)
        Ord_i = Array{Int,1}(sortslices([M[i,ind_nz_i] ind_nz_i], dims = 1)[:,2])
        Ord_dic[i] = Ord_i
    end
    Ord_dic

    Rho = zeros(Int,size(M,1))
    for i = 1:size(M,1)
        Rho[i] = length(Ord_dic[i])
    end

    ind_rem = [i for i = 1:size(M,1) if Rho[i]>=cutoff_size] ## rem = remain
    println("number of rows remained: ", length(ind_rem))
    M_rem = M[ind_rem, :];

    ## Construct the corresponding order matrices. 
    ## Note: the indices are reset to 1:length(ind_rem)
    Ord_dic_rem = Dict()
    for i = 1:length(ind_rem)
        Ord_dic_rem[i] = Ord_dic[ind_rem[i]]
    end
    Ord_dic_rem

    Rho_rem = BigInt[length(Ord_dic_rem[i]) for i = 1:length(ind_rem)];

    ## find out the filtration indices
    ind_filt = BigInt[]
    ind_top = lcm(Rho_rem)

    for i = 1:length(ind_rem)
        temp = Array{BigInt, 1}(Array{BigInt, 1}(collect(1:Rho_rem[i])).*ind_top./Rho_rem[i])
        ind_filt = BigInt[ind_filt; temp]
    end
    ind_filt = sort(collect(Set(ind_filt)));

    ## To efficiently compute the filtration, 
    ## compute the dictionary 
    ## that maps each filtration index to the matrix position
    ind_filt_2_pos = Dict{BigInt,Array{Tuple{Int64,Int64}, 1}}(ind=>[] for ind in ind_filt)
    for i = 1:length(ind_rem)
        diff_i = BigInt(ind_top/Rho_rem[i])
        ##println(diff_i)
        ##println(Rho_rem[i])
        for l = 1:Rho_rem[i]
            ##println("l:", l)
            ##println("diff_i*l: ", diff_i*l)
            ##println("diff_i*l: ", BigInt(diff_i)*l)
            push!(ind_filt_2_pos[diff_i*l], (i,Ord_dic_rem[i][l]))
        end
    end
    ind_filt_2_pos;

    ## Construct the filtration
    Birth = BigInt[]
    Faces = CodeWord[]
    column_2_face = Dict(j=>CodeWord() for j = 1:size(M_rem,2))
    for t = 1:length(ind_filt)
        ind_t = ind_filt[t]
        for pair in ind_filt_2_pos[ind_t]
            push!(Birth, ind_filt[t])
            push!(column_2_face[pair[2]], pair[1])
            push!(Faces, copy(column_2_face[pair[2]]))
        end
    end
    
    ## This block is run in case 
    ## the indices in Birth is too big to be handled by FiltrationOfSimplicialComplexes
    true_birth_dic = Dict{Int, BigInt}(1=>Birth[1])
    new_birth = [1]
    i = 1
    for j = 2:length(Birth)
        if Birth[j]!=Birth[j-1]
            i+=1
            true_birth_dic[i] = Birth[j]
        end
        push!(new_birth, i)
    end
    
    FS = FiltrationOfSimplicialComplexes(Faces, new_birth);
    return (FS, ind_rem, true_birth_dic, ind_top)
end