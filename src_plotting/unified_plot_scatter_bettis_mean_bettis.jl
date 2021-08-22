## compute estimated area below a (discrete) function curve
function est_area(inds, vals)
    if length(inds)!=length(vals)
        error("The size of indices does not match the size of values!!!")
    end
    area = inds[1]*vals[1]
    for i = 2:length(inds)
        area+=(inds[i]-inds[i-1])*vals[i]
    end
    area
    return area
end

using JLD

"""
"get_bettis_mean_bettis" loads and computes the indices and values of the betti and mean betti
Inputs:
(1) bettis_path: a .jld file, recording the betti results
(2) (Optional) d_up: an integer, deciding the upper limit for the mean betti computation

Outputs:
a tuple (filt_inds, bettis, mean_bettis), where
(1) filt_inds: the filtration indices (in [0,1])
(2) bettis: a matrix; each row is a betti curve; [b0;b1;b2;...]
(3) mean_bettis: a vector, recording mean (b1, b2, b3, ...)
"""

function get_bettis_mean_bettis(bettis_path::String, d_up = 5::Int; d_start = 1)
    D_bettis = load(bettis_path)
    bettis = D_bettis["bettis"]
    println("size of bettis:", size(bettis))
    true_birth_dic = D_bettis["true_birth_dic"]
    denom = D_bettis["denominator"]

    ## compute the (normalized) filtration index
    filt_inds = Float64[true_birth_dic[k]/denom for k = 1:size(bettis,2)]
    mean_bettis = [est_area(filt_inds, bettis[k+1,:]) for k = d_start:d_up-1]
    return (filt_inds, bettis, mean_bettis)
end

using Plots
pyplot()

##################################################################################
## functions for plotting the betti curves
##################################################################################
using ColorSchemes
colors = ColorSchemes.jet
color_inds = [1,3,5,7,9]

trans = 0.2

"""

"""
function plot_betti_curves_real(path_base, name_tail, 
        base_data_name::String, d_up = 5::Int, d_start = 1::Int)
    global colors
    global color_inds
    ## load the bettis result of the real data
    bettis_path = "$(path_base)_$(name_tail).jld"
    filt_inds, bettis_real, mean_bettis = get_bettis_mean_bettis(bettis_path, d_up, d_start = d_start)
    ind_color = 1
    println(size(bettis_real))
    for k = d_start:d_up-1
        bettis_real[k+1,:]
        plot!(filt_inds, bettis_real[k+1,:], 
            label = "$(base_data_name): betti $(k)", linewidth = 3, line = :dash, 
            color = colors[color_inds[ind_color]])
        ind_color+=1
    end
    plot!()
end

"""

"""
function plot_betti_curves_ctrl(base_path::String, name_tail::String, 
        d::Int, seed::Int, 
        dir_type::String, pt_type::String,
        base_data_name::String, 
        d_up = 5::Int, d_start = 1::Int)
    global colors
    global color_inds
    trans = 0.2
    ## load the bettis result of the controlled data
    name_ctrl = "dir_$(dir_type)_pt_$(pt_type)_d_$(d)_seed_$(seed)"
    bettis_path = "$(base_path)_$(name_ctrl)_$(name_tail).jld"
#     try
        filt_inds, bettis, mean_bettis = get_bettis_mean_bettis(bettis_path, d_up)
        ind_color = 1
        for k = d_start:d_up-1
            if seed==1
                plot!(filt_inds, bettis[k+1,:], label = "betti $k", 
                        color = colors[color_inds[ind_color]], alpha = trans, linewidth = 1)
            else
                plot!(filt_inds, bettis[k+1,:], label = "", 
                        color = colors[color_inds[ind_color]], alpha = trans, linewidth = 1)
            end
            ind_color+=1
        end
#     catch
#     end
    plt_title = "$(base_data_name) vs controlled
    d = $d, dir: $(dir_type), pt: $(pt_type)"
    try
        if dir_type[1:5]=="angle"
            angle_ratio = round(parse(Int, dir_type[7:end])/pi)/100
            plt_title = "$(base_data_name) vs controlled
            d = $d, dir: $(angle_ratio) pi, pt: $(pt_type)"            
        end
    catch
    end
    plot!(title = plt_title)
end

"""

"""
function plot_betti_layout(dpi = 300, ylim = (0,100), legend = :topleft)
    Plots.plot!(dpi = dpi, grid = false, ylim = ylim, legend = legend)
end


#####################################################################################
## functions for scattering the mean betti vectors
#####################################################################################
using Plots
pyplot()

using ColorSchemes
colors_scatter = ColorSchemes.jet
color_scatter_inds = [1,2,3,4,5,6,7,8,9]


"""

"""
function scatter_mean_betti_real_2D(path_base, name_tail, 
        label_real, betti_inds = [1,2], 
        real_markersize = 8::Int, d_up = 5::Int, trans = 0.7, 
        color = "black"; d_start = 1)
    ind1 = betti_inds[1]
    ind2 = betti_inds[2]
    bettis_path = "$(path_base)_$(name_tail).jld"
    filt_inds, bettis, mean_bettis = get_bettis_mean_bettis(bettis_path, d_up)
    scatter!([mean_bettis[ind1]], [mean_bettis[ind2]],
        label = label_real, color = color, alpha = trans, 
        markersize = real_markersize, markerstrokewidth=0)
end

"""

"""
function scatter_mean_betti_real_3D(path_base, name_tail, 
        label_real, betti_inds = [1,2,3], 
        real_markersize = 8::Int, d_up = 5::Int, trans = 0.7, 
        color = "black"; d_start = 1)
    ind1 = betti_inds[1]
    ind2 = betti_inds[2]
    ind3 = betti_inds[3]
    bettis_path = "$(path_base)_$(name_tail).jld"
    filt_inds, bettis, mean_bettis = get_bettis_mean_bettis(bettis_path, d_up)
    scatter!([mean_bettis[ind1]],[mean_bettis[ind2]], [mean_bettis[ind3]], 
        label = label_real, color = color, alpha = trans, 
        markersize = real_markersize, markerstrokewidth=0)
end

"""

"""
function scatter_mean_betti_ctrl_2D(base_path::String, name_tail::String, 
        dir_type::String, pt_type::String, num_seed::Int, 
        betti_inds = [1,2], d_up = 5::Int, ctrl_markersize = 1::Int, 
        trans = 0.7::Float64; d_start = 1)
    global colors_scatter
    global color_inds
    ind1 = betti_inds[1]
    ind2 = betti_inds[2]
    for d = 2:10
        mean_bettis_vecs = zeros(d_up-1,0)
        for seed = 1:num_seed
            try
                ## load the bettis result of the controlled data            
                name_ctrl = "dup_$(d_up)_dir_$(dir_type)_pt_$(pt_type)_d_$(d)_seed_$(seed)"            
                bettis_path = "$(base_path)_$(name_ctrl)_$(name_tail).jld"
                filt_inds, bettis, mean_bettis = get_bettis_mean_bettis(bettis_path, d_up)

                mean_bettis_vecs = [mean_bettis_vecs mean_bettis]
            catch
                continue
            end
        end
        scatter!(mean_bettis_vecs[ind1,:], mean_bettis_vecs[ind2,:],
            label = "d: $d", color = colors_scatter[color_scatter_inds[d-1]], 
            alpha = trans,
            markersize = ctrl_markersize, markerstrokewidth=0)
    end
end

"""

"""
function scatter_mean_betti_ctrl_3D(base_path::String, name_tail::String, 
        dir_type::String, pt_type::String, num_seed::Int, 
        betti_inds = [1,2,3], d_up = 5::Int, ctrl_markersize = 1::Int, 
        trans = 0.7::Float64; d_start = 1)
    global colors_scatter
    global color_inds
    ind1 = betti_inds[1]
    ind2 = betti_inds[2]
    ind3 = betti_inds[3]
    for d = 2:10
        mean_bettis_vecs = zeros(d_up-1,0)
        for seed = 1:num_seed
            try
                ## load the bettis result of the controlled data            
                name_ctrl = "dup_$(d_up)_dir_$(dir_type)_pt_$(pt_type)_d_$(d)_seed_$(seed)"            
                bettis_path = "$(base_path)_$(name_ctrl)_$(name_tail).jld"
                filt_inds, bettis, mean_bettis = get_bettis_mean_bettis(bettis_path, d_up)

                mean_bettis_vecs = [mean_bettis_vecs mean_bettis]
            catch
                continue
            end
        end
        scatter!(mean_bettis_vecs[ind1,:],mean_bettis_vecs[ind2,:], mean_bettis_vecs[ind3,:],
            label = "d: $d", color = colors_scatter[color_scatter_inds[d-1]], 
            alpha = trans,
            markersize = ctrl_markersize, markerstrokewidth=0)
    end
end

"""

"""
function scatter_mean_betti_layout(dpi = 300, legend = :topleft)
    plt_title = "fly vs controlled
    dir: $(dir_type), pt: $(pt_type)"
    
    Plots.plot!(title = plt_title,
        dpi = dpi, grid = false, legend = legend) 
end


################################################################################################
##
################################################################################################
using MultivariateStats

"""

"""
function scatter_mean_betti_PCA(d_PCA::Int, 
        base_path_real::String, name_tail_real::String,
        base_path_ctrl::String, name_tail_ctrl::String, 
        dir_type::String, pt_type::String, 
        num_seed::Int, base_title::String, label_real::String, 
        d_up = 5::Int, 
        ctrl_markersize = 1::Int, real_markersize = 5::Int, 
        ctrl_trans = 0.7::Float64, real_trans = 0.7::Float64; d_start = 1)
    
    colors = ColorSchemes.jet
    color_inds = [1,3,5,7,9]
    
    ## plotting dimension should be less than or equal to 3
    ## if not, show error information
    if d_PCA>3
        error("PCA dimension for plotting should be less than or equal to 3.")
    end
    
    ## get the betti-vector of real data
    betti_vec_real = get_bettis_mean_bettis("$(base_path_real)_$(name_tail_real).jld", d_up, d_start=d_start)[3]
    ctrl_d_to_betti_pts = Dict()
    
    ## get the betti-vectors of controlled data
    Pts_ctrl = zeros(d_up-d_start, 0)
    for d = 2:10
        Pts_ctrl_d = zeros(d_up-d_start, 0)
        for seed = 1:num_seed
            name_ctrl = "dir_$(dir_type)_pt_$(pt_type)_d_$(d)_seed_$(seed)"
            bettis_path = "$(base_path_ctrl)_$(name_ctrl)_$(name_tail_ctrl).jld"

            betti_vec_ctrl = get_bettis_mean_bettis(bettis_path, d_up, d_start=d_start)[3]
            Pts_ctrl_d = [Pts_ctrl_d betti_vec_ctrl]
        end
        ctrl_d_to_betti_pts[(dir_type, pt_type, d)] = Pts_ctrl_d
        Pts_ctrl = [Pts_ctrl Pts_ctrl_d]
    end
    Pts_ctrl
    ctrl_d_to_betti_pts
    
    ################################################################################
    # train a PCA model in 2D
    ################################################################################
    M_2D = fit(PCA, [Pts_ctrl betti_vec_real]; maxoutdim=2, pratio = 1.0)
    P_2D = projection(M_2D);

    ################################################################################
    # train a PCA model in 3D
    ################################################################################
    M_3D = fit(PCA, [Pts_ctrl betti_vec_real]; maxoutdim=3)
    P_3D = projection(M_3D);
    
    ########################################################
    ## 2D scatter plot using PCA
    ########################################################
    if d_PCA==2
        ##################################################################
        ## (static) use pyplot backend to generate and save the plot as a .png file
        ##################################################################
        pyplot()
        pt_real_2D = P_2D'*betti_vec_real
        println("P_2D", P_2D)
        scatter(dpi = 300)
        scatter!([pt_real_2D[1]], [pt_real_2D[2]], 
            label = label_real,
            markersize = 15, markerstrokewidth=0, color = :black, alpha = 0.8)
        for d = 2:10
            Pts_ctrl_2D = P_2D'*ctrl_d_to_betti_pts[(dir_type, pt_type, d)]
            scatter!(Pts_ctrl_2D[1,:], Pts_ctrl_2D[2,:], label = "d: $d",
                markerstrokewidth=0, alpha = 0.8, markersize = 5, color = colors[d-1])
        end
        
        try         
            if dir_type[1:5]=="angle"
                angle_ratio = round(parse(Int, dir_type[7:end])/pi)/100
                scatter!(title = "$(base_title) (PCA: 2D)
                    (dir: $(angle_ratio) pi, pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_2D))/100)", 
                    legend=:topright)
            else
                scatter!(title = "$(base_title) (PCA: 2D)
                    (dir: $(dir_type), pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_2D))/100)", 
                    legend=:topright)
            end
        catch
            scatter!(title = "$(base_title) (PCA: 2D)
                    (dir: $(dir_type), pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_2D))/100)", 
                    legend=:topright)
        end

        ## save the plot
        Plots.savefig("scatter_$(label_real)_dir_$(dir_type)_pt_$(pt_type)_PCA_2D.png")
        
        ##################################################################
        ## (interactive) use plotly backend to generate and save the plot as a .png file
        ##################################################################
        plotly()
        pt_real_2D = P_2D'*betti_vec_real
        scatter(dpi = 300)
        scatter!([pt_real_2D[1]], [pt_real_2D[2]], 
            label = label_real,
            markersize = 10, markerstrokewidth=0, color = :black, alpha = 0.8)
        for d = 2:10
            Pts_ctrl_2D = P_2D'*ctrl_d_to_betti_pts[(dir_type, pt_type, d)]
            scatter!(Pts_ctrl_2D[1,:], Pts_ctrl_2D[2,:], label = "d: $d",
                markerstrokewidth=0, alpha = 0.8, markersize = 4, color = colors[d-1])
        end

        try         
            if dir_type[1:5]=="angle"
                angle_ratio = round(parse(Int, dir_type[7:end])/pi)/100
                scatter!(title = "$(base_title) (PCA: 2D)
                    (dir: $(angle_ratio) pi, pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_2D))/100)", 
                    legend=:topright)
            else
                scatter!(title = "$(base_title) (PCA: 2D)
                    (dir: $(dir_type), pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_2D))/100)", 
                    legend=:topright)
            end
        catch
            scatter!(title = "$(base_title) (PCA: 2D)
                    (dir: $(dir_type), pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_2D))/100)", 
                    legend=:topright)
        end

        ## save the plot
        Plots.savefig("scatter_$(label_real)_dir_$(dir_type)_pt_$(pt_type)_PCA_2D.html")   
    end
    
    ########################################################
    ## 3D scatter plot using PCA
    ########################################################
    if d_PCA==3
        plotly()
        pt_real_3D = P_3D'*betti_vec_real
        scatter(dpi = 300)
        scatter!([pt_real_3D[1]], [pt_real_3D[2]], [pt_real_3D[3]], 
            label = label_real,
            markersize = 5, markerstrokewidth=0, color = :black, alpha = 0.8)
        for d = 2:10
            Pts_ctrl_3D = P_3D'*ctrl_d_to_betti_pts[(dir_type, pt_type, d)]
            scatter!(Pts_ctrl_3D[1,:], Pts_ctrl_3D[2,:], Pts_ctrl_3D[3,:], 
                label = "d: $d",
                markerstrokewidth=0, alpha = 0.8, markersize = 2, color = colors[d-1])
        end
               
        try         
            if dir_type[1:5]=="angle"
                angle_ratio = round(parse(Int, dir_type[7:end])/pi)/100
                scatter!(title = "$(base_title) (PCA: 3D)
                    (dir: $(angle_ratio) pi, pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_3D))/100)", 
                    legend=:topright)
            else
                scatter!(title = "$(base_title) (PCA: 3D)
                    (dir: $(dir_type), pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_3D))/100)", 
                    legend=:topright)
            end
        catch
            scatter!(title = "$(base_title) (PCA: 3D)
                    (dir: $(dir_type), pt: $(pt_type))
                    principal ratio: $(round(100*principalratio(M_3D))/100)", 
                    legend=:topright)
        end

        ## save the plot
        Plots.savefig("scatter_$(label_real)_dir_$(dir_type)_pt_$(pt_type)_PCA_3D.html")
    end
    return (ctrl_d_to_betti_pts, betti_vec_real)
end



