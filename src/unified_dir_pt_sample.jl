using Random
using LinearAlgebra
"""
dir_sample(dir_type::String, d::Int, m::Int, seed::Int)

"dir_sample" 
samples "m" directions in dimension "d" of type "dir_type"

Inputs:
(1) dir_type: a string, the direction type
(2) d: an integer, the dimension in which the directions are sampled
(3) m: an integer, the number of directions to generate
(4) seed: an integer, used for random sampling

Output:
an m x d matrix, where each row is a direction in R^d
"""
function dir_sample(dir_type::String, d::Int, m::Int, seed::Int,
    theta = pi/4)
    ## set the seed for random number generation
    Random.seed!(seed)
    if dir_type=="pos"
        directions = zeros(0,d)
        num_row = 0
        while num_row<m
            direction = randn(1,d)
            if sum(direction.>0)==d
                directions = [directions; direction]
                num_row+=1
            end
        end
        return directions
    elseif dir_type=="unif"
        directions = randn(m,d)
        return directions
    elseif dir_type=="angle"
        directions = zeros(0,d)
        num_row = 0
        vec_1 = ones(d)/sqrt(d) ## the vector [1,1,...,1]\in\R^d
        num_row = 0
        while num_row<m
            direction = randn(1,d)
            direction = direction/norm(direction)            
            if acos(dot(vec_1, direction'))<theta
                directions = [directions; direction]
                num_row+=1
            end
        end
        return directions            
    end  
end

# ## sample usage code
# dir_type = "pos" ## "unif"
# d = 2
# m = 100
# seed = 1
# A = dir_sample(dir_type::String, d::Int, m::Int, seed::Int)

# using Plots
# pyplot()

# plot()
# for i = 1:size(A,1)
#     plot!([0, A[i,1]], [0,A[i,2]], legend = false)
# end
# plot!()


"""
"UnifRandBall" 
generates a d x N matrix, whose columns represent points, 
randomly uniformly sampled in the unit ball of R^d
"""
function UnifRandBall(d::Int,N::Int)
    Collect = zeros(d,0) ## Set an (d x 0) empty matrix for storing
    i = 0 ## number of collected points
    while i<N
        pt = 2*rand(d).-1 ## (uniformly) generate a random point inside [-1,1]^d
        if norm(pt)<1 ## if the point is inside the unit ball, collect it
            i = i+1
            Collect = [Collect pt]
        end
    end
    return Collect
end

"""
pt_sample(pt_type::String, d::Int, n::Int, seed::Int)

"pt_sample" 
samples "n" points in R^d of type "pt_type"

Inputs:
(1) pt_type: a string, the point type
(2) d: an integer, the dimension in which the directions are sampled
(3) n: an integer, the number of directions to generate
(4) seed: an integer, used for random sampling
(5) r(optional): the reciprocal of the curvature of the de-sitter space, between 0 and 1, exclusive

Output:
a d x n matrix, where each row is a direction in R^d
"""
function pt_sample(pt_type::String, d::Int, n::Int, seed::Int,
        r=0.1::Float64)
    ## set the seed for random number generation
    Random.seed!(seed)
    if pt_type=="uniform"
        pts = UnifRandBall(d,n)
        return pts
    elseif pt_type=="sphere"
        pts = UnifRandBall(d,n)
        for a = 1:n
            pts[:,a] = pts[:,a]/norm(pts[:,a])
        end
        return pts
    elseif pt_type=="de-sitter"
        c = 1/r ## the constant/curvature c in x_1^2+...+x_(d-1)^2-x_d^2 = 1/c^2
        pts = zeros(d,0)
        a = 0
        while a<n
            new_pt = UnifRandBall(d-1,1)
            x_d_a2 = sum(new_pt.^2)-1/c^2
            if x_d_a2>0
                new_pt = [new_pt;(-1)^a*sqrt(x_d_a2)] ## (-1)^a is used to cover both the upper and lower piece
                pts = [pts new_pt]
                a+=1
            end
        end
        return pts
    end
end

# ## sample usage code
# pt_type = "sphere" ## "uniform"
# d = 2
# n = 100
# seed = 1
# X = pt_sample(pt_type::String, d::Int, n::Int, seed::Int)

# using Plots
# pyplot()

# scatter(X[1,:], X[2,:], legend = false)

