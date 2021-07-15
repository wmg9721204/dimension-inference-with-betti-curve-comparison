using Combinatorics, Simplicial

"""
The following highly specialized and technical function takes an array of integer vectors (faces), and an array of
births (a birth time for each face). This function implementation is extremely ugly due to speed optimization, alas!
-------------------------------------------
This function does _two_ different things:
 1) It removes all redundant faces
 2) For non-redundant faces it changes the births to the minimal possible birth in each equivalence class of faces.
------------------------------------------
This function can be called in the following two different ways:
RemoveDuplicates!(faces,births);
or, alternatively,
indx, p, births = RemoveDuplicates(faces,births);
----------------------------------------
Here,
p = p=sortperm(faces) ;
indx is the index of the redundant faces according to the new ordering of faces
births is the new vector of faces reordered according to p with modified births (and some deleted)
"""
function  RemoveDuplicates(faces2::Array{Array{TheIntegerType,1},1},  births2::Array{Int,1}) # ::Array{Int,1}
N=length(faces2)
if N != length(births2)
    error(" Mismatch of array lengths in the arguments of RemoveDuplicates!")
end
  p=sortperm(faces2); # lexicographically order faces
  faces=faces2[p]; # Do not mutate the arrays
  births=births2[p];
  # now that all the faces are sorted, we can go through the list
  # of those and since duplicated faces are adjacent, we can quickly compile the list of things to get rif of
  indx=Int[]
# now we loop all faces
i=2; last_i=1; last_face=faces[1];
while i<=N
              # here we scan for repeated faces until they are not repeated
              inner_loop_worked=false
              while isequal(faces[i],last_face)
                    inner_loop_worked=true
                    if i==N
                       break
                    else
                         i+=1
                    end
              end
              # Now we have reached the first non-repeated face, or just reached i==N
              if inner_loop_worked  # this means that during the last inner loop we had to skip some repeated faces
                                    # in particular, i points now to the first non-repeated face

                    if i<N # this is when the i-th element is distinct from the previous one
                       births[last_i]=minimum(births[last_i:i-1])
                       append!(indx,last_i+1:i-1) # add the list of indices to exclude
                    else # i.e. we ended up at  i==N
                          deletion_range_end= isequal(faces[N],last_face) ? N : N-1
                          births[last_i]=minimum(births[last_i:deletion_range_end])
                          append!(indx,last_i+1:deletion_range_end)
                          break # finally, since we reached the end we may break the loop here
                    end   # if i<N
                last_i=i-1
                else
                  last_i=i# since we know the i-th face is not duplicated we can start new checks here
                  i+=1
            end  #   if inner_loop_worked
              last_face=faces[last_i]
end # of outer while
deleteat!(births,indx);
return indx, p, births  # return all the information
end



 """
     Skeleton_a(FS::FiltrationOfSimplicialComplexes,dim::Int)::FiltrationOfSimplicialComplexes
     This Funcion takes a filtration of simplicial complexes and produces a filtration of their skeletons
     Usage:
     S=Skeleton_a(FS,max__dim);
 """
 function Skeleton_a(FS::FiltrationOfSimplicialComplexes,max__dim::Int)::FiltrationOfSimplicialComplexes

 if max__dim<=0; error("The maximal mimension needs to be positive"); end;
 if max__dim>MaximalHomologicalDimension;
   #error("This function is currently not designed to handle skeletons in dimension that is higher than $MaximalHomologicalDimension");
 end
 faces=sort.(collect.(FS.faces))     # here convert faces into an array form. This is done for improving speed
 is_below=FS.dimensions.<max__dim
 Index_of_Below_max__dim=findall(is_below)   # this array contains indices of faces of dimention< max__dim
 Index_of_high_dim=findall(.!is_below)       # this contains indices of faces of dimention>= max__dim
 N_small_faces= length(Index_of_Below_max__dim)  # the number of faces of dimention< max__dim
 There_are_large_facets= !isempty(Index_of_high_dim)


if There_are_large_facets
  # Now we create subsets of length  max__dim+1 for each Maximal face
  lengths_of_combinations= [ binomial(length(faces[i]) ,max__dim+1)  for i in Index_of_high_dim];
  N_large_faces=sum(lengths_of_combinations)

  # Now we assemble all these subsets in one array
  large_facets=Array{Array{TheIntegerType,1}}(undef,N_large_faces);
  large_facets_births= Array{Int}(undef,N_large_faces);
  counter=1
  for i=1:length(Index_of_high_dim)
    fac=faces[Index_of_high_dim[i]]
    b=FS.birth[Index_of_high_dim[i]]
    for f in combinations(fac, max__dim+1)
        large_facets[counter]=f
        large_facets_births[counter]=b
        counter+=1
    end
  end

  ###########################################
  # now we delete duplicates, while preserving correct births
  indx, p, large_facets_births = RemoveDuplicates(large_facets,large_facets_births);
  large_facets=large_facets[p];
  deleteat!(large_facets,indx);
  ##########################################
else # i.e. there are no large facets of size max__dim in this skeleton
    large_facets=Array{TheIntegerType,1}[]
    large_facets_births =Int[]
end

# construct the array of low-dim faces here
small_facets= FS.faces[Index_of_Below_max__dim]
small_facets_births= FS.birth[Index_of_Below_max__dim]
# colect the pieces together
AllFaces=vcat(small_facets,CodeWord.(large_facets));
AllBirths= vcat(small_facets_births,large_facets_births);
p=sortperm(AllBirths);
# Finally, we return the resulting filtration of simplicial complexes.
# Note that here we use "unsafe" version of the constructor, that does not check sanity of this
# filtration, since our construction is designed to produce a sane result
return FiltrationOfSimplicialComplexes(AllFaces[p], AllBirths[p],FS.vertices)
end
####################












"""
       Usage:
       TheSkeleton, f_vector=Skeleton_and_fvector(FS,max__dim);
       This Funcion takes a filtration of simplicial complexes and produces a filtration of their skeletons, as well as
       the f_vecors for each dimension in 1:max__dim (it does not outpur f_0, because it is useless in most use cases)
"""
       function Skeleton_and_fvector(FS::FiltrationOfSimplicialComplexes,max__dim::Int)::Tuple{FiltrationOfSimplicialComplexes,Vector{Vector{Int}}}
      # first, we compute the skeleton
        TheSkeleton=Skeleton_a(FS,max__dim);
       # Next, we create the  f_vector
       unique_births=unique(TheSkeleton.birth);
       N_births=length(unique_births);
       f_vector=Vector{Vector{Int}}(undef,max__dim);
       # first we compute f_vector[end]
       f_vector[max__dim]=zeros(Int,N_births);
       inx= findall(TheSkeleton.dimensions.==max__dim);
       N_facets=length(inx);
       f_vector[end]=[sum(TheSkeleton.birth[ inx].<=fill(unique_births[i], N_facets)) for i=1:N_births];
       # now we loop over the lower dimensions
       cofaces=sort.(collect.(TheSkeleton.faces[inx]));
       coface_births=TheSkeleton.birth[inx];
       for d=max__dim-1:-1:1
           ind_d_faces= findall(TheSkeleton.dimensions.==d);
           N_d_faces=length(ind_d_faces)
           # Now we append to d_dim_faces all the d-dim faces stemming from the cofaces
           N=length(ind_d_faces)+(d+2)*length(cofaces)
           faces=Array{Array{TheIntegerType,1}}(undef,N);
           births= Array{Int}(undef,N);
           faces[1:N_d_faces]=sort.(collect.(TheSkeleton.faces[ind_d_faces]));
           births[1:N_d_faces]=TheSkeleton.birth[ind_d_faces];
           counter=N_d_faces+1;
            # now we fill-out the rest of the arrays faces and births
            for i=1:length(cofaces)
              cf=cofaces[i]
              births[counter:counter+d+1].=coface_births[i]
              for j=1:d+2
                  c=copy(cf)
                  faces[counter]=deleteat!(c,j)
                  counter+=1
              end
            end
        # next we remove duplicates and adjust the births
        if !isempty(faces)
          indx,p,births = RemoveDuplicates(faces,births);
          faces=faces[p];
          deleteat!(faces,indx);
        end
        # fote that the "cleaned up" faces and births are ordered according to faces (and not births!)
        # next, we compute the d-th f-vector
        f_vector[d]=[sum(births.<=fill(unique_births[i], length(births))) for i=1:N_births];
        # Finally, if the loop is going to continue, then we need to update cofaces, and coface_births
          if d>1
          cofaces=faces
          coface_births=births
          end
        end # d=max__dim-1:-1:1
       return TheSkeleton, f_vector
       end
       ####################
























function TestComplex(m::Int,n::Int,Threshhold=0.45)
A=rand(m,n);
K,r=Simplicial.DowkerComplex(A)
filt_index= findfirst(r.>Threshhold)-1
deleteat!(K.faces, filt_index:K.depth)
deleteat!(K.dimensions, filt_index:K.depth)
deleteat!(K.birth, filt_index:K.depth)
K.depth=maximum(K.birth);
return K
end

#######
