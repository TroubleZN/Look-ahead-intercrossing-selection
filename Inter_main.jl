include("Code/breed_fun.jl")
include("Code/BGS.jl")
include("Code/LAMC.jl")
include("Code/waterpip.jl")
using StatsBase

function getPCV(elite,donors,eft,RF)
    G = cat(elite,donors,dims=(3))
    GEBVs = GEBV(G,eft)
    i_best = sortperm(GEBVs[:],rev = true)[1:min(length(GEBVs),10)]
    G_best = G[:,:,i_best]

    L,_,N = size(G)
    PCVs = zeros(N,N)
    M = getM(RF);
    temp = zeros(4)
    X = zeros(4,size(elite,1))
    W = zeros(size(X))
    @inbounds for i in 1:length(i_best)
        @inbounds for j in 1:N
            G1 = view(G_best,:,:,i)
            G2 = view(G,:,:,j)
            getXW!(X,W,G1,G2)
            PCVs[i_best[i],j],_ = getwater(M,X,W,temp)
        end
    end
    return G,PCVs
end

function check_OHV(G1,G2,loc)
    for i in loc[:]
        if sum(G1[i,:])+sum(G2[i,:]) == 0
            return false
        end
    end
    return true
end


function Inter_main(elite,donors,eft,loc,RF,select_method;n_pair = 2,k=200)
    L,_,N = size(donors)
    res = zeros(N)
    REC_rate = 0
    gen_T = 10
    
    prog = zeros(L)
    res = zeros(N,2)
    
    @showprogress for i in 1:N

        G = [elite[:,1] donors[:,1,i]]
        gen = 1.0
        REC_rate = 0
        prog = zeros(L,2,n_pair*k)

        while REC_rate < 0.99 && gen <= gen_T
            G_all,PCVs = getPCV(elite,G,eft,RF)
            if gen == 1.0
                _,pair = findmax(PCVs)
                PCVs[pair] = 0.0
                PCVs[pair[2],pair[1]] = 0.0
                prog[:,:,1:k] .= cross(G_all[:,:,pair[1]], G_all[:,:,pair[2]], RF, k)
                G = cat(G,prog,dims = 3)
            else
                for j in 1:n_pair
                    pair = [1,1]
                    while !check_OHV(G_all[:,:,pair[1]],G_all[:,:,pair[2]],loc)
                        _,pair = findmax(PCVs)
                        PCVs[pair] = 0
                        PCVs[pair[2],pair[1]] = 0.0
                    end
                    prog[:,:,(1:k).+k*(j-1)] .= cross(G_all[:,:,pair[1]], G_all[:,:,pair[2]], RF, k)
                    G = cat(G,prog,dims = 3)
                    # print(pair[1],"\t",pair[2],"\n")
                end
            end
            REC_rate = maximum([Recovery_rate(prog[:,:,i],eft,loc) for i in 1:size(prog,3)])
            gen+=1.0
            # print(gen,"\t",REC_rate,"\n")
        end
        # i_positive = findall(sum(prog[loc[:],:,:],dims = (1,2))[:].==length(loc)*2)
        # prog_pos = prog[:,:,i_positive]
        # REC_rate = maximum([Recovery_rate(prog_pos[:,:,i],eft,loc) for i in 1:size(prog_pos,3)])

        res[i,:] = [gen,REC_rate]
    end
    return res
end



method = "Inner"

N_pairs = [2,4,6,8,10]
RES = zeros(90,5)

for i in eachindex(N_pairs)
    res = Inter_main(elite,donors,effect,loc,RF,method;n_pair = N_pairs[i],k=200)
    RES[:,i] = res[:,1]
end

res_plot(RES,method)
