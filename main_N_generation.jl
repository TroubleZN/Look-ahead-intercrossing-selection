using Distributed
using StatsBase
using SharedArrays
using JLD

include("breed_fun.jl")
include("LAMC.jl")
include("waterpip.jl")
include("res_plot.jl")
include("select_method.jl")
include("select_sim.jl")
include("PCV_2way.jl")
include("BGS.jl")
include("LAS.jl")

function check_positive(G,loc)
    i_pos = findall(sum(maximum(G[loc[:],:,:],dims = 2),dims = 1)[:].==length(loc))
    return i_pos
end

function check_perf(G,loc)
    i_pos = findall(sum(G[loc[:],:,:],dims = (1,2))[:].==2*length(loc))
    return i_pos
end

using ProgressMeter
function n_generation_main_p(elite,donors,BCs,eft,loc,RF,select_method;gen_T = 5,n_cross = 2,k=200,Goal = 0.95)
    
    L,_,N = size(donors)
    
    Num_Goal = SharedArray(zeros(N,gen_T))
    RC_max = SharedArray(zeros(N,gen_T))
    RC_mean = SharedArray(zeros(N,gen_T))
    GT_1 = SharedArray(zeros(L,2,k*n_cross,N))
    GT = SharedArray(zeros(L,2,k*n_cross,N))

    @sync @showprogress @distributed for i in 1:N

        Pairs = zeros(Int,2,n_cross)
        progs = zeros(L,2,k*n_cross)
        G_all = zeros(L,2,3+k*(gen_T-2))

        G_all[:,:,1] = elite
        G_all[:,:,2] = donors[:,:,i]
        G_all[:,:,3] = [elite[:,1] donors[:,1,i]]
        RC_max[i,1] = Recovery_rate(G_all[:,:,3],eft,loc)
        Num_Goal[i,1] = 0
        RC_mean[i,1] = RC_max[i,1]

        G_new = cross(elite,[elite[:,1] donors[:,1,i]],RF,k)  
        for gen in 2:gen_T-2
            # if gen == gen_T-2
            #     GT_1[:,:,:,i] = G_new
            # end
            i_pos = check_positive(G_new,loc)
            G_pos = G_new[:,:,i_pos]
            if select_method in ["BGS_BC","PCV_BC","LAMC","LAS"]
                G_new = G_pos
            end
            G_Pool = cat(elite,G_new,dims = 3)

            REC_rates = Recovery_rate(G_pos,eft,loc)
            RC_max[i,gen] = maximum(REC_rates)
            Num_Goal[i,gen] = sum(REC_rates.>Goal)/length(REC_rates)
            RC_mean[i,gen] = sum(REC_rates)/length(REC_rates)

            Pairs .= Select_pair(G_Pool,G_new,eft,loc,RF,select_method,n_cross,k,gen,gen_T,Goal = Goal)

            # Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
            # Pairs = Int.(Pairs)
            for j in 1:size(Pairs,2)
                progs[:,:,(1:k).+(j-1)*k] = cross(G_Pool[:,:,Pairs[1,j]],G_Pool[:,:,Pairs[2,j]],RF,k)
            end
            G_new = progs
        end
        GT[:,:,:,i] = G_new
        i_pos = check_positive(G_new,loc)
        G_pos = G_new[:,:,i_pos]
        REC_rates = Recovery_rate(G_pos,eft,loc)
        RC_max[i,gen_T-1] = maximum(REC_rates)
        Num_Goal[i,gen_T-1] = sum(REC_rates.>Goal)/length(REC_rates)
        RC_mean[i,gen_T-1] = sum(REC_rates)/length(REC_rates)

        ## selfing
        res_selfing = zeros(length(i_pos))
        for j in 1:length(i_pos)
            progs = cross(G_pos[:,:,j],G_pos[:,:,j],RF,5000)
            i_perf = check_perf(progs,loc)
            REC_rates = Recovery_rate(progs[:,:,i_perf],eft,loc)
            if length(i_perf) != 0
                res_selfing[j] = maximum(REC_rates)*sum(REC_rates.>Goal)
            end
        end
        i_select = partialsortperm(res_selfing,1:n_cross,rev = true)
        Pairs[1,:] =  i_select
        Pairs[2,:] =  i_select

        progs = zeros(L,2,k*n_cross)
        for j in 1:size(Pairs,2)
            progs[:,:,(1:k).+(j-1)*k] = cross(G_pos[:,:,Pairs[1,j]],G_pos[:,:,Pairs[2,j]],RF,k)
        end
        i_perf = check_perf(progs,loc)
        REC_rates = Recovery_rate(progs[:,:,i_perf],eft,loc)
        RC_max[i,gen_T] = maximum(REC_rates)
        Num_Goal[i,gen_T] = sum(REC_rates.>Goal)/length(REC_rates)
        RC_mean[i,gen_T] = sum(REC_rates)/length(REC_rates)
    end

    return Num_Goal,RC_max,RC_mean,GT,GT_1
end


## New pipeline with 3 cross in the first generation
function n_generation_main_p2(elite,donors,BCs,eft,loc,RF,select_method;gen_T = 5,n_cross = 2,k=200,Goal = 0.95)
    
    L,_,N = size(donors)
    
    Num_Goal = SharedArray(zeros(N,gen_T))
    RC_max = SharedArray(zeros(N,gen_T))
    RC_mean = SharedArray(zeros(N,gen_T))
    # G3 = SharedArray(zeros(L,2,k*n_cross,N))
    GT_1 = SharedArray(zeros(L,2,k*n_cross,N))
    GT = SharedArray(zeros(L,2,k*n_cross,N))

    @sync @showprogress @distributed for i in 1:N

        Pairs = zeros(Int,2,n_cross)
        progs = zeros(L,2,k*n_cross)
        G_all = zeros(L,2,3+k*(gen_T-2))

        G_all[:,:,1] = elite
        G_all[:,:,2] = donors[:,:,i]
        G_all[:,:,3] = [elite[:,1] donors[:,1,i]]
        RC_max[i,1] = Recovery_rate(G_all[:,:,3],eft,loc)
        Num_Goal[i,1] = 0
        RC_mean[i,1] = RC_max[i,1]

        G_new = cross(elite,[elite[:,1] donors[:,1,i]],RF,3k)
        G_Pool = cat(elite,G_new,dims = 3)

        for gen in 2:gen_T-2
            if gen == gen_T-2
                GT_1[:,:,:,i] = G_new
            end
            i_pos = check_positive(G_new,loc)
            G_pos = G_new[:,:,i_pos]
            if select_method in ["BGS_BC","PCV_BC","LAMC","LAS"]
                G_new = G_pos
            end
            G_Pool = cat(elite,G_new,dims = 3)

            REC_rates = Recovery_rate(G_pos,eft,loc)
            RC_max[i,gen] = maximum(REC_rates)
            Num_Goal[i,gen] = sum(REC_rates.>Goal)/length(REC_rates)
            RC_mean[i,gen] = sum(REC_rates)/length(REC_rates)

            Pairs .= Select_pair(G_Pool,G_new,eft,loc,RF,select_method,n_cross,k,gen,gen_T,Goal = Goal)

            # Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
            # Pairs = Int.(Pairs)
            for j in 1:size(Pairs,2)
                progs[:,:,(1:k).+(j-1)*k] = cross(G_Pool[:,:,Pairs[1,j]],G_Pool[:,:,Pairs[2,j]],RF,k)
            end
            G_new = progs
        end

        GT[:,:,:,i] = G_new
        i_pos = check_positive(G_new,loc)
        G_pos = G_new[:,:,i_pos]
        REC_rates = Recovery_rate(G_pos,eft,loc)
        RC_max[i,gen_T-1] = maximum(REC_rates)
        Num_Goal[i,gen_T-1] = sum(REC_rates.>Goal)/length(REC_rates)        
        RC_mean[i,gen_T-1] = sum(REC_rates)/length(REC_rates)

        ## selfing
        res_selfing = zeros(length(i_pos))
        for j in 1:length(i_pos)
            progs = cross(G_pos[:,:,j],G_pos[:,:,j],RF,5000)
            i_perf = check_perf(progs,loc)
            REC_rates = Recovery_rate(progs[:,:,i_perf],eft,loc)
            if length(i_perf) != 0
                res_selfing[j] = maximum(REC_rates)*sum(REC_rates.>Goal)
            end
        end
        i_select = partialsortperm(res_selfing,1:n_cross,rev = true)
        Pairs[1,:] =  i_select
        Pairs[2,:] =  i_select

        progs = zeros(L,2,k*n_cross)
        for j in 1:size(Pairs,2)
            progs[:,:,(1:k).+(j-1)*k] = cross(G_pos[:,:,Pairs[1,j]],G_pos[:,:,Pairs[2,j]],RF,k)
        end
        i_perf = check_perf(progs,loc)
        REC_rates = Recovery_rate(progs[:,:,i_perf],eft,loc)
        RC_max[i,gen_T] = maximum(REC_rates)
        Num_Goal[i,gen_T] = sum(REC_rates.>Goal)/length(REC_rates)
        RC_mean[i,gen_T] = sum(REC_rates)/length(REC_rates)
    end

    return Num_Goal,RC_max,RC_mean,GT,GT_1
end
