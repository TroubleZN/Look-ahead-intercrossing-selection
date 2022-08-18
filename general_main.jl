using Distributed
using StatsBase
using SharedArrays

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

using ProgressMeter
function general_main_s(elite,donors,BCs,eft,loc,RF,select_method;n_cross = 2,k=200,Goal = 0.98)
    L,_,N = size(donors)
    res = zeros(N)
    gen_T = 10
    
    RC_Res = ones(N,10)
    res = zeros(N,2)
    Pairs = zeros(Int,2,n_cross)
    progs = zeros(L,2,k*n_cross)
    @showprogress for i in 1:N
        # i = sample(1:N,1)[1]
        gen = 1
        REC_rate = 0.0

        G_all = elite
        G_all = cat(G_all,[elite[:,1] donors[:,1,i]],dims = 3)
        BC1 = BCs[string(i)]
        i_pos = check_positive(BC1,loc)
        G_new = BC1[:,:,i_pos]        
        G_all = cat(G_all,G_new,dims = 3)
        REC_rates = Recovery_rate(G_new,eft,loc)
        REC_rate = maximum(REC_rates)
        RC_Res[i,Int(gen)] = REC_rate

        loc_val = sum(G_new[loc,:,:],dims = (1,2))[:]
        G_perf = Recovery_rate(G_new[:,:,loc_val.==2*length(loc)],eft,loc)
        if length(G_perf)==0
            G_perf = 0
        else
            G_perf = maximum(G_perf)
        end
        # RC_Res[i,Int(gen)] = G_perf


        while gen<gen_T && G_perf< Goal #(gen<=gen_T && REC_rate < (L*2-length(loc))/(L*2)) && (G_perf<0.97)
            gen+=1
            Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
            for j in 1:size(Pairs,2)
                progs[:,:,(1:k).+(j-1)*k] = cross(G_all[:,:,Pairs[1,j]],G_all[:,:,Pairs[2,j]],RF,k)
            end
            i_pos = check_positive(progs,loc)
            # G_new = progs
            G_new = progs[:,:,i_pos]
            G_all = cat(G_all,G_new,dims = 3)

            # REC_rate = maximum(Recovery_rate(G_new,eft,loc))
            i_pos = check_positive(G_new,loc)
            REC_rates = Recovery_rate(G_new,eft,loc)
            REC_rates[setdiff(1:size(G_new,3),i_pos)] .= 0
            REC_rate = maximum(REC_rates)
            RC_Res[i,Int(gen)] = REC_rate
            # print(REC_rate,"\t")
            
            loc_val = sum(G_new[loc,:,:],dims = (1,2))[:]
            G_perf = Recovery_rate(G_new[:,:,loc_val.==2*length(loc)],eft,loc)
            if length(G_perf)==0
                G_perf = 0
            else
                G_perf = maximum(G_perf)
            end
            display(Pairs)
            print(REC_rate," ",G_perf," ",gen,"\n")
            # RC_Res[i,Int(gen)] = G_perf
        end
        res[i,:] = [gen,REC_rate]
    end
    return res,RC_Res
end


function general_main_p(elite,donors,BCs,eft,loc,RF,select_method;n_cross = 2,k=200,Goal = 0.98)
    L,_,N = size(donors)
    gen_T = 10
    
    RC_Res = SharedArray(ones(N,gen_T))
    res = SharedArray(zeros(N,2))

    loc_score = SharedArray(zeros(N,gen_T))
    num_pos = SharedArray(zeros(N,gen_T).+k*n_cross)

    @sync @showprogress @distributed for i in 1:N

        Pairs = zeros(Int,2,n_cross)
        progs = zeros(L,2,k*n_cross)

        gen = 1
        REC_rate = 0.0

        G_all = cat(elite,[elite[:,1] donors[:,1,i]],dims = 3)
        BC1 = BCs[string(i)]
        i_pos = check_positive(BC1,loc)
        G_new = BC1[:,:,i_pos]        
        G_all = cat(elite,G_new,dims = 3)
        REC_rates = Recovery_rate(G_new,eft,loc)
        REC_rate = maximum(REC_rates)
        # REC_rate = mean(REC_rates.>0.95) 
        # REC_rate = var(REC_rates) 
                   
        RC_Res[i,Int(gen)] = REC_rate
        loc_score[i,Int(gen)] = maximum(sum(G_new[loc,:,:],dims = (1,2)))
        num_pos[i,Int(gen)] = size(G_new,3)

        loc_val = sum(G_new[loc,:,:],dims = (1,2))[:]
        G_perf = Recovery_rate(G_new[:,:,loc_val.==2*length(loc)],eft,loc)
        if length(G_perf)==0
            G_perf = 0
        else
            G_perf = maximum(G_perf)
        end
        
        # RC_Res[i,Int(gen)] = G_perf

        while gen<gen_T && G_perf< Goal #(gen<=gen_T && REC_rate < (L*2-length(loc))/(L*2)) && (G_perf<0.97)
            gen+=1
            ## Pairs: (2*n_cross)

            i_positive = check_positive(G_new,loc)
            if length(i_positive)>=0
                RCs = Recovery_rate(G_new[:,:,i_positive],eft,loc)
                if sum(RCs.>=Goal)>=n_cross #min(Goal+0.01,1)
                    Pairs[1,:] = i_positive[partialsortperm(RCs,1:n_cross,rev = true)]
                    Pairs[2,:] = Pairs[1,:]
                    # Pairs .= i_positive[findmax(RCs)[2]]
                    Pairs .+= 1
                else
                    Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
                end
            else
                Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
            end

            # Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
            # Pairs = Int.(Pairs)
            for j in 1:size(Pairs,2)
                progs[:,:,(1:k).+(j-1)*k] = cross(G_all[:,:,Pairs[1,j]],G_all[:,:,Pairs[2,j]],RF,k)
            end
            i_pos = check_positive(progs,loc)
            G_new = progs[:,:,i_pos]
            G_all = cat(elite,G_new,dims = 3)

            # REC_rate = maximum(Recovery_rate(G_new,eft,loc))
            i_pos = check_positive(G_new,loc)
            REC_rates = Recovery_rate(G_new,eft,loc)
            REC_rates[setdiff(1:size(G_new,3),i_pos)] .= 0
            # REC_rate = mean(sort(REC_rates,rev=true)[1:10])
            REC_rate = maximum(REC_rates)
            # REC_rate = mean(REC_rates.>0.95) 
            # REC_rate = var(REC_rates) 
           
            RC_Res[i,gen] = REC_rate
            loc_score[i,Int(gen)] = maximum(sum(G_new[loc,:,:],dims = (1,2)))
            num_pos[i,Int(gen)] = size(G_new,3)

            loc_val = sum(G_new[loc,:,:],dims = (1,2))[:]
            G_perf = Recovery_rate(G_new[:,:,loc_val.==2*length(loc)],eft,loc)
            if length(G_perf)==0
                G_perf = 0
            else
                G_perf = maximum(G_perf)
            end
            # display(Pairs)
            # print(REC_rate,"\t")

            # RC_Res[i,Int(gen)] = G_perf
        end
        res[i,:] = [gen,REC_rate]
    end
    return res,RC_Res,loc_score,num_pos
end

function general_main_p_inner(elite,donors,BCs,eft,loc,RF,select_method;n_cross = 2,k=200,Goal = 0.98)
    L,_,N = size(donors)
    gen_T = 10
    
    RC_Res = SharedArray(ones(N,10))
    res = SharedArray(zeros(N,2))

    loc_score = SharedArray(zeros(N,gen_T))
    num_pos = SharedArray(zeros(N,gen_T).+k*n_cross)

    @sync @showprogress @distributed for i in 1:N

        Pairs = zeros(Int,2,n_cross)
        progs = zeros(L,2,k*n_cross)

        gen = 1
        REC_rate = 0.0

        G_all = cat(elite,[elite[:,1] donors[:,1,i]],dims = 3)
        BC1 = BCs[string(i)]
        BC1 = cross(elite,[elite[:,1] donors[:,1,i]],RF,k)

        # i_pos = check_positive(BC1,loc)
        G_new = BC1
        G_all = cat(elite,G_new,dims = 3)
        REC_rates = Recovery_rate(G_new,eft,loc)
        REC_rate = maximum(REC_rates)
        # REC_rate = mean(REC_rates.>0.95) 
        # REC_rate = var(REC_rates) 
        
        RC_Res[i,Int(gen)] = REC_rate
        loc_score[i,Int(gen)] = maximum(sum(G_new[loc,:,:],dims = (1,2)))
        num_pos[i,Int(gen)] = size(G_new,3)

        loc_val = sum(G_new[loc,:,:],dims = (1,2))[:]
        G_perf = Recovery_rate(G_new[:,:,loc_val.==2*length(loc)],eft,loc)
        if length(G_perf)==0
            G_perf = 0
        else
            G_perf = maximum(G_perf)
        end
        
        # RC_Res[i,Int(gen)] = G_perf

        while gen<gen_T && G_perf< Goal #(gen<=gen_T && REC_rate < (L*2-length(loc))/(L*2)) && (G_perf<0.97)
            gen+=1
            ## Pairs: (2*n_cross)

            i_positive = check_positive(G_new,loc)
            if length(i_positive)>0
                RCs = Recovery_rate(G_new[:,:,i_positive],eft,loc)
                if sum(RCs.>=Goal)>=n_cross #min(Goal+0.01,1)
                    Pairs[1,:] = i_positive[partialsortperm(RCs,1:n_cross,rev = true)]
                    Pairs[2,:] = Pairs[1,:]
                    # Pairs .= i_positive[findmax(RCs)[2]]
                    Pairs .+= 1
                else
                    Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
                end
            else
                Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
            end

            # Pairs .= Select_pair(G_all,G_new,eft,loc,RF,select_method,n_cross,k,gen,Goal = Goal)
            # Pairs = Int.(Pairs)
            for j in 1:size(Pairs,2)
                progs[:,:,(1:k).+(j-1)*k] = cross(G_all[:,:,Pairs[1,j]],G_all[:,:,Pairs[2,j]],RF,k)
            end
            i_pos = check_positive(progs,loc)
            G_new = progs
            G_all = cat(elite,G_new,dims = 3)

            # REC_rate = maximum(Recovery_rate(G_new,eft,loc))
            i_pos = check_positive(G_new,loc)
            REC_rates = Recovery_rate(G_new,eft,loc)
            REC_rates[setdiff(1:size(G_new,3),i_pos)] .= 0
            # REC_rate = mean(sort(REC_rates,rev=true)[1:10])
            REC_rate = maximum(REC_rates)
            # REC_rate = mean(REC_rates.>0.95) 
            # REC_rate = var(REC_rates) 
           
            RC_Res[i,gen] = REC_rate
            loc_score[i,Int(gen)] = maximum(sum(G_new[loc,:,:],dims = (1,2)))
            num_pos[i,Int(gen)] = size(G_new,3)

            loc_val = sum(G_new[loc,:,:],dims = (1,2))[:]
            G_perf = Recovery_rate(G_new[:,:,loc_val.==2*length(loc)],eft,loc)
            if length(G_perf)==0
                G_perf = 0
            else
                G_perf = maximum(G_perf)
            end
            # display(Pairs)
            # print(REC_rate,"\t")

            # RC_Res[i,Int(gen)] = G_perf
        end
        res[i,:] = [gen,REC_rate]
    end
    return res,RC_Res,loc_score,num_pos
end

function general_main(elite,donors,BCs,eft,loc,RF,select_method;n_cross = 2,k=200,Goal = 0.98,multi_core=false)
    if select_method in ["PCV","PCV_nway","select_sim","LAS"]
        return general_main_p_inner(elite,donors,BCs,eft,loc,RF,select_method;n_cross = n_cross,k=k,Goal = Goal)
    end
    if multi_core == false
        return general_main_s(elite,donors,BCs,eft,loc,RF,select_method;n_cross = n_cross,k=k,Goal = Goal)
    else
        return general_main_p(elite,donors,BCs,eft,loc,RF,select_method;n_cross = n_cross,k=k,Goal = Goal)
    end
end
