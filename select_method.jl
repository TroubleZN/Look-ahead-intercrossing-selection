
function Select_pair(G_all,G_new,eft,loc,RF,method,n_cross,k,gen,gen_T;Goal = 0.95)
    L,_,N = size(G_all)
    N_new = size(G_new,3)
    Pairs = zeros(Int,2,n_cross)
    elite = G_all[:,:,1]

    n_simu = 50

    eft_temp = copy(eft)
    eft_temp = ones(L)
    eft_temp[loc] .= (1-Goal)*2L/length(loc)

    if method == "Random"
        Pairs .= reshape(sample(1:N,2*n_cross),2,:)
    elseif method == "BGS-IC"
        Pairs .= BGS(cat(elite,G_new,dims = 3),eft,loc,n_cross)
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "PCV-IC"
        PCVs = getPCV(cat(elite,G_new,dims = 3),eft,RF)
        for i in 1:n_cross
            _,pair = findmax(PCVs)
            Pairs[1,i] = pair[1]
            Pairs[2,i] = pair[2]

            i_sel = [pair[1],pair[2]]
            i_sel = setdiff(i_sel,[1])
            PCVs[i_sel,:] .= 0
            PCVs[:,i_sel] .= 0
        end
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "PCV2"
        PCVs = getPCV2(cat(elite,G_new,dims = 3),RF,Goal,loc)
        for i in 1:n_cross
            _,pair = findmax(PCVs)
            Pairs[1,i] = pair[1]
            Pairs[2,i] = pair[2]

            i_sel = [pair[1],pair[2]]
            i_sel = setdiff(i_sel,[1])
            PCVs[i_sel,:] .= 0
            PCVs[:,i_sel] .= 0
        end
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "PCV_nway"
        if maximum(Recovery_rate(G_new,eft_temp,loc))<=Goal-0.05
            Pairs = PCV_nway(cat(elite,G_new,dims = 3),n_cross,RF,loc,Goal,gen_T-gen-1;n_elite_cross=15,eft = eft_temp)
            Pairs[Pairs.!=1] .+= (N-N_new-1)
        else
            PCVs = getPCV2(cat(elite,G_new,dims = 3),eft_temp,RF,Goal,loc)
            for i in 1:n_cross
                _,pair = findmax(PCVs)
                Pairs[1,i] = pair[1]
                Pairs[2,i] = pair[2]

                i_sel = [pair[1],pair[2]]
                i_sel = setdiff(i_sel,[1])
                PCVs[i_sel,:] .= 0
                PCVs[:,i_sel] .= 0
            end
            Pairs[Pairs.!=1] .+= (N-N_new-1)
        end
    elseif method == "PCV_nway_imperf"
        Pairs = PCV_nway_imperf(cat(elite,G_new,dims = 3),n_cross,RF,loc,Goal,gen_T-gen-1;n_elite_cross=15,eft = eft_temp)
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "select_sim"
        Pairs = select_sim(cat(elite,G_new,dims = 3),eft_temp,RF,loc,Goal,k=k,n_cross=n_cross,n_simu = n_simu,gen = gen_T-gen-1)
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "LAIC"
        Pairs = select_sim2(cat(elite,G_new,dims = 3),eft,RF,loc,Goal,k=k,n_cross=n_cross,n_simu = n_simu,gen = gen,gen_T = gen_T)
        # display(Pairs)
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "PCV_sim"
        Pairs = PCV_sim(cat(elite,G_new,dims = 3),eft_temp,RF,loc,Goal,gen_T-gen-1;k=k,n_cross=n_cross,n_simu = n_simu)
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "PCV_sim2"
        Pairs = PCV_sim2(cat(elite,G_new,dims = 3),eft_temp,RF,Goal,k=k,n_cross=n_cross,n_simu = 50)
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "Random_BC"
        Pairs[1,:] = sample((1:N_new).+(N-N_new),n_cross)
        Pairs[2,:] .= 1
    elseif method == "BGS_BC"
        Pairs[1,:] = BGS(G_new,eft,loc)[1:n_cross]
        Pairs .+= (N-N_new)
        Pairs[2,:] .= 1
    elseif method == "PCV_BC"
        elite = G_all[:,:,1]
        PCVs = zeros(N_new)
        for i in 1:N_new
            PCVs[i],_ = PCV(G_new[:,:,i],elite,RF)
        end
        Pairs[1,:] = sortperm(PCVs,rev=true)[1:n_cross]
        Pairs[1,:] .+= (N-N_new)
        Pairs[2,:] .= 1
    elseif method == "PCV2_BC"
        elite = G_all[:,:,1]
        PCVs = zeros(N_new)
        for i in 1:N_new
            PCVs[i],_ = PCV2(G_new[:,:,i],elite,RF,Goal)
        end
        Pairs[1,:] = sortperm(PCVs,rev=true)[1:n_cross]
        Pairs[1,:] .+= (N-N_new)
        Pairs[2,:] .= 1
    elseif method == "LABC"
        Pairs .= LMC(G_new,elite,RF,loc;k=k,n_cross=n_cross,gen = gen_T-gen-1,n_simu = n_simu,Goal = Goal)
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    elseif method == "LAS"
        selection= LAS_Complete(cat(elite,G_new,dims = 3), eft_temp, RF, n_cross, n_cross*k, gen_T-gen-1,Goal)
        Pairs = selection
        Pairs[Pairs.!=1] .+= (N-N_new-1)
    end
    return Pairs
end

# function OHV(G::Matrix,eft,RF;B=10)
#     G = G.*eft
#     i_pick = sort(sortperm(RF,rev = true)[1:(B-1)])
#     i_pick = [i_pick; size(G,1)]
#     ohv = maximum(sum(G[1:i_pick[1],:],dims = 1))
#     for i in 1:length(i_pick)-1
#         ohv+=maximum(sum(G[i_pick[i]+1:i_pick[i+1],:],dims = 1))
#     end
#     return ohv*2
# end

# function OHV(G::Array{Float64, 3},eft,RF;B=20)
#     L,_,N= size(G)
#     G = G.*eft
#     i_pick = sort(sortperm(RF,rev = true)[1:(B-1)])
#     i_pick = [i_pick; L]
#     OHVs = zeros(N)
#     for j in 1:N
#         OHVs[j] += maximum(sum(G[1:i_pick[1],:,j],dims = 1))
#         for i in 1:length(i_pick)-1
#             OHVs[j] += maximum(sum(G[i_pick[i]+1:i_pick[i+1],:,j],dims = 1))
#         end
#     end
#     return OHVs*2
# end

# function select_OHV(G1,G2,gen_T,eft,loc,RF;Rec_r = 1,n_cross = 2,k=400,n_simu=5)
#     L,_ = size(G1)    
#     res = zeros(n_simu)
#     Pairs = zeros(Int,2,n_cross)
#     for i in 1:n_simu
#         if gen_T==1
#             G_new = cross(G1,G2,RF,k)
#             res[i] = maximum(Recovery_rate(G_new,eft,loc))
#         else
#             prog = cross(G1,G2,RF,k)
#             for gen in 2:gen_T
#                 G_new = prog[:,:,check_positive(prog,loc)]
#                 # OHVs = OHV(G_new,eft,RF)
#                 OHVs = GEBV(G_new,eft)
#                 Pairs[:].=sortperm(OHVs,rev = true)[1:2*n_cross]
#                 prog = cross(G_new[:,:,Pairs[1,1]],G_new[:,:,Pairs[2,1]],RF,k)
#                 for j in 2:n_cross
#                     prog = cat(prog,cross(G_new[:,:,Pairs[1,j]],G_new[:,:,Pairs[2,j]],RF,k),dims = 3)
#                 end
#             end
#             res[i] = maximum(Recovery_rate(prog,eft,loc))
#         end
#     end
#     return mean(res)
# end

# function select_OHV(G1,G2,elite,gen_T,eft,loc,RF;Rec_r = 1,n_cross = 2,k=400,n_simu=5)
#     L,_ = size(G1)    
#     res = zeros(n_simu)
#     Pairs = zeros(Int,2,n_cross)
#     for i in 1:n_simu
#         if gen_T==1
#             G_new = cross(G1,G2,RF,k)
#             res[i] = maximum(Recovery_rate(G_new,eft,loc))
#         else
#             prog = cross(G1,G2,RF,k)
#             for gen in 2:gen_T
#                 G_new = prog[:,:,check_positive(prog,loc)]
#                 G_all = cat(elite,G_new,dims = 3)
#                 _,PCVs = getPCV(G_all,eft,RF)
#                 for i in 1:n_cross
#                     _,pair = findmax(PCVs)
#                     Pairs[1,i] = pair[1]
#                     Pairs[2,i] = pair[2]
#                     PCVs[pair] = 0
#                     PCVs[pair[2],pair[1]] = 0
#                 end
#                 prog = cross(G_all[:,:,Pairs[1,1]],G_all[:,:,Pairs[2,1]],RF,k)
#                 for j in 2:n_cross
#                     prog = cat(prog,cross(G_all[:,:,Pairs[1,j]],G_all[:,:,Pairs[2,j]],RF,k),dims = 3)
#                 end
#             end
#             res[i] = maximum(Recovery_rate(prog,eft,loc))
#         end
#     end
#     return mean(res)
# end
