include("Code/breed_fun.jl")
include("Code/BGS.jl")
include("Code/LAMC.jl")

using StatsBase

function BC_main(elite,donors,eft,loc,RF,select_method;n_pair = 2,k=200)
    L,_,N = size(donors)
    res = zeros(N)
    REC_rate = 0
    gen_T = 10
    
    prog = zeros(L)
    res = zeros(N,2)
    
    @showprogress for i in 1:N
        F1 = [elite[:,1] donors[:,1,i]]
        gen = 1.0
        REC_rate = 0
        while REC_rate < 0.99 && gen <= gen_T
            if gen == 1
                Rand_g = random_garmets(F1,RF,k)
            else
                Rand_g = random_garmets(prog[1],RF,k)
                for j in 2:n_pair
                    Rand_g = [Rand_g random_garmets(prog[j],RF,k)]
                end
            end
            Positive_g = Rand_g[:,findall(sum(Rand_g[loc[:],:],dims=1)[:].==length(loc))]
            g_selected = Select_g(Positive_g,elite,eft,loc,RF,method=select_method,n_pair=n_pair,T = gen_T-gen+1)
            prog = [[elite[:,1] g_selected[:,i]] for i in 1:size(g_selected,2)]

            REC_rate = maximum([Recovery_rate(prog[i],eft,loc) for i in 1:length(prog)])
            gen+=1.0
        end
        res[i,:] = [gen,REC_rate]
        # print(i,"\n")
    end
    return res
end

method = "Random"
method = "BGS"
method = "LAMC"
method = "PCV"

N_pairs = [2,4,6,8,10]
RES = zeros(90,5)

for i in eachindex(N_pairs)
    res = BC_main(elite,donors,effect,loc,RF*0.2,method;n_pair = N_pairs[i],k=400)
    RES[:,i] = res[:,1]
end

p = res_plot(RES,method)
