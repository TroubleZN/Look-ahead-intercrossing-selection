using Distributed

if nworkers()==1
    addprocs(12)
end

# include("data_pre.jl")

@everywhere begin
    cd("Code_fixed_generation")
    include("main_N_generation.jl")
    cd("../")

    # using MAT
    # Data = matread("Data\\data_3QTL.mat")
    # elite = Data["elite"]
    # loc = Int.(Data["loc"][:])
    # eft = Data["effect"][:]
    # donors = Data["donors"]
    # RF = matread("Data\\RF.mat")["RF"]
    Data = load("Data/Data_3QTL.jld")
    elite = Data["elite"]
    loc = Data["loc"]
    eft = Data["effect"][:]
    donors = Data["donors"]
    RF = Data["RF"]
    RF[RF.==0] .= 0.5
    BCs = Data["BCs"]

    # Methods = ["BGS_BC","PCV_BC","BGS","PCV","PCV_sim","PCV_sim2"]
    Methods = ["BGS_BC","PCV_BC","BGS","PCV","PCV_nway","PCV_sim","PCV_sim2","select_sim2"]
    Methods = ["BGS_BC","PCV_BC","PCV","PCV_nway"]
    # Methods = ["BGS_BC","PCV_BC","BGS","PCV","PCV_nway","select_sim2","select_sim3","select_sim"]
    # Methods = ["BGS_BC","select_sim3","PCV_nway"]
    Methods = ["BGS_BC","PCV_BC","PCV","PCV2","PCV_sim","PCV_sim2","PCV_nway","PCV_nway_imperf","LAMC"]
    Methods = ["BGS_BC","PCV_BC","LAMC"]
    Methods = ["BGS-IC","PCV_IC"] 
    
    Goal = 0.96
    k = 200
end

plot_dir = string("./server2/res_N_generation_",Goal,"_",k,"/")
if !isdir(plot_dir)
    mkdir(plot_dir)
end
using JLD

for method in Methods
    print("\n",method,"\n")
    N_pairs = 2:6
    for i in eachindex(N_pairs)
        print("N_cross = ",N_pairs[i],"\n")
        Num_Goal,RC_max,RC_mean,GT,GT_1 = n_generation_main_p(elite, donors, BCs, eft, loc, RF, method;gen_T = 5,n_cross = N_pairs[i],k=k,Goal = Goal)
        save(string(plot_dir,method,"_",N_pairs[i],"_Num_Goal.jld"),"Num_Goal",Array(Num_Goal))
        save(string(plot_dir,method,"_",N_pairs[i],"_RC_max.jld"),"RC_max",Array(RC_max))
        save(string(plot_dir,method,"_",N_pairs[i],"_RC_mean.jld"),"RC_mean",Array(RC_mean))
        save(string(plot_dir,method,"_",N_pairs[i],"_GT.jld"),"GT",Array(GT))
        save(string(plot_dir,method,"_",N_pairs[i],"_GT_1.jld"),"GT_1",Array(GT_1))
        display(mean(RC_max,dims = 1))
        display(mean(Num_Goal,dims = 1))
    end
end


# M1 = "PCV_nway"
# M2 = "select_sim3"
# print("\n",M1,"\n")
# res1 = load(string(plot_dir,M1,"_RC.jld"))["RC"]
# display(res1)
# print("\n",M2,"\n")
# res2 = load(string(plot_dir,M2,"_RC.jld"))["RC"]
# display(res2)
# res2-res1

