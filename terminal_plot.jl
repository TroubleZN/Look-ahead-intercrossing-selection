
    Goal = 0.96
    N_pairs = 2:6
    Methods = ["BGS-BC","PCV-BC","LAS-BC","BGS-IC","PCV-IC","LAS-IC"]
    res_1k = string("G://My Drive//Lab//Trait introgression//server2//res_N_generation_",Goal,"_200//")
    # res_1k = "G://My Drive//Lab//Trait introgression//server2//res_stable//"
    plot_dir = string("G://My Drive//Lab//Trait introgression//server2//")

    i = 1
    res_gen = Array(load(string(res_1k,"terminal_gen_",Methods[i],".jld"),"res_gen"))
    p1 = res_plot(res_gen,Methods[i],N_pairs)
    savefig(p1,string(plot_dir,"terminal_gen_",Goal,"_",Methods[i],".png"))

    i = 2
    res_gen = load(string(res_1k,"terminal_gen_",Methods[i],".jld"),"res_gen")
    p2 = res_plot(res_gen,Methods[i],N_pairs,ylab = nothing,leg = nothing,fsize = (500,500))
    savefig(p2,string(plot_dir,"terminal_gen_",Goal,"_",Methods[i],".png"))

    i = 3
    res_gen = load(string(res_1k,"terminal_gen_",Methods[i],".jld"),"res_gen")
    p3 = res_plot(res_gen,Methods[i],N_pairs,ylab = nothing,leg = nothing,fsize = (500,500))
    savefig(p3,string(plot_dir,"terminal_gen_",Goal,"_",Methods[i],".png"))

    i = 4
    res_gen = load(string(res_1k,"terminal_gen_",Methods[i],".jld"),"res_gen")
    p4 = res_plot(res_gen,Methods[i],N_pairs,ylab = nothing,leg = nothing,fsize = (500,500))
    savefig(p4,string(plot_dir,"terminal_gen_",Goal,"_",Methods[i],".png"))

    i = 5
    res_gen = load(string(res_1k,"terminal_gen_",Methods[i],".jld"),"res_gen")
    p5 = res_plot(res_gen,Methods[i],N_pairs,ylab = nothing,leg = nothing,fsize = (500,500))
    savefig(p4,string(plot_dir,"terminal_gen_",Goal,"_",Methods[i],".png"))

    i = 6
    res_gen = load(string(res_1k,"terminal_gen_",Methods[i],".jld"),"res_gen")
    p6 = res_plot(res_gen,Methods[i],N_pairs,ylab = nothing,fsize = (670,500))
    p6 = res_plot(res_gen,Methods[i],N_pairs,ylab = nothing,leg = nothing,fsize = (500,500))
    savefig(p4,string(plot_dir,"terminal_gen_",Goal,"_",Methods[i],".png"))

    p = plot(p1, p4, p2, p5, p3, p6,size=(850,1000),layout = grid(3, 2), dpi=600)
    savefig(p,string(plot_dir,"terminal_gen_",Goal,".png"))


labels = string.(collect((1:5).+4)')
labels[end] = "failed"
reverse!(labels)

p_legend = groupedbar(
    [NaN NaN NaN NaN NaN],
    label=labels,
    legend=:outerright,
    color = reverse(collect(cgrad(:roma, 6, categorical = true))'),
    framestyle = :none
)

plot(p1,p2,p3,p4,p_legend,layout = (1,5),size = (1200,400))


method = Methods[i]
RES = res_gen
N_simu,_ = size(RES)
# N_pairs = [2,4,6,8,10]
v_min = minimum(RES)
v_min = 5
v_max = maximum(RES)

RES[RES.>v_min+5-1] .= v_min+5-1
RES = Int.(RES)
RES_trans = zeros(length(N_pairs),5)

for k in 1:size(RES_trans,1)
    for j in 1:size(RES_trans,2)
        RES_trans[k,j] = count(i->i==j+v_min-1,RES[:,k])/size(RES,1)
    end
end

reverse!(RES_trans,dims = 2)

labels = string.(collect((1:size(RES_trans,2)).+(v_min-1))')
labels[end] = "failed"
reverse!(labels)
p = groupedbar(
    RES_trans,
    bar_position = :stack,
    bar_width=0.7,
    xticks=(1:length(N_pairs), string.("(",1:length(N_pairs),")")),
    yticks=(0:0.1:1, string.(0:10:100,"%")),
    
    # xlabel = "Resource allocation plan",
    # ylabel = "Probability of success",
    label=labels,
    legend=:outerright,
    # legend=nothing,

    color = reverse(collect(cgrad(:roma, 6, categorical = true))'),
    title = string(method)
)

for i = 1:5
    if RES_trans[i,5] >0.1
        annotate!(i,RES_trans[i,5]/2,"5")
    end
    if RES_trans[i,4] >0.1
        annotate!(i,RES_trans[i,5]+RES_trans[i,4]/2,"6")
    end
    if RES_trans[i,4] >0.1
        annotate!(i,RES_trans[i,5]+RES_trans[i,4]+RES_trans[i,3]/2,"7")
    end
end
