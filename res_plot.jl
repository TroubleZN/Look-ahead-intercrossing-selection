# RES is a matrix N_simu*length(N_pairs), each cell of which is imputed by T_gen
# N_simu is the number of total simulations
# N_pairs is the number of crosses in each generation
# T_gen is number of generations needed to reach the 95% recovery rate for each simulation

using Plots
using StatsPlots

function res_plot(RES,method,N_pairs;xlab = "Resource allocation plan",ylab = "Probability of success",tit = method, leg = :outerright,fsize = (700,500))
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

    if method == "LABC"
        method = "LAS-BC"
    elseif method == "LAIC"
        method = "LAS-IC"
    elseif method == "BGS_BC"
        method = "BGS-BC"
    elseif method == "PCV_BC"
        method = "PCV-BC"
    end

    labels = string.(collect((1:size(RES_trans,2)).+(v_min-1))')
    labels[end] = "failed"
    reverse!(labels)
    p = groupedbar(
        RES_trans,
        bar_position = :stack,
        bar_width=0.7,
        tickfontsize = 12,
        labelfontsize = 12,
        legendfontsize = 12,
        annotatefontsize = 12,
        titlefontsize = 15,
        xticks=(1:length(N_pairs), string.("(",1:length(N_pairs),")")),
        yticks=(0:0.1:1, string.(0:10:100,"%")),
        label=labels,
        legend=leg,
        size = fsize,
        color = reverse(collect(cgrad(:roma, 6, categorical = true))'),
    )
    if !isnothing(tit)
        title!(tit)
    end
    if !isnothing(xlab)
        xlabel!(xlab)
    end
    if !isnothing(ylab)
        ylabel!(ylab)
    end

    frontsize_anno = 12

    for i = 1:5
        if RES_trans[i,5] >0.08
            annotate!(i,RES_trans[i,5]/2,text("5",frontsize_anno))
        end
        if RES_trans[i,4] >0.08
            annotate!(i,RES_trans[i,5]+RES_trans[i,4]/2,text("6",frontsize_anno))
        end
        if RES_trans[i,3] >0.08
            annotate!(i,RES_trans[i,5]+RES_trans[i,4]+RES_trans[i,3]/2,text("7",frontsize_anno))
        end
    end

    display(p)
    return p
end

