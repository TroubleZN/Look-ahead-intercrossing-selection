###
using VCFTools
using CSV
using DataFrames
using JLD

N = 384
N_QTL = 3

###############
## load data ##
###############
vcf_path = string("Data/New_09172021/introgression_",N,"/phased_SNPs_",N,".vcf")
eft_path = string("Data/New_09172021/introgression_",N,"/thresholded_volume_effects_",N,".csv")
pheno_path = string("Data/New_09172021/introgression_",N,"/phenotypes_",N,".csv")
map_path = string("Data/New_09172021/introgression_",N,"/genetic_map_",N,".csv")

vcf_path = string("Data/New_09232021/introgression_",N,"_enrich/phased_SNPs_",N,"_enriched.vcf")
eft_path = string("Data/New_09232021/introgression_",N,"_enrich/thresholded_volume_effects_",N,"_enriched.csv")
pheno_path = string("Data/New_09232021/introgression_",N,"_enrich/phenotypes_",N,"_enriched.csv")
map_path = string("Data/New_09232021/introgression_",N,"_enrich/genetic_map_",N,"_enriched.csv")

## Genoa
H = convert_ht(Float64, vcf_path)'
G = reshape(H,size(H,1),2,:)

## effectaa
effect = CSV.read(eft_path,DataFrame)
names(effect)
sum(effect.Major1)
eft = effect.Effect
Tracked = effect.Tracked

## pheno
pheno = CSV.read(pheno_path,DataFrame)
names(pheno)
EXPVP = pheno.EXPVP

## genetic map
map = CSV.read(map_path,DataFrame)
names(map)
RF = map.RecombRate[2:end]

(effect.SNP == map.SNP)

##########
## main ##
##########
# eft = eft[Tracked.==1]
# G = G[Tracked.==1,:,:]
G[eft.<0,:,:] = abs.(G[eft.<0,:,:].-1)
eft = abs.(eft)

elites=G[:,:,EXPVP.==1]
donors=G[:,:,EXPVP.==0]

loc = partialsortperm(eft,1:N_QTL,rev = true)

donors_sub=donors[loc,:,:]
elites_sub=elites[loc,:,:]

donors_N_QTL=donors[:,:,findall(sum(donors_sub,dims = (1,2))[:].==N_QTL*2)]
elites_0_QTL=elites[:,:,findall(sum(elites_sub,dims = (1,2))[:].==0)]
elites_N_QTL=elites[:,:,findall(sum(elites_sub,dims = (1,2))[:].==N_QTL*2)]

val,ind=findmax(sum(elites_0_QTL,dims = (1,2))[:])

elite=elites_0_QTL[:,:,ind]

zero_locs=findall(elite[:,1].==0)
deleteat!(zero_locs,zero_locs .∈ (loc,))

elite[zero_locs,:] = 1 .-elite[zero_locs,:]


RF = RF./2 ##################

donors = donors_N_QTL
effect = eft
rf = RF

# donors = zeros(size(donors)) ######################
# donors[loc,:,:] .= 1 ######################

L,_,N = size(donors)
BCs = Dict()
for i in 1:N
    BC1 = cross(elite,[elite[:,1] donors[:,1,i]],RF,200)
    BCs[string("BC_",i)] = BC1
end

Data_3QTL = Dict(
    "elite" => elite,
    "loc" => loc[:],
    "effect" => eft,
    "donors" => donors,
    "RF" => RF,
)

save("Data/Data_3QTL.jld",Data_3QTL)

using MAT
matwrite("matfile.mat", Dict(
    "elite" => elite,
    "loc" => loc,
    "effect" => eft,
    "donors" => donors,
    "RF" => RF
    ); compress = true)

matwrite("Data_3QTL.mat", Data_3QTL)

matwrite("BCs.mat", BCs)

## load data
# using MAT

# Data = matread("./Data/data_4QTL.mat")

# elite = Data["elite"]
# loc = Int.(Data["loc"])
# effect = Data["effect"]
# donors = Data["donors"]

# RF = matread("./Data/RF.mat")["RF"]

