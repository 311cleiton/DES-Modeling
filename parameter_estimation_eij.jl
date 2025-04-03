print(@__FILE__,"\n")
@time begin
# ESTIMATION OF εᵢⱼ (and kᵢⱼ)
# https://github.com/ClapeyronThermo/Clapeyron.jl
# https://clapeyronthermo.github.io/Clapeyron.jl/dev/api/estimation/
# https://github.com/jmejia8/Metaheuristics.jl

# STEP 1 - PACKAGES
using Clapeyron
using Metaheuristics
using Statistics
using Printf
using Dates
print(now(),"\n")

# STEP 2 - GENERATE MODEL
species_all = [
"ChCl_EG_2"     # 1
"TPPPBr_Ph_6"   # 2
"TBPBr_Ph_4"    # 3
"TPPPBr_Ph_4"   # 4
"TBACl_CA_2"    # 5
"MTACl_CA_2"    # 6
"MTABr_CA_2"    # 7
"TOACl_CA_2"    # 8
"TOABr_CA_2"    # 9
"TOACl_CA_1.5"  # 10
] # CSV Clapeyron database
choose_species = 10
species2 = species_all[choose_species]
species1 = "CO2_unified" 
# model_eij = PCSAFT([species1,species2])
model_eij = SAFTVRMie([species1,species2])
print(model_eij,"\n")

# PARAMETERS TO BE FITTED
toestimate = [
    Dict(
        :param => :epsilon,
        :indices => (1,2), # (1,1)=species1(CO₂); (2,2)=species2(DES)
        :symmetric => true, # true → epsilon(1,2)=epsilon(2,1)
        # :recombine => true, # true → use combining rules
        :recombine => false,
        :lower => model_eij.params.epsilon[1,2]-10.,
        :upper => model_eij.params.epsilon[1,2]+10.,
        :guess => model_eij.params.epsilon[1,2]
    ),
];

# PROPERTIES TO FIT
function bubble_point(model_eij::EoSModel,T,x)
    bub = bubble_pressure(model_eij,T,[x,1-x])
    return bub[1]
end

# ESTIMATOR
    # EXPERIMENTAL DATA CSV
    file_path_VLE = "C:/Users/beral/My Drive/DTU/PROJECTS/BS5/CSV/kij/"*species2*".csv"
estimator,objective,initial,upper,lower = Estimation(model_eij,toestimate,[file_path_VLE]);
method = ECA(;options=Options(iterations=100));

local_params, local_model = optimize(objective, estimator, method);

eij_0 = model_eij.params.epsilon[1,2]
eij = local_params[1]
kij = 1 - eij / eij_0

# PRINT ALL
print("εᵢⱼ_0, εᵢⱼ_fit, kᵢⱼ:","\n")
@printf("%.4f",eij_0)
print("\n")
@printf("%.4f",eij)
print("\n")
@printf("%.4f",kij)
print("\n")

clipboard(eij)

print("Runtime: ")
end#@time
