# ESTIMATION OF εᵢⱼ (and kᵢⱼ)
# https://github.com/ClapeyronThermo/Clapeyron.jl
# https://clapeyronthermo.github.io/Clapeyron.jl/dev/api/estimation/
# https://github.com/jmejia8/Metaheuristics.jl

# STEP 1 - PACKAGES
using Clapeyron
using Metaheuristics
using Printf
using Dates
print(now())
print("\n")

# STEP 2 - GENERATE MODEL
species_1 = "carbon dioxide" # PCSAFT
species_2 = "TOACl_CA_1.5"
model_eij = PCSAFT([species_1,species_2])
print(model_eij)
print("\n")

# STEP 3 - PARAMETERS TO FIT
toestimate = [
    Dict(
        :param => :epsilon,
        :indices => (1,2), # (1,1)=species_1; (2,2)=species_2
        :symmetric => true, # true → epsilon(1,2)=epsilon(2,1)
        :recombine => true, # true → use combining rules
        :lower => 100.0,
        :guess => 150.0,
        :upper => 200.0
    ),
];

# STEP 4 - PROPERTIES TO FIT
function bubble_point(model_eij::EoSModel,T,x)
    bub = bubble_pressure(model_eij,T,[x,1-x])
    return bub[1]
end

# STEP 5 - ESTIMATOR
path_VLE = "C:/Users/beral/My Drive/DTU/PROJECTS/BS5/CSV/kij/TOACl_CA_1.5.csv"
estimator,objective,initial,upper,lower = Estimation(model_eij,toestimate,[path_VLE]);

method = ECA(;options=Options(iterations=50));

local_params, local_model = optimize(objective, estimator, method);

# PRINT εᵢⱼ
print("### εᵢⱼ of ",species_1,"(i)–",species_2,"(j)")
print("\n")
eij = local_params[1]
@printf("%.4f",eij)
print("\n")

# STEP 6 - kᵢⱼ CALCULATION (Lorentz-Berthelot combining rules)
ei = local_model.params.epsilon[1,1]
ej = local_model.params.epsilon[2,2]
eiej = ei*ej
sqrt_eiej = sqrt(eiej)
eij_sqrt_eiej = eij / sqrt_eiej
kij = 1 - eij_sqrt_eiej

# PRINT kᵢⱼ
print("### kᵢⱼ of ",species_1,"(i)–",species_2,"(j)")
print("\n")
@printf("%.4f",kij)
print("\n")
