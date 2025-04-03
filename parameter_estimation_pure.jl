print(@__FILE__)
print("\n")
# ESTIMATION OF segment(m),sigma(σ),epsilon(ε)
# https://github.com/ClapeyronThermo/Clapeyron.jl
# https://clapeyronthermo.github.io/Clapeyron.jl/dev/api/estimation/
# https://github.com/jmejia8/Metaheuristics.jl

# PACKAGES
using Clapeyron
using Metaheuristics
using Statistics
using Printf
using Dates
print(now())
print("\n")

# MODEL
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
choose_species = 9
species_pure = species_all[choose_species]
# model_pure = PCSAFT([species_pure]);
model_pure = SAFTVRMie([species_pure]);

# PARAMETERS TO FIT
segmentlower = 1.0
segmentupper = 7.0
segmentguess = median([segmentlower,segmentupper])
sigmalower = 1.0
sigmaupper = 7.0
sigmaguess = median([sigmalower,sigmaupper])
epsilonlower = 100.00
epsilonupper = 700.00
epsilonguess = median([epsilonlower,epsilonupper])
toestimate = [
    Dict(
        :param => :segment,
        :lower => segmentlower,
        :upper => segmentupper,
        :guess => segmentguess
    ),
    Dict(
        :param => :sigma, # [Å]
        :factor => 1E-10, # convert [Å] to [m]
        :lower => sigmalower,
        :upper => sigmaupper,
        :guess => sigmaguess
    ),
    Dict(
        :param => :epsilon, # [K]
        :lower => epsilonlower,
        :upper => epsilonupper,
        :guess => epsilonguess
    ),
];

# PROPERTY USED FOR ESTIMATION
function mass_rho(model_pure::EoSModel,p,T)
    md = mass_density(model_pure,p*1.0E6,T) # p = [Pa]; T = [K]
    return md[1] # md = [kg/m³]
end

# ESTIMATOR
path_rho = "C:/Users/beral/My Drive/DTU/PROJECTS/BS5/CSV/PE/"*species_pure*".csv"
estimator,objective,initial,upper,lower = Estimation(model_pure,toestimate,[path_rho]);

local_method = ECA(;options=Options(iterations=100))

local_params, local_model = optimize(objective, estimator, local_method);

# PRINT PARAMETERS
print(" ESTIMATED PARAMS:","\n")
@printf("%.4f",local_params[1])
print(",")
@printf("%.4f",local_params[2])
print(",")
@printf("%.2f",local_params[3])
print("\n")

