# ESTIMATION OF segment(m),sigma(σ),epsilon(ε)
# https://github.com/ClapeyronThermo/Clapeyron.jl
# https://clapeyronthermo.github.io/Clapeyron.jl/dev/api/estimation/
# https://github.com/jmejia8/Metaheuristics.jl

# STEP 1 - PACKAGES
using Clapeyron
using Metaheuristics
using Statistics
using Printf
using Dates
print(now())
print("\n")

# STEP 2 - MODEL
species = "TBACl_CA_2"
# model_pure = SAFTVRMie([species]);
model_pure = PCSAFT([species]);
print(model_pure)
print("\n")

# STEP 3 - PARAMETERS TO FIT
segmentlower = 3.0
segmentupper = 5.0
segmentguess = median([segmentlower,segmentupper])
sigmalower = 3.0
sigmaupper = 5.0
sigmaguess = median([sigmalower,sigmaupper])
epsilonlower = 200.00
epsilonupper = 400.00
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

# STEP 4 - PROPERTY USED FOR ESTIMATION
function mass_rho(model_pure::EoSModel,p,T)
    md = mass_density(model_pure,p*1.0E6,T) # p = [Pa]; T = [K]
    return md[1] # md = [kg/m³]
end

# STEP 5 - ESTIMATOR
path_rho = "C:/Users/beral/My Drive/DTU/PROJECTS/BS5/CSV/PE/TBACl_CA_2.csv"
estimator,objective,initial,upper,lower = Estimation(model_pure,toestimate,[path_rho]);

local_method = ECA(;options=Options(iterations=50))

print(local_method)
print("\n")

local_params, local_model = optimize(objective, estimator, local_method);

# PRINT PARAMETERS
print(species," [segment,sigma,epsilon] =")
print("\n")
@printf("%.4f",local_params[1])
print(",")
@printf("%.4f",local_params[2])
print(",")
@printf("%.2f",local_params[3])
print("\n")
