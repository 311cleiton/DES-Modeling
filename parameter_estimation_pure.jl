# ESTIMATION OF segment(m),sigma(σ),epsilon(ε)
# https://clapeyronthermo.github.io/Clapeyron.jl/dev/api/estimation/

# STEP 1 - PACKAGES
using Clapeyron # https://github.com/ClapeyronThermo/Clapeyron.jl
using BlackBoxOptim # https://github.com/robertfeldt/BlackBoxOptim.jl
using Statistics
using CSV
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

# STEP 3 - PARAMETERS TO BE FITTED
segmentlower = 3.0
segmentupper = 5.0
segmentguess = median([segmentlower,segmentupper])
sigmalower = 3.0
sigmaupper = 5.0
sigmaguess = median([sigmalower,sigmaupper])
epsilonlower = 200.0
epsilonupper = 400.0
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
test_md = []
function mass_rho(model_pure::EoSModel,p,T)
    md = mass_density(model_pure,p*1.0E6,T)
    append!(test_md,md[1])
    return md[1]
end

# STEP 5 - ESTIMATOR
    # EXPERIMENTAL PATH CSV
    file_path_rho = "C:/Users/beral/My Drive/DTU/PROJECTS/BS5/CSV/PE/TBACl_CA_2.csv"
estimator,objective,initial,upper,lower = Estimation(model_pure,toestimate,[file_path_rho]);

nparams = length(initial)
bounds  = [(lower[i],upper[i]) for i in 1:nparams]

result = BlackBoxOptim.bboptimize(objective; 
        SearchRange = bounds, 
        NumDimensions = nparams,
        MaxSteps=10000,
        PopulationSize = 1000,
        # TraceMode=:silent
        )

local_params = BlackBoxOptim.best_candidate(result); # !!!

# PRINT PARAMETERS
print(species," [segment,sigma,epsilon] =")
print("\n")
@printf("%.4f",local_params[1])
print(",")
@printf("%.4f",local_params[2])
print(",")
@printf("%.2f",local_params[3])
print("\n")
