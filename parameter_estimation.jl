print(@__FILE__)
print("\n")
@time begin
# PARAMETER ESTIMATION

# EXPERIMENTAL PATH CSV
full_file_path_rho = "C:/Users/clsobe/My Drive/DTU/PROJECTS/BS5/CSV/PE/TOABr_CA_2.csv"
# full_file_path_rho = "C:/Users/cleiton/My Drive/DTU/PROJECTS/BS5/CSV/PE/TOABr_CA_2.csv"

# PACKAGES
using Clapeyron
using BlackBoxOptim
using Statistics
using Printf

# GENERATE MODEL YOU WISH TO PARAMETRIZE
species = "TOABr_CA_2"
model = SAFTVRMie([species]);

# ESTIMATION PARAMETERS
segmentlower = 4.5543
segmentupper = 4.5543
segmentguess = median([segmentlower,segmentupper])
sigmalower = 4.6441
sigmaupper = 4.7691
sigmaguess = median([sigmalower,sigmaupper])
epsilonlower = 400.00
epsilonupper = 400.00
epsilonguess = median([epsilonlower,epsilonupper])

# ESTIMATION FRAMEWORK
toestimate = [
    Dict(
        :param => :segment,
        # :indices => 1,
        # :cross_assoc => true,
        :lower => segmentlower,
        :upper => segmentupper,
        :guess => segmentguess
    ),
    Dict(
        :param => :sigma, # [Å]
        # :indices => 1,
        # :cross_assoc => true,
        :factor => 1E-10, # convert [Å] to SI unit [m]
        :lower => sigmalower,
        :upper => sigmaupper,
        :guess => sigmaguess
    ),
    # Dict(
    #     :param => :epsilon, # [K]
    #     # :indices => 1,
    #     # :cross_assoc => true,
    #     :lower => epsilonlower,
    #     :upper => epsilonupper,
    #     :guess => epsilonguess
    # ),
    # Dict(
    #     :param => :lambda_r,
    #     # :indices => 1,
    #     # :cross_assoc => true,
    #     :lower => 17.0,
    #     :upper => 25.0,
    #     :guess => 21.0
    # ),
];

# PROPERTY FUNCTION YOU WISH TO ESTIMATE
function mass_rho(model::EoSModel,p,T)
    md = mass_density(model,p*1.0E6,T)
    return md[1]
end

# ESTIMATOR FUNCTION
estimator,objective,initial,upper,lower = Estimation(model,toestimate,[full_file_path_rho]);

# OPTMIZATION
nparams = length(initial)
bounds  = [(lower[i],upper[i]) for i in 1:nparams]

result = BlackBoxOptim.bboptimize(objective; 
        SearchRange = bounds, 
        NumDimensions = nparams,
        # MaxSteps=10000,
        MaxSteps=20000,
        # MaxSteps=42000,
        PopulationSize = 1000,
        # PopulationSize = 2000,
        TraceMode=:silent)

local_params = BlackBoxOptim.best_candidate(result);

# PRINT PARAMETERS
print(species," [segment,sigma,epsilon] =")
print("\n")
# @printf("%.2f",params[1],",","%.2f",params[2],",","%.2f",params[3])
@printf("%.4f",local_params[1])
print(",")
@printf("%.4f",local_params[2])
# print(",")
# @printf("%.2f",local_params[3])
print("\n")

end#@time