print(@__FILE__)
print("\n")
# BUBBLE PRESSURE PLOTTING

################
### PACKAGES ###
################
using Clapeyron
using CSV
using DataFrames
using Statistics
using Printf
using PyCall
import PyPlot

@time begin
#############
### INPUT ###
#############
species2 = "TOACl_CA_1.5"
species1 = "CO2_Dufal_etal_2015" # Table 7 [DOI: 10.1080/00268976.2015.1029027]

# cute_species = "[TBA][Cl]:CA (1:2)"
# fig_title = "a) [TBA][Cl]:CA (1:2)"
# cute_species = "[MTA][Cl]:CA (1:2)"
# fig_title = "b) [MTA][Cl]:CA (1:2)"
# cute_species = "[MTA][Br]:CA (1:2)"
# fig_title = "c) [MTA][Br]:CA (1:2)"
# cute_species = "[TOA][Cl]:CA (1:2)"
# fig_title = "d) [TOA][Cl]:CA (1:2)"
# cute_species = "[TOA][Br]:CA (1:2)"
# fig_title = "e) [TOA][Br]:CA (1:2)"
cute_species = "[TOA][Cl]:CA (1:1.5)"
fig_title = "f) [TOA][Cl]:CA (1:1.5)"

####################
### EXPERIMENTAL ###
####################
exp_path = ["C:/Users/clsobe/My Drive/DTU/PROJECTS/BS5/CSV/CO2/TOACl_CA_1.5.csv"]
# exp_path = ["C:/Users/cleiton/My Drive/DTU/PROJECTS/BS5/CSV/CO2/TOACl_CA_1.5.csv"]

exp_data = CSV.read(exp_path,DataFrame)

data_P1 = collect(skipmissing(exp_data[:,1])) # [MPa]
data_T1 = collect(skipmissing(exp_data[:,2])) # [K]
data_x1 = collect(skipmissing(exp_data[:,3])) # [mol/mol]
data_P2 = collect(skipmissing(exp_data[:,4])) # [MPa]
data_T2 = collect(skipmissing(exp_data[:,5])) # [K]
data_x2 = collect(skipmissing(exp_data[:,6])) # [mol/mol]
data_P3 = collect(skipmissing(exp_data[:,7])) # [MPa]
data_T3 = collect(skipmissing(exp_data[:,8])) # [K]
data_x3 = collect(skipmissing(exp_data[:,9])) # [mol/mol]

model_x1 = Clapeyron.FractionVector.(data_x1) # creates a vector [x, 1-x]
model_x2 = Clapeyron.FractionVector.(data_x2)
model_x3 = Clapeyron.FractionVector.(data_x3)

data_x = 0:0.04:0.4
model_x = Clapeyron.FractionVector.(data_x)

################
### MODELING ###
################
model_DES = SAFTVRMie([species1,species2])

# https://clapeyronthermo.github.io/Clapeyron.jl/dev/properties/multi/
# bubble_pressure(model::EoSModel, T, x, method = ChemPotBubblePressure())
# bubble_pressure calculates 4 properties in 4 tuples: [1] bub pressure, [2] V_l, [3] V_v, [4] y
model_VLE1 = bubble_pressure.(model_DES, data_T1[1], model_x)
model_VLE2 = bubble_pressure.(model_DES, data_T2[1], model_x)
model_VLE3 = bubble_pressure.(model_DES, data_T3[1], model_x)

model_P1 = []
model_P2 = []
model_P3 = []
for i in 1:length(model_x)
    append!(model_P1,append!([model_VLE1[i][1]]))
    append!(model_P2,append!([model_VLE2[i][1]]))
    append!(model_P3,append!([model_VLE3[i][1]]))
end

################
### PLOTTING ###
################
PyPlot.clf()
plot_font = "times new roman"
PyPlot.rc("font", family=plot_font)
PyPlot.figure(dpi=311)
# EXPERIMENTAL
PyPlot.plot(data_x1,data_P1,label=string(data_T1[1]," ","K"),linestyle="",marker="o",color="black")
PyPlot.plot(data_x2,data_P2,label=string(data_T2[1]," ","K"),linestyle="",marker="s",color="royalblue")
PyPlot.plot(data_x3,data_P3,label=string(data_T3[1]," ","K"),linestyle="",marker="^",color="forestgreen")
# MODEL
PyPlot.plot(data_x,1e-6*model_P1,label="",linestyle="solid",marker="",color="black")
PyPlot.plot(data_x,1e-6*model_P2,label="",linestyle="solid",marker="",color="royalblue")
PyPlot.plot(data_x,1e-6*model_P3,label="",linestyle="solid",marker="",color="forestgreen")
# 
# PyPlot.legend(loc="best",frameon=false,fontsize=16)
PyPlot.legend(loc=4,frameon=true,fontsize=13)
PyPlot.xlabel("CO₂ mole fraction",fontsize=16)
PyPlot.ylabel("Pressure (MPa)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)
PyPlot.xticks(collect(0.0:0.05:1.0))
xlim_max = 0.3
PyPlot.xlim([0.0,xlim_max])
PyPlot.yticks(collect(0.0:0.5:5.0))
ylim_max = 2.0
PyPlot.ylim([0.0,ylim_max])
PyPlot.text(0.05*xlim_max,0.9*ylim_max, fig_title, fontsize=19)
#
display(PyPlot.gcf())

end#@time
