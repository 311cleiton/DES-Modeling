print(@__FILE__)
print("\n")
# RHO
# PLOTTING

################
### PACKAGES ###
################
using Clapeyron
# using BlackBoxOptim
using CSV
using DataFrames
using Statistics
using Printf
using PyCall
import PyPlot

#############
### INPUT ###
#############
species_DES = [
   "TBACl_CA_2"
   "MTACl_CA_2"
   "MTABr_CA_2"
   "TOACl_CA_2"
   "TOABr_CA_2"
   "TOACl_CA_1.5"
]
length_DES = length(species_DES)

cute_species = [
    "[TBA][Cl]:CA (1:2)"
    "[MTA][Cl]:CA (1:2)"
    "[MTA][Br]:CA (1:2)"
    "[TOA][Cl]:CA (1:2)"
    "[TOA][Br]:CA (1:2)"
    "[TOA][Cl]:CA (1:1.5)"
]
# fig_title = species_des

################
### MODELING ###
################
username = ENV["USERNAME"]
# exp_path = ["C:/Users/"*username*"/My Drive/PROJECTS/DTU/BS5/CSV/RHO/TOACl_CA_1.5.csv"]
# exp_path = ["C:/Users/clsobe/My Drive/DTU/PROJECTS/BS5/CSV/RHO/TOACl_CA_1.5.csv"]
exp_path = [] 
exp_data = [] # 1 - rho [kg/m³] | 2 – pressure [MPa] | 3 – temperature [K] 
data_rho = [] # specific density [kg/m³]
# data_P = [] # pressure [MPa]
model_P = 0.1*1E6 # Pa
data_T = [] # K
model_T = collect(280:10:330) # temperature [K]
model_DES = []
model_rho = [] # specific density [kg/m³]
for i in 1:length_DES
    # READ EXPERIMENTAL DATA
    # append!(exp_path,["C:/Users/cleiton/My Drive/DTU/PROJECTS/BS5/CSV/RHO/"*species_DES[i]*".csv"])
    append!(exp_path,["C:/Users/"*username*"/My Drive/DTU/PROJECTS/BS5/CSV/RHO/"*species_DES[i]*".csv"])
    append!(exp_data,[CSV.read(exp_path[i],DataFrame)])
    # COLLECT AND ORGANIZE EXPERIMENTAL DATA
    append!(data_rho,[collect(skipmissing(exp_data[i][:,1]))])
    # append!(data_P,[collect(skipmissing(exp_data[i][:,2]))])
    append!(data_T,[collect(skipmissing(exp_data[i][:,3]))])

    # GENERATE MODEL
    append!(model_DES,[SAFTVRMie([species_DES[i]])])
    # CALCULATE PROPERTY
    append!(model_rho,[mass_density.(model_DES[i],model_P,model_T; phase=:liquid)])
end

################
### PLOTTING ###
################
PyPlot.clf()
plot_font = "times new roman"
PyPlot.rc("font", family=plot_font)
# PyPlot.figure(figsize=(6,6), dpi = 311)
PyPlot.figure(dpi=311)
# EXPERIMENTAL
# PyPlot.plot(data_T,data_rho,label=species_des,linestyle="",marker="o",color="black")
PyPlot.plot(data_T[1],data_rho[1],label=cute_species[1],linestyle="",marker="o",color="black")
PyPlot.plot(data_T[2],data_rho[2],label=cute_species[2],linestyle="",marker="s",color="royalblue")
PyPlot.plot(data_T[3],data_rho[3],label=cute_species[3],linestyle="",marker="^",color="forestgreen")
PyPlot.plot(data_T[4],data_rho[4],label=cute_species[4],linestyle="",marker="D",color="orange")
PyPlot.plot(data_T[5],data_rho[5],label=cute_species[5],linestyle="",marker="*",color="red")
PyPlot.plot(data_T[6],data_rho[6],label=cute_species[6],linestyle="",marker="x",color="indigo")
# MODEL
PyPlot.plot(model_T,model_rho[1],label="",linestyle="solid",marker="",color="black")
PyPlot.plot(model_T,model_rho[2],label="",linestyle="solid",marker="",color="royalblue")
PyPlot.plot(model_T,model_rho[3],label="",linestyle="solid",marker="",color="forestgreen")
PyPlot.plot(model_T,model_rho[4],label="",linestyle="solid",marker="",color="orange")
PyPlot.plot(model_T,model_rho[5],label="",linestyle="solid",marker="",color="red")
PyPlot.plot(model_T,model_rho[6],label="",linestyle="solid",marker="",color="indigo")
# 
PyPlot.legend(loc=2,frameon=true,fontsize=13,ncol=2)
# PyPlot.legend(loc=(1.01,0.25),frameon=true,fontsize=16)
PyPlot.xlabel("Temperature (K)",fontsize=16)
PyPlot.ylabel("Density (kg/m³)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)
xlim_min = 280
xlim_max = 330
PyPlot.xticks(collect(xlim_min:10:xlim_max))
PyPlot.xlim([xlim_min,xlim_max])
# ylim_min = 860
# ylim_max = 960
ylim_min = 850
ylim_max = 1000
# PyPlot.yticks(collect(ylim_min:20:ylim_max))
PyPlot.yticks(collect(ylim_min:25:ylim_max))
PyPlot.ylim([ylim_min,ylim_max])
# PyPlot.text(property_2B[150],700, "298.15 K, 313.15 K, 323.15 K, 333.15 K, 343.15 K, 353.15 K, 363.16 K", fontsize=14)
#
display(PyPlot.gcf())