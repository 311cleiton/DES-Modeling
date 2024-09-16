print(@__FILE__)
print("\n")
# Reduced Helmholtz energy (dadT and d²adT² plot)
# Comparison of each contribution (at constant V,N)

################
### PACKAGES ###
################
using Clapeyron
using Dierckx
using CSV
using DataFrames
using Statistics
using Printf
using LaTeXStrings
using PyCall
import PyPlot

#############
### INPUT ###
#############
species_des = [
"TPPPBr_Ph_6"   # 1
"TPPPBr_Ph_4"   # 2
"TPPPBr_TEG_16" # 3
"TPPPBr_TEG_10" # 4
"TPPPBr_TEG_4"  # 5
"TBACl_CA_2"    # 6
"MTACl_CA_2"    # 7
"MTABr_CA_2"    # 8
"TOACl_CA_2"    # 9
"TOABr_CA_2"    # 10
"TOACl_CA_1.5"  # 11
]
species_evaluated = species_des[11]
species_1 = species_evaluated
fig_title_1 = species_1

################
### MODELING ###
################
# DEFINE MODEL
model_1 = SAFTVRMie([species_1])

# MODEL INPUTS
fixed_T = 298.15 # temperature [K]
input_P = 0.1E+06 # pressure [Pa]
input_T = collect(290:1:340) # temperature [K]
input_v = volume.(model_1, input_P, fixed_T; phase=:liquid) # volume [m³]

# MODEL OUTPUTS
output_a_hs = Clapeyron.a_hs.(model_1, input_v, input_T, [1.0])
output_a_disp = Clapeyron.a_disp.(model_1, input_v, input_T, [1.0])
output_a_mono = Clapeyron.a_mono.(model_1, input_v, input_T, [1.0])
output_a_chain = Clapeyron.a_chain.(model_1, input_v, input_T, [1.0])
output_a_assoc = Clapeyron.a_assoc.(model_1, input_v, input_T, [1.0])
output_a_res = Clapeyron.a_res.(model_1, input_v, input_T, [1.0])
output_a_test = output_a_res
output_a_dispchain = Clapeyron.a_dispchain.(model_1, input_v, input_T, [1.0])
output_a_residual = output_a_hs + output_a_disp + output_a_chain + output_a_assoc

# https://github.com/kbarbary/Dierckx.jl
spl_res = Spline1D(input_T,output_a_res,k=2)
spl_hs = Spline1D(input_T,output_a_hs,k=2)
spl_disp = Spline1D(input_T,output_a_disp,k=2)
spl_chain = Spline1D(input_T,output_a_chain,k=2)
spl_assoc = Spline1D(input_T,output_a_assoc,k=2)

d_res = derivative(spl_res,input_T)
d_hs = derivative(spl_hs,input_T)
d_disp = derivative(spl_disp,input_T)
d_chain = derivative(spl_chain,input_T)
d_assoc = derivative(spl_assoc,input_T)

spl_d_res = Spline1D(input_T,d_res,k=2)
spl_d_hs = Spline1D(input_T,d_hs,k=2)
spl_d_disp = Spline1D(input_T,d_disp,k=2)
spl_d_chain = Spline1D(input_T,d_chain,k=2)
spl_d_assoc = Spline1D(input_T,d_assoc,k=2)

d2_res = derivative(spl_d_res,input_T)
d2_hs = derivative(spl_d_hs,input_T)
d2_disp = derivative(spl_d_disp,input_T)
d2_chain = derivative(spl_d_chain,input_T)
d2_assoc = derivative(spl_d_assoc,input_T)

####################
### PLOTTING dadT###
####################
PyPlot.clf()
plot_font = "times new roman"
PyPlot.rc("font", family=plot_font)
PyPlot.figure(dpi=311)
#
label_res = "a"*L"\rm{_r}"*L"\rm{_e}"*L"\rm{_s}"
label_hs = "a"*L"\rm{_H}"*L"\rm{_S}"
label_disp = "a"*L"\rm{_d}"*L"\rm{_i}"*L"\rm{_s}"*L"\rm{_p}"
label_mono = "a"*L"\rm{_m}"*L"\rm{_o}"*L"\rm{_n}"*L"\rm{_o}"
label_chain = "a"*L"\rm{_c}"*L"\rm{_h}"*L"\rm{_a}"*L"\rm{_i}"*L"\rm{_n}"
label_assoc = "a"*L"\rm{_a}"*L"\rm{_s}"*L"\rm{_s}"*L"\rm{_o}"*L"\rm{_c}"
#
PyPlot.plot(input_T,d_res,label=label_res,linestyle="solid",marker="",color="black")
PyPlot.plot(input_T,d_hs,label=label_hs,linestyle="dotted",marker="",color="royalblue")
PyPlot.plot(input_T,d_disp,label=label_disp,linestyle="dashed",marker="",color="hotpink")
PyPlot.plot(input_T,d_chain,label=label_chain,linestyle="dashdot",marker="",color="forestgreen")
PyPlot.plot(input_T,d_assoc,label=label_assoc,linestyle=(5, (10, 3)),marker="",color="orange")
#
PyPlot.legend(loc=6,frameon=true,fontsize=13,ncol=2)
PyPlot.xlabel("Temperature (K)",fontsize=16)
PyPlot.ylabel("∂a/∂T (K⁻¹)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)
PyPlot.title(fig_title_1,fontsize=19)
display(PyPlot.gcf())

#######################
### PLOTTING  d²adT²###
#######################
PyPlot.clf()
plot_font = "times new roman"
PyPlot.rc("font", family=plot_font)
PyPlot.figure(dpi=311)
#
label_res = "a"*L"\rm{_r}"*L"\rm{_e}"*L"\rm{_s}"
label_hs = "a"*L"\rm{_H}"*L"\rm{_S}"
label_disp = "a"*L"\rm{_d}"*L"\rm{_i}"*L"\rm{_s}"*L"\rm{_p}"
label_mono = "a"*L"\rm{_m}"*L"\rm{_o}"*L"\rm{_n}"*L"\rm{_o}"
label_chain = "a"*L"\rm{_c}"*L"\rm{_h}"*L"\rm{_a}"*L"\rm{_i}"*L"\rm{_n}"
label_assoc = "a"*L"\rm{_a}"*L"\rm{_s}"*L"\rm{_s}"*L"\rm{_o}"*L"\rm{_c}"
#
PyPlot.plot(input_T,d2_res*1000,label=label_res,linestyle="solid",marker="",color="black")
PyPlot.plot(input_T,d2_hs*1000,label=label_hs,linestyle="dotted",marker="",color="royalblue")
PyPlot.plot(input_T,d2_disp*1000,label=label_disp,linestyle="dashed",marker="",color="hotpink")
PyPlot.plot(input_T,d2_chain*1000,label=label_chain,linestyle="dashdot",marker="",color="forestgreen")
PyPlot.plot(input_T,d2_assoc*1000,label=label_assoc,linestyle=(5, (10, 3)),marker="",color="orange")
#
PyPlot.legend(loc=6,frameon=true,fontsize=13,ncol=2)
PyPlot.xlabel("Temperature (K)",fontsize=16)
PyPlot.ylabel("∂²a/∂T² (1∙10³ K⁻²)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)
PyPlot.title(fig_title_1,fontsize=19)
display(PyPlot.gcf())
