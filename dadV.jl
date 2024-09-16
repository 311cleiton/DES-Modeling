print(@__FILE__)
print("\n")
# Reduced Helmholtz energy da/dV (at constant T,N) (plot with who)
# https://clapeyronthermo.github.io/Clapeyron.jl/dev/properties/basic/
# https://clapeyronthermo.github.io/Clapeyron.jl/dev/properties/bulk/

################
### PACKAGES ###
################
using Clapeyron
using Dierckx
using LaTeXStrings
using PyCall
import PyPlot
using Statistics

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
species_evaluated = species_des[6]
species_1 = species_evaluated
fig_title_1 = "[TBA][Cl]:CA (1:2)"

################
### MODELING ###
################
model_1 = SAFTVRMie([species_1])

# Model inputs
input_P = collect(0.1E+06:1.0E+06:150.1E+06) # pressure [Pa]
fixed_P = 0.1E+06 # pressure [Pa]
input_T = collect(270:1:330) # temperature [K]
fixed_T = 298.15 # temperature [K]
fixed_z = [1.0] # mole number [mol]

input_v = volume.(model_1, input_P, fixed_T; phase=:liquid) # volume [m³]
# mass_density(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true)
model_rho = mass_density.(model_1, input_P, fixed_T; phase=:liquid) # specific density [kg/m³]

r_v = reverse(input_v)
r_rho = reverse(model_rho)

# a_res(model::EoSModel, V, T, z,args...) (Reduced residual Helmholtz free energy)
output_a_res = Clapeyron.a_res.(model_1, r_v, fixed_T, fixed_z)
output_a_hs = Clapeyron.a_hs.(model_1, r_v, fixed_T, fixed_z)
output_a_disp = Clapeyron.a_disp.(model_1, r_v, fixed_T, fixed_z)
output_a_chain = Clapeyron.a_chain.(model_1, r_v, fixed_T, fixed_z)
output_a_assoc = Clapeyron.a_assoc.(model_1, r_v, fixed_T, fixed_z)

# https://github.com/kbarbary/Dierckx.jl
spl_res = Spline1D(r_v,output_a_res,k=1)
spl_hs = Spline1D(r_v,output_a_hs,k=1)
spl_disp = Spline1D(r_v,output_a_disp,k=1)
spl_chain = Spline1D(r_v,output_a_chain,k=1)
spl_assoc = Spline1D(r_v,output_a_assoc,k=1)

d_res = derivative(spl_res,r_v)
d_hs = derivative(spl_hs,r_v)
d_disp = derivative(spl_disp,r_v)
d_chain = derivative(spl_chain,r_v)
d_assoc = derivative(spl_assoc,r_v)

d_res = d_res*(-1e-6)
d_hs = d_hs*(-1e-6)
d_disp = d_disp*(-1e-6)
d_chain = d_chain*(-1e-6)
d_assoc = d_assoc*(-1e-6)

abs_d_res = abs.(d_res)
abs_d_hs = abs.(d_hs)
abs_d_disp = abs.(d_disp)
abs_d_chain = abs.(d_chain)
abs_d_assoc = abs.(d_assoc)

d_overall = abs_d_hs + abs_d_disp + abs_d_chain + abs_d_assoc

abs_assoc_VS_overall = 100 .* (abs_d_assoc ./ d_overall)

#############################
### PLOTTING da/dV vs rho ###
#############################
PyPlot.clf()
plot_font = "times new roman"
PyPlot.rc("font", family=plot_font)
PyPlot.figure(dpi=311)
#
label_res = "a"*L"\rm{_r}"*L"\rm{_e}"*L"\rm{_s}"
label_hs = "a"*L"\rm{_H}"*L"\rm{_S}"
label_disp = "a"*L"\rm{_d}"*L"\rm{_i}"*L"\rm{_s}"*L"\rm{_p}"
label_chain = "a"*L"\rm{_c}"*L"\rm{_h}"*L"\rm{_a}"*L"\rm{_i}"*L"\rm{_n}"
label_assoc = "a"*L"\rm{_a}"*L"\rm{_s}"*L"\rm{_s}"*L"\rm{_o}"*L"\rm{_c}"
#
PyPlot.plot(r_rho,d_hs,label=label_hs,linestyle="dashed",marker="",color="royalblue")
PyPlot.plot(r_rho,d_disp,label=label_disp,linestyle="dotted",marker="",color="forestgreen")
PyPlot.plot(r_rho,d_chain,label=label_chain,linestyle="dashdot",marker="",color="orange")
PyPlot.plot(r_rho,d_assoc,label=label_assoc,linestyle=(5, (10, 3)),marker="",color="red")
PyPlot.plot(r_rho,d_res,label=label_res,linestyle="solid",marker="",color="black")
#
PyPlot.legend(loc=(0.03,0.5),frameon=true,fontsize=16,ncol=3)
PyPlot.xlabel("Mass density (kg/m³)",fontsize=16)
PyPlot.ylabel("–∂a/∂V (1∙10⁻⁶ m⁻³)",fontsize=16)
PyPlot.xticks(fontsize=13)
PyPlot.yticks(fontsize=13)
xlim_min = 930
xlim_max = 970
PyPlot.xlim([xlim_min,xlim_max])
ylim_min = -0.1
ylim_max = 0.2
PyPlot.ylim([ylim_min,ylim_max])
PyPlot.text(1.001*xlim_min,0.90*ylim_max, fig_title_1, fontsize=19, style="normal",font="times new roman")
display(PyPlot.gcf())