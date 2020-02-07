
# dual decomposition + cutting plane for the one option example

using JuMP, Gurobi, DataFrames
global const SysDir = "C:/Users/Yuting/OneDrive - UCL/Demand Response Business Model/Julia/PSPToyExample/OneOptionExampleADMMNewDatav2"
RHO = ones(Float64, 2)

# bigger rho 100
RHO = RHO*100
# small rho 10
#RHO = RHO*10

TotalCostScale = 5  # scale down totalCost by 5 times
Profits = 6
NOptions = 1
NGenerators = 2
Thours = 2
NConsumers =  10
Dl = 1  #1  average demand
Vmax = 10
Vmin = 0
D = 5   #total amount of subscription
VB = [6, Vmax]
DemandScaling = [1 1]
Vl =  [i for i = 1:10]
MC = [1, 7]
Pmax = [2 0;
        10 0]
Omega = [0.8 0.2]
# parameters
OptionsSet = 1:NOptions
GeneratorsSet = 1:NGenerators
TimeSet = 1:Thours
ConsumersSet = 1:NConsumers

PriceMax = 10 # max price of the whole slice, the value is assigned by us to use McCormick envelope
MaxIterationADMM = 50

rho = zeros(Float64, 2)
r_sol = zeros(Float64, NOptions)
p_sol = zeros(Float64, NOptions)

supplyR_X_sol = zeros(Float64, NOptions, Thours)
supplyR_X_sol_vec = zeros(Float64, Thours, MaxIterationADMM)  # only one option

totalCost_X_sol = zeros(Float64, Thours)
totalCost_X_sol_vec = zeros(Float64, Thours, MaxIterationADMM)  # only record the first scenario


production_X_sol_vec = zeros(Float64, NGenerators, MaxIterationADMM)


supplyR_Z_sol = zeros(Float64, NOptions, Thours)
supplyR_Z_sol_vec = zeros(Float64, Thours, MaxIterationADMM)  # only record the first scenario

totalCost_Z_sol = zeros(Float64, Thours)
totalCost_Z_sol_vec = zeros(Float64, Thours, MaxIterationADMM)  # only record the first scenario


supplyR_U_sol = zeros(Float64, NOptions, Thours)
supplyR_U_sol_vec = zeros(Float64, Thours, MaxIterationADMM)  # only record the first scenario

totalCost_U_sol = zeros(Float64, Thours)
totalCost_U_sol_vec = zeros(Float64, Thours, MaxIterationADMM)  # only record the first scenario

r_sol_vec =  zeros(Float64, NOptions, MaxIterationADMM)
p_sol_vec = zeros(Float64, NOptions, MaxIterationADMM)

# used to calculate the dual residual

AugmentedLagrangianValue =  zeros(Float64, MaxIterationADMM)
CMUpdateValue = zeros(Float64, MaxIterationADMM)
PMUpdateValue =  zeros(Float64, MaxIterationADMM)

# residuals in ADMM iterations
#residual_primal = zeros(Float64, MaxIterationADMM)
#residual_dual = zeros(Float64, MaxIterationADMM)
# relative residuals in ADMM iterations
#residual_primal_rel = zeros(Float64, MaxIterationADMM)
#residual_dual_rel = zeros(Float64, MaxIterationADMM)


#include("PRM.jl")
#include("PRMProducer.jl")
include("ADMMPMOneOptionNewDataV2.jl")
include("ADMMCMOneOptionNewDataV2.jl")
println("###########################################")

niterationsADMM = 1
while niterationsADMM <= MaxIterationADMM
    println("##################### niterationsADMM: ", niterationsADMM)
    #println("dual_vb_k: ", dual_vb_k)
    if niterationsADMM == 1
        println("X-Update")
        flush(STDOUT)
        ## run  PM
        objvalue_PM, supplyR_X_sol, totalCost_X_sol,production_X_sol =
            ADMMPM(niterationsADMM, rho, supplyR_Z_sol, totalCost_Z_sol, supplyR_U_sol, totalCost_U_sol)
        rho = RHO
        ## run CM
        objvalue_CM, r_sol, p_sol,supplyR_Z_sol,totalCost_Z_sol =
            ADMMCM(rho, supplyR_X_sol, totalCost_X_sol, supplyR_U_sol, totalCost_U_sol)
        ## update dual
        include("ADMMDualUpdateOneOptionNewDataV2.jl")
    else
        rho = RHO
        objvalue_PM, supplyR_X_sol, totalCost_X_sol,production_X_sol =
            ADMMPM(niterationsADMM, rho, supplyR_Z_sol, totalCost_Z_sol, supplyR_U_sol, totalCost_U_sol)
        objvalue_CM, r_sol, p_sol,supplyR_Z_sol,totalCost_Z_sol =
            ADMMCM(rho, supplyR_X_sol, totalCost_X_sol, supplyR_U_sol, totalCost_U_sol)
        include("ADMMDualUpdateOneOptionNewDataV2.jl")
    end



    println("objvalue_PM: ", objvalue_PM)
    println("production_X_sol: ", production_X_sol)
    println("supplyR_X_sol: ", supplyR_X_sol)
    println("totalCost_X_sol: ", totalCost_X_sol)
    println("supplyR_Z_sol: ", supplyR_Z_sol)
    println("totalCost_Z_sol: ", totalCost_Z_sol)

    println("totalCost_U_sol: ", totalCost_U_sol)
    println("supplyR_U_sol: ", supplyR_U_sol)

    production_X_sol_vec[:,niterationsADMM]  =    production_X_sol[:,1]
    supplyR_X_sol_vec[:,niterationsADMM] = supplyR_X_sol[1,:]
    totalCost_X_sol_vec[:,niterationsADMM] = totalCost_X_sol*TotalCostScale

    supplyR_Z_sol_vec[:,niterationsADMM] = supplyR_Z_sol[1,:]
    totalCost_Z_sol_vec[:,niterationsADMM] = totalCost_Z_sol*TotalCostScale

    supplyR_U_sol_vec[:,niterationsADMM] = supplyR_U_sol[1,:]
    totalCost_U_sol_vec[:,niterationsADMM] = totalCost_U_sol*TotalCostScale


    r_sol_vec[:, niterationsADMM] = r_sol
    p_sol_vec[:, niterationsADMM] = p_sol



    AugmentedLagrangianValue[niterationsADMM] = objvalue_CM + objvalue_PM
    CMUpdateValue[niterationsADMM] = objvalue_CM
    PMUpdateValue[niterationsADMM] = objvalue_PM



    #println("AugmentedLagrangianValue: ", AugmentedLagrangianValue)
    #println("CMUpdateValue: ", CMUpdateValue)
    #println("PMUpdateValue: ", PMUpdateValue)
    #println("d1_vec: ", d1_vec)
    #println("r_sol_vec: ", r_sol_vec)
    #println("p_sol_vec: ", p_sol_vec)

    niterationsADMM = niterationsADMM +1

end

DDDateTime = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")

println("production_X_sol_vec: ", production_X_sol_vec)
#println("r_sol_vec: ", r_sol_vec)
#println("p_sol_vec: ", p_sol_vec)


#=
AugmentedLagrangianValueDF = DataFrame(AugmentedLagrangianValue=AugmentedLagrangianValue)
writetable(string(SysDir,"/AugmentedLagrangianValue", DDDateTime,".csv"), AugmentedLagrangianValueDF; header=true)

objvalue_master_vecDF = DataFrame(objvalue_master_vec=objvalue_master_vec)
writetable(string(SysDir,"/objvalue_master_vec", DDDateTime,".csv"), objvalue_master_vecDF; header=true)



PMUpdateValueDF = DataFrame(PMUpdateValue=PMUpdateValue)
writetable(string(SysDir,"/PMUpdateValue", DDDateTime,".csv"), PMUpdateValueDF; header=true)

CMUpdateValueDF = DataFrame(CMUpdateValue=CMUpdateValue)
writetable(string(SysDir,"/CMUpdateValue", DDDateTime,".csv"),CMUpdateValueDF; header=true)
=#
if false
r_solDF = DataFrame(r_sol_vec)
writetable(string(SysDir,"/r_sol_vec", DDDateTime,".csv"), r_solDF; header=true)

p_solDF = DataFrame(p_sol_vec)
writetable(string(SysDir,"/p_sol_vec", DDDateTime,".csv"), p_solDF; header=true)

production_X_sol_vecDF = DataFrame(production_X_sol_vec)
writetable(string(SysDir,"/production_X_sol_vec", DDDateTime,".csv"), production_X_sol_vecDF; header=true)

#residual_primal_relDF = DataFrame(residual_primal_rel)
#writetable(string(SysDir,"/residual_primal_rel", DDDateTime,".csv"), residual_primal_relDF; header=true)

#residual_dual_relDF = DataFrame(residual_dual_rel)
#writetable(string(SysDir,"/residual_dual_rel", DDDateTime,".csv"), residual_dual_relDF; header=true)
end
if false

production_X_sol_vecDF = DataFrame(production_X_sol_vec)
writetable(string(SysDir,"/production_X_sol_vec", DDDateTime,".csv"), production_X_sol_vecDF; header=true)

supplyR_X_sol_vecDF = DataFrame(supplyR_X_sol_vec)
writetable(string(SysDir,"/supplyR_X_sol_vec.csv", DDDateTime,".csv"), supplyR_X_sol_vecDF; header=true)
totalCost_X_sol_vecDF = DataFrame(totalCost_X_sol_vec)
writetable(string(SysDir,"/totalCost_X_sol_vec.csv", DDDateTime,".csv"), totalCost_X_sol_vecDF; header=true)

supplyR_Z_sol_vecDF = DataFrame(supplyR_Z_sol_vec)
writetable(string(SysDir,"/supplyR_Z_sol_vec.csv", DDDateTime,".csv"), supplyR_Z_sol_vecDF; header=true)
totalCost_Z_sol_vecDF = DataFrame(totalCost_Z_sol_vec)
writetable(string(SysDir,"/totalCost_Z_sol_vec.csv", DDDateTime,".csv"), totalCost_Z_sol_vecDF; header=true)

supplyR_U_sol_vecDF = DataFrame(supplyR_U_sol_vec)
writetable(string(SysDir,"/supplyR_U_sol_vec.csv", DDDateTime,".csv"), supplyR_U_sol_vecDF; header=true)
totalCost_U_sol_vecDF = DataFrame(totalCost_U_sol_vec)
writetable(string(SysDir,"/totalCost_U_sol_vec.csv", DDDateTime,".csv"), totalCost_U_sol_vecDF; header=true)
end
