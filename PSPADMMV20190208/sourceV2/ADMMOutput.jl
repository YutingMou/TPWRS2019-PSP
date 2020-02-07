## creat output for ADMM in each iteration
function ADMMOutput(NSubproblemsX, niterationsADMM, SlaveX_solutions, r_sol, p_sol,
    primalResidual_rel_vec, dualResidual_rel_vec, supplyR_U, ADMMDateTime)
    # output into files by retrieving solutions
    # define the solution
    totalCost_sol = zeros(Float64, Thours,NScenarios)
    supplyR_sol = zeros(Float64, NOptions, Thours,NScenarios)
    objval_sol = zeros(Float64, NSubproblemsX)
    u_sol = zeros(Int, NGenerators, Thours,NScenarios)
    v_sol = zeros(Int, NGenerators, Thours,NScenarios)
    z_sol = zeros(Int, NGenerators, Thours,NScenarios)
    production_sol = zeros(Float64, NGenerators, Thours,NScenarios)
    productionCost_sol = zeros(Float64, NGenerators, Thours,NScenarios)
    startupCost0_sol = zeros(Float64, NGenerators, Thours,NScenarios)
    miniLoadCost_sol = zeros(Float64, NGenerators, Thours,NScenarios)
    s_sol = zeros(Float64, Thours,NScenarios)
    p_pump_sol = zeros(Float64, Thours,NScenarios)
    p_prod_sol = zeros(Float64, Thours,NScenarios)
    loadShedding_sol = zeros(Float64, Thours,NScenarios)
    renewableShedding_sol = zeros(Float64, Thours,NScenarios)
    # from slave x
    for subproblemIdx = 1:NSubproblemsX
        MonthIdx, ScenarioIdx= returnIDX(subproblemIdx,NSubhoriozns)

        totalCost_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][2]/cSCALE

        for option = OptionsSet
            supplyR_sol[option,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
                SlaveX_solutions[subproblemIdx][3][option,:]/sSCALE[option]   # scale it back
        end

        objval_sol[subproblemIdx] =  SlaveX_solutions[subproblemIdx][4]
        u_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][6]
        v_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][7]
        z_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][8]
        production_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][9]
        productionCost_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][10]
        startupCost0_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][11]
        miniLoadCost_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][12]
        s_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][13]
        p_pump_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][14]
        p_prod_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][15]
        loadShedding_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][16]
        renewableShedding_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][17]
    end

    # subscription profile of each option
    subsProfile = Array{Float64}(NOptions, Thours)
    subsProfileSystem = Array{Float64}(NOptions, Thours)

    for option = OptionsSet
         subsProfile[option, :] = sum(
             u_consumer_initial[consumer, option]*Dl[consumer]*ConsumerProfile[consumer,:]
         for consumer = ConsumersSet)

         subsProfileSystem[option, :] = sum(
             u_consumer_initial[consumer, option]*Dl[consumer]*SystemProfile'
         for consumer = ConsumersSet)
     end


     # calculate realized reliability, need to consider real subscription profile and supply
     # real supply to each option, related to real subscription profile
     residentialD_sol = sum(supplyR_sol,1)
     residentialD_sol = reshape(residentialD_sol, Thours, NScenarios )
     supplyRealized = Array{Float64}(NOptions, Thours, NScenarios)
     for scenario = ScenariosSet
         residentialD_temp = residentialD_sol[:,scenario]
         for option = NOptions:-1:1
             supplyRealized[option,:,scenario] = min.(subsProfile[option, :], residentialD_temp)
             residentialD_temp = residentialD_temp - supplyRealized[option,:,scenario]
         end
     end

     # realized reliablity due to unsynchronized profiles
     r_realized =sum(sum(supplyRealized[:,:,scenario],2)*ScenarioProb[scenario] for scenario = ScenariosSet)./sum(subsProfile,2)
     r_realized = r_realized[:]
     ReliabilityRealizedHourlyPSP = sum(supplyRealized[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet)./subsProfile

     r_realized2 =sum(sum(supplyR_sol[:,:,scenario],2)*ScenarioProb[scenario] for scenario = ScenariosSet)./sum(subsProfileSystem,2)
     r_realized2 = r_realized2[:]



    # output to files of the results from every iterations
    if test_5D == 1
        OutputDir = string(SysDir,"/output/ADMM5D", ADMMDateTime ,"/$(niterationsADMM)")
    else
        OutputDir = string(SysDir,"/output/ADMMFullYear", ADMMDateTime ,"/$(niterationsADMM)")
    end
    mkpath(OutputDir)

    Residual_rel_vec_DF = DataFrame(primalResidual_rel = primalResidual_rel_vec,
    dualResidual_rel = dualResidual_rel_vec)
    writetable(string(OutputDir,"/Residual_rel.csv"),
        Residual_rel_vec_DF; header=true)

    productionCost_sum = sum(sum(productionCost_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    startupCost_sum = sum(sum(startupCost0_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    miniLoadCost_sum = sum(sum(miniLoadCost_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    totalCost_sum = sum(sum(totalCost_sol[:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    totalCostR = totalCost_sum - ICcostTotal        # total cost of residential consumers

    benefitsR = sum(supplyR_sol[option,t,scenario]*(VB[option]+VB[option+1])*0.5*ScenarioProb[scenario] for scenario = ScenariosSet, t = TimeSet, option = OptionsSet)
      # expected residential benefits
    socialWelfareR = benefitsR - totalCostR
    ProfitsRealized = sum(p_sol[option]*Di[option] for option = OptionsSet)*sum(SystemProfile) - totalCostR
    #sum(y_consumer_sol[consumer,option]*Dl[consumer] for option = OptionsSet, consumer = ConsumersSet) *
    #   sum(SystemProfile[t] for t = TimeSet)- totalCostR
    revenueCompany = totalCostR + ProfitsRealized
    netBenefitsR = benefitsR - revenueCompany

   # creat dataframe
    WelfareComparsionDF = DataFrame(
       SocialWelfare = socialWelfareR,
       ConsumerBenefits = benefitsR,
       ConsumerNetBenefits = netBenefitsR,
       CompanyRevenue = revenueCompany,
       CompanyProfits = ProfitsRealized,
       CompanyCosts = totalCostR,
       ProductionCost = productionCost_sum,
       MiniLoadCost = miniLoadCost_sum,
       StartUpCost = startupCost_sum,
       ICCost = ICcostTotal
   )
    writetable(string(OutputDir,"/WelfareComparsionPSP_","$(niterationsADMM)",".csv"),  WelfareComparsionDF; header=true)

    # output subsProfile
    subsProfileDF = DataFrame(subsProfile)
    writetable(string(OutputDir,"/subsProfile.csv"),  subsProfileDF; header=true)

    ReliabilityRealizedHourlyPSPDF = DataFrame(ReliabilityRealizedHourlyPSP)
    writetable(string(OutputDir,"/ReliabilityRealizedHourlyPSP.csv"),  ReliabilityRealizedHourlyPSPDF; header=false)

    # create DataFrame and output as a csv file
#=
    # from u update
    totalCost_UDF = DataFrame(totalCost_U)
    writetable(string(OutputDir,"/totalCost_U.csv"),  totalCost_UDF; header=true)
    supplyR_U_sol = supplyR_U
    supplyR_U_sol = reshape(supplyR_U_sol, NOptions, Thours*NScenarios)
    for option = OptionsSet
        supplyR_U_sol[option,:] = supplyR_U_sol[option,:]/sSCALE[option]
    end
    supplyR_UDF = DataFrame(supplyR_U_sol)
    writetable(string(OutputDir,"/supplyR_U_sol.csv"),  supplyR_UDF; header=true)
=#

    # from x update
    totalCostDF = DataFrame(totalCost_sol)
    writetable(string(OutputDir,"/totalCost.csv"),  totalCostDF; header=true)
    supplyR_sol = reshape(supplyR_sol, NOptions, Thours*NScenarios)
    supplyRDF = DataFrame(supplyR_sol)
    writetable(string(OutputDir,"/supplyR.csv"),  supplyRDF; header=true)
    objvalDF = DataFrame(objval=objval_sol)
    writetable(string(OutputDir,"/objval.csv"),  objvalDF; header=true)
    u_sol = reshape(u_sol, NGenerators, Thours*NScenarios)
    uDF = DataFrame(u_sol)
    writetable(string(OutputDir,"/u.csv"),  uDF; header=true)
    v_sol = reshape(v_sol, NGenerators, Thours*NScenarios)
    vDF = DataFrame(v_sol)
    writetable(string(OutputDir,"/v.csv"),  vDF; header=true)
    z_sol = reshape(z_sol, NGenerators, Thours*NScenarios)
    zDF = DataFrame(z_sol)
    writetable(string(OutputDir,"/z.csv"),  zDF; header=true)
    production_sol = reshape(production_sol, NGenerators, Thours*NScenarios)
    productionDF = DataFrame(production_sol)
    writetable(string(OutputDir,"/production.csv"),  productionDF; header=true)
    productionCost_sol = reshape(productionCost_sol, NGenerators, Thours*NScenarios)
    productionCostDF = DataFrame(productionCost_sol)
    writetable(string(OutputDir,"/productionCost.csv"),  productionCostDF; header=true)
    startupCost0_sol = reshape(startupCost0_sol, NGenerators, Thours*NScenarios)
    startupCost0DF = DataFrame(startupCost0_sol)
    writetable(string(OutputDir,"/startupCost0.csv"),  startupCost0DF; header=true)
    miniLoadCost_sol = reshape(miniLoadCost_sol, NGenerators, Thours*NScenarios)
    miniLoadCostDF = DataFrame(miniLoadCost_sol)
    writetable(string(OutputDir,"/miniLoadCost.csv"),  miniLoadCostDF; header=true)
    s_sol = reshape(s_sol, Thours*NScenarios)
    p_pump_sol = reshape(p_pump_sol, Thours*NScenarios)
    p_prod_sol = reshape(p_prod_sol, Thours*NScenarios)
    PumpedHydroDF = DataFrame(s=s_sol, pump=p_pump_sol, prod=p_prod_sol)
    writetable(string(OutputDir,"/PumpedHydro.csv"),  PumpedHydroDF; header=true)
    loadShedding_sol = reshape(loadShedding_sol, Thours*NScenarios)
    renewableShedding_sol = reshape(renewableShedding_sol, Thours*NScenarios)
    sheddingDF =  DataFrame(loadShedding=loadShedding_sol, renewableShedding=renewableShedding_sol)
    writetable(string(OutputDir,"/shedding.csv"),  sheddingDF; header=true)

    PriceMenuDF2 = DataFrame(Reliability = r_sol,
        ReliabilityRealized = r_realized,
        ReliabilityRealized2 = r_realized2,
        Price = p_sol,
        Subscription = Di)
    writetable(string(OutputDir,"/PriceMenuFinal_","$(niterationsADMM)",".csv"),  PriceMenuDF2; header=true)

end
