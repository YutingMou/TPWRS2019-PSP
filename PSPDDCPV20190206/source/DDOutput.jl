## creat output for DD in each iteration
function DDOutput(NSubproblemsX, niterationsDD, SlaveX_solutions,
    r_sol, r_sol_vec, p_sol, p_sol_vec, r_realized_vec, ProfitsRealized_vec, residual_profits_vec,
    dual_vb_vec, dual_reliability_vec, dual_profits_vec,
    u_consumer_sol, y_consumer_sol, UB_vec, DDDateTime)

    println("Output data in each DD iteration")
    flush(STDOUT)

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
    #u_partition_chosen_sol = zeros(Int, NSubproblemsX)
    # from slave x
    for subproblemIdx = 1:NSubproblemsX
        MonthIdx, ScenarioIdx= returnIDX(subproblemIdx,NSubhoriozns)

        totalCost_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][2]

        for option = OptionsSet
            supplyR_sol[option,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
                SlaveX_solutions[subproblemIdx][3][option,:]/sSCALE[option]   # scale it back
        end

        u_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][4]
        v_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][5]
        z_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][6]
        production_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][7]
        productionCost_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][8]
        startupCost0_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][9]
        miniLoadCost_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][10]
        s_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][11]
        p_pump_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][12]
        p_prod_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][13]
        loadShedding_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][14]
        renewableShedding_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][15]
        #u_partition_chosen_sol[subproblemIdx] =  SlaveX_solutions[subproblemIdx][16]
    end


    # consumers average valuation is actually different at each hour if profiles are not the same,
    # because the demand prfiles would change
    Vit =  zeros(Float64, NOptions, Thours) # average valuation of each option at each hour
    for t = TimeSet
        for option = OptionsSet
            Vit[option, t] =
            sum(u_consumer_sol[consumer, option]*Dl[consumer]*Vl[consumer]*ConsumerProfile[consumer, t] for consumer = ConsumersSet)/sum(u_consumer_sol[consumer, option]*Dl[consumer]*ConsumerProfile[consumer, t] for consumer = ConsumersSet)
        end
    end


    # subscription profile of each option
    subsProfile = Array{Float64}(NOptions, Thours)
    subsProfileSystem = Array{Float64}(NOptions, Thours)

    for option = OptionsSet
        subsProfile[option, :] = sum(
            u_consumer_sol[consumer, option]*Dl[consumer]*ConsumerProfile[consumer,:]
        for consumer = ConsumersSet)

        subsProfileSystem[option, :] = sum(
            u_consumer_sol[consumer, option]*Dl[consumer]*SystemProfile'
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
        OutputDir = string(SysDir,"/output/DD5D", DDDateTime ,"/$(niterationsDD)")
    else
        OutputDir = string(SysDir,"/output/DDFullYear", DDDateTime ,"/$(niterationsDD)")
    end


    #OutputDir = string(SysDir,"/output/MIP")
    mkpath(OutputDir)

    UB_vec_DF = DataFrame(UB_vec[1:niterationsDD])
    writetable(string(OutputDir,"/UB_vec_","$(niterationsDD)",".csv"),
        UB_vec_DF; header=true)

    r_realized_vec_DF = DataFrame(r_realized_vec[:, 1:niterationsDD])
    writetable(string(OutputDir,"/r_realized_","$(niterationsDD)",".csv"),
        r_realized_vec_DF; header=true)

    r_sol_vec_DF = DataFrame(r_sol_vec[:, 1:niterationsDD])
    writetable(string(OutputDir,"/r_sol_","$(niterationsDD)",".csv"),
        r_sol_vec_DF; header=true)
    p_sol_vec_DF = DataFrame(p_sol_vec[:, 1:niterationsDD])
    writetable(string(OutputDir,"/p_sol_","$(niterationsDD)",".csv"),
        p_sol_vec_DF; header=true)

    ProfitsRealized_vec_DF = DataFrame(ProfitsRealized_vec[1:niterationsDD])
    writetable(string(OutputDir,"/ProfitsRealized_","$(niterationsDD)",".csv"),
        ProfitsRealized_vec_DF; header=true)

    residual_profits_vec_DF = DataFrame(residual_profits_vec[1:niterationsDD])
    writetable(string(OutputDir,"/residual_profits_vec_","$(niterationsDD)",".csv"),
            residual_profits_vec_DF; header=true)

    dual_vb_vec_DF = DataFrame(dual_vb_vec[:, 1:niterationsDD])
    writetable(string(OutputDir,"/dual_vb_vec_","$(niterationsDD)",".csv"),
        dual_vb_vec_DF; header=true)

    dual_reliability_vec_DF = DataFrame(dual_reliability_vec[:, 1:niterationsDD])
    writetable(string(OutputDir,"/dual_reliability_vec_","$(niterationsDD)",".csv"),
        dual_reliability_vec_DF; header=true)

    dual_profits_vec_DF = DataFrame(dual_profits_vec[1:niterationsDD])
    writetable(string(OutputDir,"/dual_profits_vec_","$(niterationsDD)",".csv"),
            dual_profits_vec_DF; header=true)


#    Residual_rel_DD_DF = DataFrame(primalResidual_rel = primalResidual_rel_vec[1:niterationsDD],
#        alpha = Alpha_DD[1:niterationsDD], dual_profits = dual_profits_vec[1:niterationsDD])
#    writetable(string(OutputDir,"/Residual_rel_","$(niterationsDD)",".csv"),
#        Residual_rel_DD_DF; header=true)
if true
    productionCost_sum = sum(sum(productionCost_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    startupCost_sum = sum(sum(startupCost0_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    miniLoadCost_sum = sum(sum(miniLoadCost_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    totalCost_sum = sum(sum(totalCost_sol[:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    totalCostR = totalCost_sum + ICcostTotal        # total cost of residential consumers ICcostTotal is negative
    #benefitsR = transpose(sum(supplyR_sol[:,t,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet, t = TimeSet))*Vi    # expected residential benefits
    benefitsRealized = sum( transpose(sum(supplyRealized[:,t,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))*Vit[:, t] for t = TimeSet )
    #socialWelfareR = benefitsR - totalCostR
    socialWelfareRealized = benefitsRealized - totalCostR
    ProfitsRealized = sum(y_consumer_sol[consumer,option]*Dl[consumer] for option = OptionsSet, consumer = ConsumersSet) *
       sum(SystemProfile[t] for t = TimeSet)- (ICcostTotal +  totalCost_sum)
    revenueCompany = totalCostR + ProfitsRealized
    #netBenefitsR = benefitsR - revenueCompany
    netBenefitsRealized = benefitsRealized - revenueCompany

#=
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
       ICCost = -ICcostTotal
   )
    writetable(string(OutputDir,"/WelfareComparsionPSP_","$(niterationsDD)",".csv"),  WelfareComparsionDF; header=true)
=#
    WelfareRealizedComparsionDF = DataFrame(
       SocialWelfare = socialWelfareRealized,
       ConsumerBenefits = benefitsRealized,
       ConsumerNetBenefits = netBenefitsRealized,
       CompanyRevenue = revenueCompany,
       CompanyProfits = ProfitsRealized,
       CompanyCosts = totalCostR,
       ProductionCost = productionCost_sum,
       MiniLoadCost = miniLoadCost_sum,
       StartUpCost = startupCost_sum,
       ICCost = -ICcostTotal
    )
    writetable(string(OutputDir,"/WelfareRealizedComparsionPSP_","$(niterationsDD)",".csv"),  WelfareRealizedComparsionDF; header=true)

    # output subsProfile
    subsProfileDF = DataFrame(subsProfile)
    writetable(string(OutputDir,"/subsProfile.csv"),  subsProfileDF; header=true)
    ReliabilityRealizedHourlyPSPDF = DataFrame(ReliabilityRealizedHourlyPSP)
    writetable(string(OutputDir,"/ReliabilityRealizedHourlyPSP.csv"),  ReliabilityRealizedHourlyPSPDF; header=false)




    # create DataFrame and output as a csv file
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

    #u_partitionDF =  DataFrame(u_partition_chosen = u_partition_chosen_sol)
    #writetable(string(OutputDir,"/u_partition.csv"),  u_partitionDF; header=true)



    PriceMenuDF2 = DataFrame(Reliability = r_sol,
        #ReliabilityRealized = r_realized,
        #ReliabilityRealized2 = r_realized2,
        Price = p_sol,
        #Subscription = Di
        )
    writetable(string(OutputDir,"/PriceMenuFinal_","$(niterationsDD)",".csv"),  PriceMenuDF2; header=true)
end
    return  r_realized

end
