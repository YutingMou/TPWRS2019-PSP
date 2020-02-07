# can be used to solve for price and reliability when u_consumer is given
# UC can be decoupled by month or week or day
# different consumers follow different profiles
# The welfare calculated is the welfare when the aggregator knows only the system-level profile
# 2018-10-26
# change the ReadSystemTables function
# 2018-11-04
# Add different scenarios
#2019-03-13
# adapted for the paper revised, to do comparison with ADDD using a few scenaros

HOST = "LENOVO"   #LENOVO, LEMAITRE
SubHorizonChoice = "MONTH"#"DAY"
MIPGapSlave = 1e-3
GasHourlyDerating = 1
WithProfitsConstraint = 1
using JuMP, DataFrames,Gurobi #，CSV

if HOST == "LEMAITRE"
    global test_5D = 0
    NTHREADS = 24
    global const SysDir = "/home/users/y/m/ymou/PSPADMMV20190208"
    SimDir = string(SysDir,"/sourceV2")
    cd(SimDir)
    include("./Constants.jl")
    include("./ReadSystemTables.jl")
    include("./time_cycling.jl")

    PSPMIPSolverChoice = GurobiSolver(
        OutputFlag=1,
        MIPGap=MIPGapSlave,
        #Method = 2,
        Threads=NTHREADS,
        TimeLimit=36000,
        )

elseif HOST == "LENOVO"
    global test_5D = 0
    NTHREADS = 4
    global const SysDir = "C:/Users/Yuting/OneDrive - UCL/Demand Response Business Model/Julia/PSPADMMV20190208"
    SimDir = string(SysDir,"/sourceV2")
    cd(SimDir)
    include("./Constants.jl")
    include("./ReadSystemTables.jl")
    include("./time_cycling.jl")
    PSPMIPSolverChoice = GurobiSolver(
        OutputFlag=1,
        MIPGap=MIPGapSlave,
        #MIPGapAbs = 3e4,
        #Threads=NTHREADS,
        Method = 2,
        #BarConvTol = 1e-6,
        #Crossover = 0,
        TimeLimit=54000,
        )

else
    println("Host is wrongly defined")
end


# load data from the table
Generators, HeatRateCurves, FuelPrice, DlVl,
    SystemProfile, ConsumerProfile,	GasAvailability,
    HistoricalProfiles, ICCost,ICcostHourly,
    WindProfiles, SolarProfiles, u_consumer_sol = ReadSystemTables(InputDirectory)
u_consumer_sol = convert(Array, u_consumer_sol)
DlVl = convert(Array, DlVl)
 Dl = DlVl[:, 1]
 Vl = DlVl[:, 2]
 Vmax = 487
 Vmin = 0
 vb_LB = [Vmin + (Vmax - Vmin)*(n-1)/NPartitions for n = 1:NPartitions]    # lower bound of vb
 vb_UB = [Vmin + (Vmax - Vmin)*n/NPartitions for n = 1:NPartitions] # upper bound of vb
 #vb_LB = vb_UB   # set lower bound and upper bound to be the same
 VB = Array{Float64}(NOptions+1)   # valuation breakpoints parameters
 VB[1] = vb_LB[1]    # the first consumer is the first breakpoint
 for option = 2:NOptions+1
     VB[option] = vb_UB[sum(u_consumer_sol[:,1:option-1])]
 end
 Di = Array{Float64}(NOptions)
 for option = OptionsSet
     Di[option] = sum(u_consumer_sol[consumer,option]*Dl[consumer] for consumer = ConsumersSet)
 end
 Vi = VB[1:NOptions]*0.5+VB[2:NOptions+1]*0.5
 ConsumerProfile = convert(Array, ConsumerProfile)
 SystemProfile = convert(Array, SystemProfile)
# process data
MinRunCapacityDerated = Array{Float64}(NGenerators, Thours)
MaxRunCapacityDerated = Array{Float64}(NGenerators, Thours)
MinRunCapacityDerated .= Generators[:MinRunCapacity]
MaxRunCapacityDerated .= Generators[:MaxRunCapacity]
# gas units are to be derated
GasGeneratorsSet = view(GeneratorsSet, Generators[:FuelGenerator].=="NAT_GAS")
if GasHourlyDerating == 1
    for t = TimeSet
        MinRunCapacityDerated[GasGeneratorsSet,t] =
            MinRunCapacityDerated[GasGeneratorsSet,t]*GasAvailability[:Ratio][t]*GasRatioScaleUp
        MaxRunCapacityDerated[GasGeneratorsSet,t] =
            MaxRunCapacityDerated[GasGeneratorsSet,t]*GasAvailability[:Ratio][t]*GasRatioScaleUp
    end
else
    GasAvailabilityAverage = mean(GasAvailability[:Ratio][TimeSet])
    for t = TimeSet
        MinRunCapacityDerated[GasGeneratorsSet,t] =
            MinRunCapacityDerated[GasGeneratorsSet,t]*GasAvailabilityAverage*GasRatioScaleUp
        MaxRunCapacityDerated[GasGeneratorsSet,t] =
            MaxRunCapacityDerated[GasGeneratorsSet,t]*GasAvailabilityAverage*GasRatioScaleUp
    end

end



ICcostTotal = convert(Array, ICCost)
ICcostTotal = sum(ICcostTotal[:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet)
ICcostTotal = - ICcostTotal[1,1]
############ the constant term in Dual decompostion
Profits = Profits - ICcostTotal # ic cost must be shifted
#############################################################


ICcostHourly = convert(Array, ICcostHourly)
WindProfiles = convert(Array, WindProfiles)
SolarProfiles = convert(Array, SolarProfiles)

# in the 5d version, scale the production level of renewable
if test_5D == 1
    WindEnergyRef = sum(HistoricalProfiles[:Wind][TimeSet])
    SolarEnergyRef = sum(HistoricalProfiles[:Solar][TimeSet])
    WindEnergySce = sum(WindProfiles[TimeSet,:],1)
    SolarEnergySce = sum(SolarProfiles[TimeSet,:],1)
    WindRatio = WindEnergySce/WindEnergyRef
    SolarRatio = SolarEnergySce/SolarEnergyRef
    for scenario = ScenariosSet
        WindProfiles[:,scenario] = WindProfiles[:,scenario]/WindRatio[scenario]
        SolarProfiles[:,scenario] = SolarProfiles[:,scenario]/SolarRatio[scenario]
    end
end


# consumers average valuation is actually different at each hour,
# because the demand prfiles would change
Vit =  zeros(Float64, NOptions, Thours) # average valuation of each option at each hour
for t = TimeSet
    for option = OptionsSet
        Vit[option, t] =
        sum(u_consumer_sol[consumer, option]*Dl[consumer]*Vl[consumer]*ConsumerProfile[consumer, t] for consumer = ConsumersSet)/sum(u_consumer_sol[consumer, option]*Dl[consumer]*ConsumerProfile[consumer, t] for consumer = ConsumersSet)
    end
end


PSPMIP = Model(solver = PSPMIPSolverChoice)
println("Starts at ", now())
flush(STDOUT)

## define decision varioables
@variables PSPMIP begin
    # conventional generators
    u[GeneratorsSet, TimeSet, ScenariosSet], Bin   # on/off binary
    v[GeneratorsSet, TimeSet, ScenariosSet], Bin   # start-up binary
    z[GeneratorsSet, TimeSet, ScenariosSet], Bin   # shut-down binary
    production[GeneratorsSet, TimeSet, ScenariosSet] >= 0
    productionCost[GeneratorsSet, TimeSet, ScenariosSet] >= 0
    startupCost0[GeneratorsSet, TimeSet, ScenariosSet] >=0    # startup cost + startup fuel cost
    miniLoadCost[GeneratorsSet, TimeSet, ScenariosSet] >=0    # regardless of the load, this cost is applicable if the unit is on
    totalCost[TimeSet, ScenariosSet] >= 0 # total cost of all generators at t
    # pumped hydro storage
    0 <= s[TimeSet, ScenariosSet] <= PSEnergy             # energy in pumped storage in MWH
    0 <= p_pump[TimeSet, ScenariosSet] <= PSPumpingMax    # power pumped into pumped storage
    0 <= p_prod[TimeSet, ScenariosSet] <= PSProducingMax  # power production of pumped sorage

    loadShedding[TimeSet, ScenariosSet] >= 0  # in case of IC demand is too high
    renewableShedding[TimeSet, ScenariosSet] >= 0

    p[OptionsSet] >= 0   # price in the menu
    r[OptionsSet] >= 0   # reliability in the menu
    residentialD[TimeSet, ScenariosSet] >= 0   # total supply to residential consuemrs
    supplyR[OptionsSet, TimeSet, ScenariosSet] >= 0   # total supply to residential consuemrs in each option
    z_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for r*u
    y_consumer[ConsumersSet, OptionsSet] >= 0  # auxilary variable for p*u
    surplus[ConsumersSet] >= 0  # gamma in the note
end
println("variables defined at ", now())
flush(STDOUT)
# promised reliability is delivered, the aggregator is only aware of system profile
for option = OptionsSet
    Expression_1 = AffExpr()
    for consumer = ConsumersSet, t = TimeSet
        append!(Expression_1, z_consumer[consumer, option]*Dl[consumer]*SystemProfile[t])
    end
    Expression_2 = AffExpr()
    for scenario = ScenariosSet
        for t = TimeSet
            append!(Expression_2, supplyR[option,t, scenario]*ScenarioProb[scenario])
        end
    end
    @constraints(PSPMIP, begin
        Expression_1 == Expression_2
        #sum(z_consumer[consumer, option]*Dl[consumer]*SystemProfile[t]
        #for consumer = ConsumersSet, t = TimeSet) == sum(supplyR[option,t] for t = TimeSet)
    end)
end
# optimality conditions
@constraints(PSPMIP, begin
    [consumer = ConsumersSet], surplus[consumer] <=
        sum((Vl[consumer]*z_consumer[consumer,option] - y_consumer[consumer,option])
            for option = OptionsSet)
end)
for option = OptionsSet
    for consumer = ConsumersSet
        @constraints(PSPMIP, begin
            surplus[consumer] >= Vl[consumer]*r[option]-p[option]
            # McCormick envelopes
            z_consumer[consumer, option] == r[option]*u_consumer_sol[consumer, option]
            y_consumer[consumer, option] == p[option]*u_consumer_sol[consumer, option]
            #z_consumer[consumer, option] <= u_consumer_sol[consumer, option]
            #z_consumer[consumer, option] <= r[option]
            #z_consumer[consumer, option] >= u_consumer_sol[consumer, option] + r[option] - 1
            #y_consumer[consumer, option] <= u_consumer_sol[consumer, option]*PriceMax
            #y_consumer[consumer, option] <= p[option]
            #y_consumer[consumer, option] >= u_consumer_sol[consumer, option]*PriceMax + p[option] - PriceMax
        end)

    end
end

for consumer = 1:NConsumers-1
    for optionStart = OptionsSet
        @constraints(PSPMIP, begin
            # Valid cuts on z_consumer and y_consumer
            sum(z_consumer[consumer, option] for option = optionStart:NOptions)  <=
    		    sum(z_consumer[consumer+1, option] for option = optionStart:NOptions)
            sum(y_consumer[consumer, option] for option = optionStart:NOptions)  <=
          		sum(y_consumer[consumer+1, option] for option = optionStart:NOptions)
        end)
    end
    @constraints(PSPMIP, begin
        # surplus of each consumer should also be increasing
        surplus[consumer] <= surplus[consumer+1]
    end)
end

@constraints(PSPMIP, begin
    # price and reliability in the menu should be (nondecreasing) increasing
    [option = 1:NOptions-1],
       r[option] <= r[option+1]
    [option = 1:NOptions-1],
        p[option] <= p[option+1]
    # PS empty everyday at mightnight
    [d = 1:Tdays, scenario = ScenariosSet], s[d*24, scenario] <= 1
end)

println("consumers constraints are defined at ", now())
flush(STDOUT)


if SubHorizonChoice == "MONTH"
    for scenario = ScenariosSet, month = 1:NMonth
        for t = StartHour[month]:EndHour[month]
            for g = GeneratorsSet
                @constraints(PSPMIP, begin
                    u[g, time_cycling(t-1, StartHour[month], EndHour[month]), scenario]+v[g, t, scenario] -
                        z[g, t, scenario] == u[g, t, scenario]
                    sum(v[g, time_cycling(y, StartHour[month], EndHour[month]), scenario]
                        for y = t - Generators[:UT][g] + 1:t ) <= u[g, t, scenario]
                    sum(v[g,time_cycling(yy, StartHour[month], EndHour[month]), scenario]
                        for yy = t + 1 : t + Generators[:DT][g]) <= 1 - u[g, t, scenario]
                end)
            end
        end
    end
elseif SubHorizonChoice == "WEEK"
    for scenario = ScenariosSet, week = 1:NWeek
        for t = StartHourX[week]:EndHourX[week]
            for g = GeneratorsSet
                @constraints(PSPMIP, begin
                    u[g, time_cycling(t-1, StartHourX[week], EndHourX[week]), scenario]+v[g, t, scenario] - z[g, t, scenario] == u[g, t, scenario]
                    sum(v[g, time_cycling(y, StartHourX[week], EndHourX[week]), scenario]
                    for y = t - Generators[:UT][g] + 1:t ) <= u[g, t, scenario]
                    sum(v[g,time_cycling(yy, StartHourX[week], EndHourX[week]), scenario]
                    for yy = t + 1 : t + Generators[:DT][g]) <= 1 - u[g, t, scenario]
                end)
            end
        end
    end
elseif SubHorizonChoice == "DAY"
    for scenario = ScenariosSet,  d = 1:Tdays
        for t = StartHourX[d]:EndHourX[d]
            for g = GeneratorsSet
                @constraints(PSPMIP, begin
                    u[g, time_cycling(t-1, StartHourX[d], EndHourX[d]), scenario]+v[g, t, scenario] - z[g, t, scenario] == u[g, t, scenario]
                    sum(v[g, time_cycling(y, StartHourX[d], EndHourX[d]), scenario]
                    for y = t - Generators[:UT][g] + 1:t ) <= u[g, t, scenario]
                    sum(v[g,time_cycling(yy, StartHourX[d], EndHourX[d]), scenario]
                    for yy = t + 1 : t + Generators[:DT][g]) <= 1 - u[g, t, scenario]
                end)
            end
        end
    end
elseif SubHorizonChoice == "YEAR"
    for scenario = ScenariosSet, t = TimeSet
        for g = GeneratorsSet
            @constraints(PSPMIP, begin
                u[g, time_cycling(t-1, 1, Thours), scenario]+v[g, t, scenario] - z[g, t, scenario] == u[g, t, scenario]
                sum(v[g, time_cycling(y, 1, Thours), scenario]
                for y = t - Generators[:UT][g] + 1:t ) <= u[g, t, scenario]
                sum(v[g,time_cycling(yy, 1, Thours), scenario]
                for yy = t + 1 : t + Generators[:DT][g]) <= 1 - u[g, t, scenario]
            end)
        end
    end
end


## PS
for scenario = ScenariosSet, month = 1:NMonth
    @constraints(PSPMIP, begin
    # initial of pumped hydro
        s[StartHour[month], scenario] ==
        0 + p_pump[StartHour[month], scenario]*PSEff - p_prod[StartHour[month], scenario]
    end)

    for t = StartHour[month]+1:EndHour[month]
        @constraints(PSPMIP, begin
            s[t,scenario] == s[t-1,scenario] + p_pump[t,scenario]*PSEff - p_prod[t,scenario]
        end)
    end
end
for scenario = ScenariosSet, t = TimeSet
    # conventional gnerators' production bounds
    for g = GeneratorsSet
        @constraints(PSPMIP, begin
            production[g,t,scenario] >= MinRunCapacityDerated[g,t]*u[g,t,scenario]
            production[g,t,scenario] <= MaxRunCapacityDerated[g,t]*u[g,t,scenario]
        end)
    end
    # supply cannot exceeds consumer Subscription
    for option = OptionsSet
        @constraints(PSPMIP, begin
        supplyR[option,t,scenario] <= sum(Dl[consumer]*u_consumer_sol[consumer,option]*SystemProfile[t]
            for consumer = ConsumersSet)
        end)
    end

    @constraints(PSPMIP, begin
        # total supply to all options = residential demand
        residentialD[t,scenario] == sum(supplyR[option,t,scenario] for option = OptionsSet)
        # supply demand balance
        sum(production[g,t,scenario] for g = GeneratorsSet)  + loadShedding[t,scenario] +
        WindProfiles[t,scenario] + SolarProfiles[t,scenario] - renewableShedding[t,scenario] +
            HistoricalProfiles[:Import][t] + p_prod[t,scenario] - p_pump[t,scenario] + 77 ==
				residentialD[t,scenario] + HistoricalProfiles[:ICDemand][t]
    end)
end

## model the costs of conventional generators
for g = GeneratorsSet
    ## collect heat rate curve and fuel price data for generator g
	hrc = view(HeatRateCurves, HeatRateCurves[:SeriesHRC] .==
		Generators[:HRSGenerator][g])
	fprows = .&(FuelPrice[:Fuel] .== Generators[:FuelGenerator][g],
		FuelPrice[:FuelHub] .== Generators[:FHGenerator][g])
    ## cost curve is piecewise Linear
    FuelPrice_g = view(FuelPrice[:Price], fprows)[1]
    for scenario = ScenariosSet, t = TimeSet
        @constraints(PSPMIP, begin
            # productionCost
            productionCost[g,t,scenario] .>=
            FuelPrice_g*(hrc[:SlopeHRC]*production[g,t,scenario] + hrc[:InterceptHRC]*u[g,t,scenario])
            # miniLoadCost
            miniLoadCost[g,t,scenario] ==
            Generators[:NoLoadConsumption][g]*MinRunCapacityDerated[g,t]*FuelPrice_g*u[g,t,scenario]
            # startup costs
            startupCost0[g,t,scenario] ==
            (Generators[:StartupCost][g] + Generators[:StartupFuel][g]*FuelPrice_g)*MinRunCapacityDerated[g,t]*v[g,t,scenario]
        end)
    end
end
println("Producer constraints are defined at ", now())
flush(STDOUT)
## profits requirement
for scenario = ScenariosSet
    @constraints(PSPMIP, begin
        [t = TimeSet],
        totalCost[t,scenario] == sum( (productionCost[g,t,scenario] + startupCost0[g,t,scenario] + miniLoadCost[g,t,scenario])
            for g = GeneratorsSet)
        #sum(y_consumer[consumer,option]*Dl[consumer] for option = OptionsSet, consumer = ConsumersSet)  ==
        #    (Profits + ICcostTotal +  sum(totalCost[t] for t = TimeSet))/sum(SystemProfile[t] for t = TimeSet)
    end)
end
## profits requirement - compile faster
if WithProfitsConstraint == 1
    ProfitsExpression_1 = AffExpr()
    for option = OptionsSet, consumer = ConsumersSet
        append!(ProfitsExpression_1, y_consumer[consumer,option]*Dl[consumer])
    end
    ProfitsExpression_2 = AffExpr()
    for scenario = ScenariosSet
        for t = TimeSet
            append!(ProfitsExpression_2, totalCost[t,scenario]*ScenarioProb[scenario])
        end
    end
    @constraints(PSPMIP, begin
        ProfitsExpression_1 == (Profits +  ProfitsExpression_2)/sum(SystemProfile)
    end)
end


# the sum takes too much time to compile, so I use append
objExpression_1 = AffExpr()
for scenario = ScenariosSet, t = TimeSet
    append!(objExpression_1, ScenarioProb[scenario]*(totalCost[t,scenario] + loadShedding[t,scenario]*VOLL))
end
objExpression_2 = AffExpr()
for scenario = ScenariosSet, t = TimeSet, option = OptionsSet
    append!(objExpression_2, ScenarioProb[scenario]*supplyR[option,t,scenario]*Vi[option])
end
@objective(PSPMIP, Max, - objExpression_1 + objExpression_2  )
#@objective(PSPMIP, Max, - sum((totalCost[t] + loadShedding[t]*VOLL)  for t=TimeSet)
#    + sum(supplyR[option,t]*Vi[option] for t = TimeSet, option = OptionsSet) + ICcostTotal  )

println("Solvs at ", now())
flush(STDOUT)

status_PSPMIP = solve(PSPMIP, relaxation = false)
println("Solved at ", now())
flush(STDOUT)
#if status_PSPMIP == :Optimal
    #currentslaveobjval = getobjectivevalue(PSPMIP)
    #println("Objective value of PSPMIP: ", currentslaveobjval)
    println("Price: ", getvalue(p[:]))
    println("Reliablity: ", getvalue(r[:]))

    #### output everything from it
    # create directory
    if test_5D == 1
        OutputDir = string(SysDir,"/output/5D",Dates.format(now(), "yyyy-mm-dd-HH-MM-SS"))
    else
        OutputDir = string(SysDir,"/output/FullYear",Dates.format(now(), "yyyy-mm-dd-HH-MM-SS"))
    end
    mkpath(OutputDir)

    r_sol = getvalue(r[:])
    p_sol = getvalue(p[:])
    u_sol = Int.(round.(getvalue(u[:,:,:])))
    v_sol = Int.(round.(getvalue(v[:,:,:])))
    z_sol = Int.(round.(getvalue(z[:,:,:])))
    production_sol = getvalue(production[:,:,:]);
    productionCost_sol = getvalue(productionCost[:,:,:])
    startupCost0_sol = getvalue(startupCost0[:,:,:])
    miniLoadCost_sol = getvalue(miniLoadCost[:,:,:])
    s_sol = getvalue(s[:,:])
    p_pump_sol = getvalue(p_pump[:,:])
    p_prod_sol = getvalue(p_prod[:,:])
    loadShedding_sol = getvalue(loadShedding[:,:]);
    renewableShedding_sol = getvalue(renewableShedding[:,:]);
    totalCost_sol = getvalue(totalCost[:,:])
    residentialD_sol = getvalue(residentialD[:,:])
    supplyR_sol = getvalue(supplyR[:,:,:]);
    y_consumer_sol = getvalue(y_consumer[:,:])

    # subscription profile of each option
    subsProfile = Array{Float64}(NOptions, Thours)
    for option = OptionsSet
        subsProfile[option, :] = sum(
            u_consumer_sol[consumer, option]*Dl[consumer]*ConsumerProfile[consumer,:]
        for consumer = ConsumersSet)
    end

    # calculate realized reliability, need to consider real subscription profile and supply
    # real supply to each option, related to real subscription profile
    supplyRealized = Array{Float64}(NOptions, Thours, NScenarios)
    for scenario = ScenariosSet
        residentialD_temp = residentialD_sol[:,scenario]
        for option = NOptions:-1:1
            supplyRealized[option,:,scenario] = min.(subsProfile[option, :], residentialD_temp)
            residentialD_temp = residentialD_temp - supplyRealized[option,:,scenario]
        end
    end


    r_realized =sum(sum(supplyRealized[:,:,scenario],2)*ScenarioProb[scenario] for scenario = ScenariosSet)./sum(subsProfile,2)
    r_realized = r_realized[:]
    ReliabilityRealizedHourlyPSP = sum(supplyRealized[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet)./subsProfile

    println("ReliabilityRealized: ", r_realized)


    productionCost_sum = sum(sum(productionCost_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    startupCost_sum = sum(sum(startupCost0_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    miniLoadCost_sum = sum(sum(miniLoadCost_sol[:,:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    totalCost_sum = sum(sum(totalCost_sol[:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    totalCostR = totalCost_sum - ICcostTotal        # total cost of residential consumers ICcostTotal is negative
    benefitsR = transpose(sum(supplyR_sol[:,t,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet, t = TimeSet))*Vi    # expected residential benefits
    benefitsRealized = sum( transpose(sum(supplyRealized[:,t,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))*Vit[:, t] for t = TimeSet )
    socialWelfareR = benefitsR - totalCostR
    socialWelfareRealized = benefitsRealized - totalCostR
    ProfitsRealized = sum(y_consumer_sol[consumer,option]*Dl[consumer] for option = OptionsSet, consumer = ConsumersSet) *
        sum(SystemProfile[t] for t = TimeSet)- (totalCost_sum - ICcostTotal )
    revenueCompany = totalCostR + ProfitsRealized
    netBenefitsR = benefitsR - revenueCompany
    netBenefitsRealized = benefitsRealized - revenueCompany

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
    writetable(string(OutputDir,"/WelfareComparsionPSP.csv"),  WelfareComparsionDF; header=true)

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
        ICCost = ICcostTotal
    )
    writetable(string(OutputDir,"/WelfareRealizedComparsionPSP.csv"),  WelfareRealizedComparsionDF; header=true)

    # output subsProfile
    subsProfileDF = DataFrame(subsProfile)
    writetable(string(OutputDir,"/subsProfile.csv"),  subsProfileDF; header=true)
    ReliabilityRealizedHourlyPSPDF = DataFrame(ReliabilityRealizedHourlyPSP)
    writetable(string(OutputDir,"/ReliabilityRealizedHourlyPSP.csv"),  ReliabilityRealizedHourlyPSPDF; header=false)
    # output price menu
    PriceMenuDF = DataFrame(Reliability = r_sol,
                ReliabilityRealized = r_realized,
                Price = p_sol,
                Subscription = Di)
    writetable(string(OutputDir,"/PriceMenu.csv"),  PriceMenuDF; header=true)

    totalCostDF = DataFrame(totalCost_sol)
    writetable(string(OutputDir,"/totalCost.csv"),  totalCostDF; header=true)

    supplyR_sol = reshape(supplyR_sol, NOptions, Thours*NScenarios)
    supplyRDF = DataFrame(supplyR_sol)
    writetable(string(OutputDir,"/supplyR.csv"),  supplyRDF; header=true)

    supplyRealized = reshape(supplyRealized, NOptions, Thours*NScenarios)
    supplyRealizedDF = DataFrame(supplyRealized)
    writetable(string(OutputDir,"/supplyRealized.csv"),  supplyRealizedDF; header=true)

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

    residentialDDF = DataFrame(residentialD_sol)
    writetable(string(OutputDir,"/residentialD.csv"),  residentialDDF; header=true)


#=
    # calculate hourly socialWelfare
    WelfareRealizedHourly = zeros(Float64, Thours)
    for t = TimeSet
        WelfareRealizedHourly[t] = transpose(supplyRealized[:,t])*Vit[:, t] - totalCost_sol[t] + ICcostHourly[t]
    end

    WelfareRealizedHourlyDF = DataFrame(WelfareRealizedHourly)
    writetable(string(OutputDir,"/WelfareRealizedHourlyPSP.csv"),  WelfareRealizedHourlyDF; header=true)
=#
    #println("supplyR_sol: ", transpose(sum(supplyR_sol,2)))
    #println("supplyRealized: ", transpose(sum(supplyRealized,2)))
    #println("supplyR_sol_sum: ", sum(transpose(sum(supplyR_sol,2))))
    #println("supplyRealized_sum: ", sum(transpose(sum(supplyRealized,2))))

    println("socialWelfareR: ", socialWelfareR)
    println("socialWelfareRealized: ", socialWelfareRealized)
    println("Ends at ", now())


    flush(STDOUT)


#end
