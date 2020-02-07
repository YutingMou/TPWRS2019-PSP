tic()   # record time to load moduels and data
HOST = "LENOVO"   #LENOVO, LEMAITRE

if HOST == "LEMAITRE"
    println("Date/time: ", now(), "\n\n")
    println("N.O. of proceses: ", nprocs())
    flush(STDOUT)
    @everywhere begin
        global  test_5D = 0
        global  SubHorizonChoice = "MONTH"  # "DAY", "WEEK", "MONTH"
        global  MIPGapSlave = 0.005
        global  MIPGapMaster = 0.001
        global  MaxIterationDD = 50   # max iterations of dual Decomposition in the X update
        global  GapTarget = 1e-1   # % in percentage
        println(gethostname(), " is started")
        using JuMP, DataFrames, Gurobi
        NTHREADS = 1
        global const SysDir = "/home/users/y/m/ymou/PSPDDCPV20190206"
        SimDir = string(SysDir,"/source")
        cd(SimDir)
        include("./Constants.jl")
        include("./ReadSystemTables.jl")
        include("./MISC.jl")
        include("./DDOutput.jl")
        SlaveXSolverChoice = GurobiSolver(
            OutputFlag=0,
            Threads=NTHREADS
            )
        SlaveXMIPSolverChoice = GurobiSolver(
            OutputFlag=0, MIPGap=MIPGapSlave, Threads=NTHREADS,
                TimeLimit=36000)
    end

elseif HOST == "LENOVO"
    println("Date/time: ", now(), "\n\n")
    addprocs(6- nprocs())
    println("N.O. of proceses: ", nprocs())
    flush(STDOUT)
    @everywhere begin
        global  test_5D = 1
        global  SubHorizonChoice = "DAY"  # "DAY", "WEEK", "MONTH"
        global  MIPGapSlave = 0.001
        global  MIPGapMaster = 0.001
        global  MaxIterationDD = 50   # max iterations of DD in the X update
        global  GapTarget = 1e-1  # % in percentage
        println(gethostname(), " is started")
        using JuMP, DataFrames, Gurobi
        NTHREADS = 1
        global const SysDir = "C:/Users/Yuting/OneDrive - UCL/Demand Response Business Model/Julia/PSPDDCPV20190206"
        SimDir = string(SysDir,"/source")
        cd(SimDir)
        include("./Constants.jl")
        include("./ReadSystemTables.jl")
        include("./MISC.jl")
        include("./DDOutput.jl")
        SlaveXSolverChoice = GurobiSolver(
            OutputFlag=0,
            Threads=NTHREADS)
        SlaveXMIPSolverChoice = GurobiSolver(
            OutputFlag=0, MIPGap=MIPGapSlave, Threads=NTHREADS,
                TimeLimit=36000)
    end

else
    println("Host is wrongly defined")
end

@eval @everywhere begin
    GasHourlyDerating = 1
    Generators, HeatRateCurves, FuelPrice, DlVl,
        SystemProfile, ConsumerProfile,	GasAvailability,
        HistoricalProfiles, ICCost,ICcostHourly,
        WindProfiles, SolarProfiles, u_consumer_initial = ReadSystemTables(InputDirectory)
    u_consumer_initial = convert(Array, u_consumer_initial)
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
        VB[option] = vb_UB[sum(u_consumer_initial[:,1:option-1])]
    end

    Di = Array{Float64}(NOptions)
    for option = OptionsSet
        Di[option] = sum(u_consumer_initial[consumer,option]*Dl[consumer] for consumer = ConsumersSet)
    end

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

##################################
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
    NSubproblemsX = length(DaysX)   # number of subproblems in the X updates in each scenario
    NSubproblemsX = NSubproblemsX*NScenarios    # total number of subproblems
end

# Slave X is solved in each process,
@everywhere begin
    include("./SlaveXScaled.jl")
end
include("./isDDfinished.jl")
include("./SlaveZScaled.jl")
include("./Master.jl")


println("Loading moduels and data takes: ", round(toq(),2))
flush(STDOUT)
tic()

############################################################################
println("*** Solving  Using Dual Decomposition ***")
flush(STDOUT)
println("*** dual initialzied as 0 ***")
# three sets of dual variables
global dual_vb_k = zeros(Float64, NSubproblemsX+1)
global dual_vb_vec = zeros(Float64, NSubproblemsX+1, MaxIterationDD)
global dual_reliability_k = zeros(Float64, NOptions)
global dual_reliability_vec = zeros(Float64, NOptions, MaxIterationDD)
global dual_profits_k = 0
global dual_profits_vec = zeros(Float64, MaxIterationDD)
# thres sets of residuals, used as subgradient
global residual_vb = zeros(Float64, NSubproblemsX+1)
global residual_vb_vec = zeros(Float64, NSubproblemsX+1, MaxIterationDD)
global residual_reliability = zeros(Float64, NOptions)
global residual_reliability_vec = zeros(Float64, NOptions, MaxIterationDD)
global residual_profits = 0
global residual_profits_vec = zeros(Float64, MaxIterationDD)


# gap between upper bound and lower bound
#global gap_current = 100
#global gap_vec = 100*ones(Float64, MaxIterationDD)
# record upper bound of the problem
#global UB_vec =  zeros(Float64, MaxIterationDD)
# record the lower bound of problem, also the dual function value
LB_vec =  zeros(Float64, MaxIterationDD)
objvalue_master_vec =  zeros(Float64, MaxIterationDD)

#global ProfitsRealized_vec = 100*ones(Float64, MaxIterationDD)
#global r_realized_last = zeros(Float64, NOptions, 1)
#global r_realized_vec = zeros(Float64, NOptions, MaxIterationDD)

#global r_sol_vec = zeros(Float64, NOptions, MaxIterationDD)
#global p_sol_vec = zeros(Float64, NOptions, MaxIterationDD)

slaveZ, p, r,u_consumer, z_consumer, y_consumer, surplus = DefineSlaveZ()


niterationsDD = 1
println("DualDecomposition starts at ", now())
DDDateTime = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
flush(STDOUT)

# while DD is not finished, keep running it
while ~isDDfinished(niterationsDD, 1)

    println("***** DualDecomposition iterations: ", niterationsDD)
    flush(STDOUT)
    #println("dual_vb_k: ", dual_vb_k)
    if niterationsDD >1
         status_master = solve(master)

         objvalue_master = getobjectivevalue(master)
         println("master is  solved,  Objective: ", objvalue_master)
         objvalue_master_vec[niterationsDD] = objvalue_master
         dual_profits_k = getvalue(dual_profits_master)
         dual_reliability_k = getvalue(dual_reliability_master[:])
         #dual_vb_k =  getvalue(dual_vb_master)
     else
         println("start DDCP with initilized dual variables")
     end
     println("dual_profits_k: ", dual_profits_k)
     println("dual_reliability_k: ", dual_reliability_k)
     flush(STDOUT)

     # store the evolution of dual variables
     #dual_vb_vec[:, niterationsDD] = dual_vb_k
     dual_reliability_vec[:, niterationsDD] = dual_reliability_k
     dual_profits_vec[niterationsDD] = dual_profits_k

################################################################################
################################################################################
    println("***** Preparing input for DD X (each monthly subproblem)")
    flush(STDOUT)

    inputSlaveXVec = []
    for subproblemIdx = 1:NSubproblemsX
        # calculate the month index and scenario index based on
        MonthIdx, ScenarioIdx= returnIDX(subproblemIdx,NSubhoriozns)

        push!(inputSlaveXVec, [subproblemIdx, MonthIdx,ScenarioIdx,
            niterationsDD, dual_vb_k, dual_reliability_k, dual_profits_k])
    end

    println("***** running DualDecomposition X Update in Parallel")
    flush(STDOUT)

    #global SlaveX_solutions
    SlaveX_solutions = pmap(SolveSlaveX, inputSlaveXVec)
    #println("***** All subproblems of X solved, Updating X1")
    # define the solution of X, collecting them together
    totalCost_X_sol =  zeros(Float64, Thours, NScenarios)#Array{Float64}(Thours)
    supplyR_X_sol = zeros(Float64, NOptions, Thours, NScenarios) #Array{Float64}(NOptions, Thours)
    vb_X_sol = zeros(Float64, NSubproblemsX)
    objval_slaveX_sol = 0.0
    # solution from the X updates, (could be in shared array)
    # after convergence of DD, we take the average
    # recover the solution form X SubProblems
    for subproblemIdx = 1:NSubproblemsX
        MonthIdx, ScenarioIdx= returnIDX(subproblemIdx,NSubhoriozns)

        totalCost_X_sol[StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][2]
        supplyR_X_sol[:,StartHourX[MonthIdx]:EndHourX[MonthIdx],ScenarioIdx] =
            SlaveX_solutions[subproblemIdx][3]
        #vb_X_sol[subproblemIdx] =  SlaveX_solutions[subproblemIdx][16]
        objval_slaveX_sol = objval_slaveX_sol + ScenarioProb[ScenarioIdx]*SlaveX_solutions[subproblemIdx][16]

    end

#println("vb_X_sol: ", vb_X_sol)
################################################################################
################################################################################

    println("***** Solving the consumer model")
    flush(STDOUT)

    # call SolveSlaveZ
    r_sol, p_sol,u_consumer_sol, y_consumer_sol, objval_slaveZ_sol =
        SolveSlaveZ(u_consumer_initial, slaveZ, p, r, u_consumer, z_consumer, y_consumer, surplus,
        niterationsDD, dual_reliability_k, dual_profits_k)

    LB_current = objval_slaveX_sol + objval_slaveZ_sol - dual_profits_k*Profits
    LB_vec[niterationsDD] = LB_current

    println("***** LB_vec: ", LB_vec)
    println("***** objvalue_master_vec: ", objvalue_master_vec)
    flush(STDOUT)



################################################################################
################################################################################
    println("***** calculate residuals and generate new cuts")
    flush(STDOUT)

    ## update values of Z
    ProfitsRealized = sum(p_sol[option]*Di[option] for option = OptionsSet)*sum(SystemProfile)-
         sum(sum(totalCost_X_sol[:,scenario]*ScenarioProb[scenario] for scenario = ScenariosSet))
    # residual of profits
    residual_profits =   ProfitsRealized - Profits
    residual_profits_vec[niterationsDD] = residual_profits

    # energySubscribed in each option, to be used to get residual
    energySubscribed = zeros(Float64, NOptions)
    for option = OptionsSet
        energySubscribed[option] = sum(SystemProfile)*Di[option]*r_sol[option]
    end
    # residual of reliability constraints
    energySupplied = zeros(Float64, NOptions)
    for option = OptionsSet
        energySupplied[option] =
            sum(supplyR_X_sol[option,t,scenario]*ScenarioProb[scenario] for t = TimeSet, scenario = ScenariosSet)
    end
    residual_reliability = energySubscribed - energySupplied
    #println("energySubscribed: ", energySubscribed)
    #println("energySupplied: ", energySupplied)
    #println("ProfitsRealized: ", ProfitsRealized)

    #println("residual_vb: ", residual_vb)
    println("residual_reliability: ", residual_reliability)
    println("residual_profits: ", residual_profits)
    flush(STDOUT)

    # add a new cut to the master problem
    @constraint(master, z_master/1e6 <= LB_vec[niterationsDD]/1e6 +
        residual_profits/1e6*(dual_profits_master - dual_profits_k) +
        sum(residual_reliability[option]*(dual_reliability_master[option] - dual_reliability_k[option])/1e6
            for option = OptionsSet)
        )

    # output the results from this DD iterations
    #=
    DDOutput(NSubproblemsX, niterationsDD, SlaveX_solutions,
        r_sol, r_sol_vec, p_sol, p_sol_vec, r_realized_vec, ProfitsRealized_vec,residual_profits_vec,
        dual_vb_vec, dual_reliability_vec, dual_profits_vec,
        u_consumer_sol, y_consumer_sol, UB_vec, DDDateTime)

        =#
    # update iteration number of DD
    println("DD iteration until ",niterationsDD, " runs to ", now())
    flush(STDOUT)

    if test_5D == 1
        OutputDir = string(SysDir,"/output/DD5D", DDDateTime ,"/$(niterationsDD)")
    else
        OutputDir = string(SysDir,"/output/DDFullYear", DDDateTime ,"/$(niterationsDD)")
    end
    ### output the bounds
    if true
        mkpath(OutputDir)
    dual_reliability_vecDF = DataFrame(dual_reliability_vec)
    writetable(string(OutputDir,"/dual_reliability_vec.csv"), dual_reliability_vecDF; header=true)

    dual_profits_vecDF = DataFrame(dual_profits_vec)
    writetable(string(OutputDir,"/dual_profits_vec.csv"), dual_profits_vecDF; header=true)

    BoundsDF = DataFrame(
                        LB_vec = LB_vec,
                        objvalue_master_vec = objvalue_master_vec
                        )
    writetable(string(OutputDir,"/Bounds.csv"),  BoundsDF; header=true)
    #LB_vecDF = DataFrame(LB_vec)
    #writetable(string(OutputDir,"/LB_vec.csv"), LB_vecDF; header=true)

    #objvalue_master_vecDF = DataFrame(objvalue_master_vec)
    #writetable(string(OutputDir,"/objvalue_master_vec.csv"), objvalue_master_vecDF; header=true)

    end


    niterationsDD = niterationsDD + 1
end


println("It", " takes ", round(toq()/3600,2), " hours in total")
flush(STDOUT)


println("*** ALL Done ! ***")
println("Date/time: ", now(), "\n\n")
flush(STDOUT)

### output the bounds
if false
DDDateTime = Dates.format(now(), "yyyy-mm-dd-HH-MM-SS")
dual_reliability_vecDF = DataFrame(dual_reliability_vec)
writetable(string(SysDir,"/dual_reliability_vec", DDDateTime,".csv"), dual_reliability_vecDF; header=true)

dual_profits_vecDF = DataFrame(dual_profits_vec)
writetable(string(SysDir,"/dual_profits_vec", DDDateTime,".csv"), dual_profits_vecDF; header=true)

LB_vecDF = DataFrame(LB_vec)
writetable(string(SysDir,"/LB_vec", DDDateTime,".csv"), LB_vecDF; header=true)

objvalue_master_vecDF = DataFrame(objvalue_master_vec)
writetable(string(SysDir,"/objvalue_master_vec", DDDateTime,".csv"), objvalue_master_vecDF; header=true)

end
