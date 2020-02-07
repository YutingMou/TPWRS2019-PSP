global const InputDirectory = string(SysDir,"/input")

if test_5D == 1
    global  ALPHA = 1e-7  # step size
    global  GasRatioScaleUp = 1.08
    global const NMonth = 1
    global const Tdays = 5#1#5
    global const Thours = Tdays*24
    global const TimeSet = 1:Thours
    global const Days = [5]#[1]#[5]
    global const StartHour = [1]
    global const EndHour = [120]#[24]#[120]
    global  Profits = 0#20.8184e6 #22.613e6#
    global const DaysX = ones(Int64,5)
    global const StartHourX = [1+24*i for i = 0:4]
    global const EndHourX = [24*i for i = 1:5]
    global  ySCALE = 1 #1e4/2#  # scale up y_consumer
    global  zSCALE = 1 #1e5/2#1e4   # scale up z_consumer
    global  sSCALE = [1, 1, 1, 1, 1]#1e1*[12, 13, 12, 8, 12]#0.5*1e2*[12, 6, 8, 5, 20]  #1e2/2#1e1   # scale up supplyR
    #global SCALE_dual_vb = 1e-4
    #global SCALE_dual_profits = 1e1
    #global SCALE_dual_reliability = 1
    global const PriceMax = 200                     # upper bound for price
    global const NSubhoriozns = 5

else
    global const ySCALE = 1   # scale up y_consumer
    global const zSCALE = 1   # scale up z_consumer
    global const sSCALE = [1, 1, 1, 1, 1]#1e2*[5, 5, 5, 9, 5]    # scale up supplyR
    global const PriceMax = 200                     # upper bound for price

    global const GasRatioScaleUp = 1.25
    global  Profits = 466.0552e6
    global const NMonth = 12
    global const NWeek = 52
    global const Tdays = 365
    global const Thours = Tdays*24
    global const TimeSet = 1:Thours
    global const Days = [31,28,31,30,31,30,31,31,30,31,30,31]
    global const StartHour = [1,745,1417,2161,2881,3625,4345,5089,5833,6553,7297,8017]
    global const EndHour = [744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]
    if SubHorizonChoice == "MONTH"
        global const ALPHA = 1e-8
        global const DaysX = [31,28,31,30,31,30,31,31,30,31,30,31]
        global const StartHourX = [1,745,1417,2161,2881,3625,4345,5089,5833,6553,7297,8017]
        global const EndHourX = [744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]
        global const NSubhoriozns = 12
    elseif SubHorizonChoice == "WEEK"
        global const ALPHA = 1e-6
        global const DaysX = 7*ones(Int64,51)   # how many days are there in the subproblem
                    push!(DaysX, 8)  # the last 'week' has 8 days
        global const StartHourX = [1+168*i for i = 0:51]
        global const EndHourX = [168*i for i = 1:51]
                    push!(EndHourX, 8760)
		global const NSubhoriozns = 52
    else SubHorizonChoice == "DAY"
        global const ALPHA = 1e-6
        global const DaysX = ones(Int64,365)
        global const StartHourX = [1+24*i for i = 0:364]
        global const EndHourX = [24*i for i = 1:365]
    end



end

global const VOLL = 3000.0			# value of lost load, for load shedding
#global const ZEROTOL = 1E-8

# specifications of pumped hydro in Belgium
#(Paper: An overview of large-scale stationary electricity storage plants in Europe)
#PSProducingRamp = 1166		# hourly
#PSPumpingRamp = 1103		# hourly
global const PSPumpingMax = 1196
global const PSProducingMax = 1301
global const PSEnergy = 5710             # MWh
global const PSEff = 0.765
global const NGenerators = 55            # number of generators
global const NOptions = 5                # number of options in the price menu
global  NConsumers = 1000
global  NScenarios = 1            # number of scenarios for renewable production

global NPartitions = 1000
# it is not necessary to use 1000 partitions for the first breakpoint
global NPartitionsChosen = 1
global PartitionsSet = 1:NPartitionsChosen#NPartitions

global supplyR_U = 10000    # upper bound on the supply to each option
global const GeneratorsSet = 1:NGenerators
global const OptionsSet = 1:NOptions
global const ConsumersSet = 1:NConsumers

global  ScenariosSet = 1:NScenarios
global  ScenarioProb = 1/NScenarios*ones(1, NScenarios) # probality



println("CONSTANTS are defined")
