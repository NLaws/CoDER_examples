


lat, lon = 30.2672, -97.7431  # Austin, TX
η = 0.95
years = 20
discountrate = 0.05
M = 10_000  # default bound for decision variables
CHILLER_COP = 4.55
#***********************#
Sbase=1e3  # best to use Sbase of 1kW because thermal model assumes kW and kJ
#***********************#

#   cost components
pwf = REoptLite.annuity(years, 0.0, discountrate)

clmp = vec(readdlm("./data/cleaned_ercot2019_RTprices.csv", ',', Float64, '\n'));
clmp = abs.(clmp) / 1e3;  # problem is unbounded with negative prices, convert from $/MWh to $/kWh
cpv = 1400
cbkW = 700
cbkWh = 350
ci = repeat([0.25], 8760);

# outdoor temperature and PV production factor
tamb = REoptLite.get_ambient_temperature(lat, lon);
prod_factor = REoptLite.get_pvwatts_prodfactor(lat, lon);

# power flow inputs
LDFinputs = LDF.singlephase38linesInputs(Sbase=Sbase);
loadnodes = collect(keys(LDFinputs.Pload))
# remove some loadnodes to make problem smaller and keep voltage further from lower limit
# loadnodes = loadnodes[3:end-3]  # removes "32", "12", "25",  "18", "30", "3"
# for n in ["32", "12", "25",  "18", "30", "3"]
#     delete!(LDFinputs.Pload, n)
#     delete!(LDFinputs.Qload, n)
# end
LLnodes_withPV = ["9", "22", "31", "34", "17"]
LLnodes_warehouse = ["5", "10", "15"]  # (price responsive refrigeration), can add any of the loadnodes
LLnodes = union(LLnodes_withPV, LLnodes_warehouse)  # all nodes in LL model (that have decisions)

ULnodes_withBESS = ["2", "7", "24"]

profile_names = ["FastFoodRest", "FullServiceRest", "Hospital", "LargeHotel", "LargeOffice", 
    "MediumOffice", "MidriseApartment", "Outpatient", "PrimarySchool", "RetailStore", 
    "SecondarySchool", "SmallHotel", "SmallOffice", "StripMall", "Supermarket", "Warehouse"]
doe_profiles = Dict{String, Vector{Float64}}()
for name in profile_names
    doe_profiles[name] = REoptLite.BuiltInElectricLoad("", name, Real(lat), Real(lon), 2017, nothing, Real[])  
end

rand_names = rand(profile_names, length(loadnodes))

R = 0.00025  # K/kW
C = 1e5   # kJ/K
A = reshape([-1/(R*C)], 1,1)
B = [1/(R*C) 1/C]
u = [tamb zeros(8760)]';  # could replace the zeros vector with endogenous heat input
J = size(B,2)
T_hi = 0
T_lo = -20

function get_LDFinputs(T; v_lolim=0.95)

    LDFinputs.Ntimesteps = T
    load_scaler = 0.7
    for (i, node) in enumerate(loadnodes)
        LDFinputs.Pload[node] = load_scaler * doe_profiles[rand_names[i]][1:T] / Sbase;
        LDFinputs.Qload[node] = load_scaler * doe_profiles[rand_names[i]][1:T] / Sbase * 0.1;
    end
    # pepper some pv into the system
    PVkW = 1e3   # TODO more baseline PV ?
    LDFinputs.Pload["5"] .-= PVkW * prod_factor[1:T] / Sbase;
    LDFinputs.Qload["5"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
    LDFinputs.Pload["31"] .-= PVkW * prod_factor[1:T] / Sbase;
    LDFinputs.Qload["31"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
    LDFinputs.Pload["28"] .-= PVkW * prod_factor[1:T] / Sbase;
    LDFinputs.Qload["28"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
    LDFinputs.Pload["27"] .-= PVkW * prod_factor[1:T] / Sbase;
    LDFinputs.Qload["27"] .-= PVkW * prod_factor[1:T] / Sbase * 0.1;
    
    LDFinputs.v_lolim = v_lolim
    
    peak_load = maximum(sum(values(LDFinputs.Pload)))
    peak_single_load = maximum(maximum(values(LDFinputs.Pload)))
    LDFinputs.P_up_bound =  peak_load * 100
    LDFinputs.Q_up_bound =  peak_load * 10
    LDFinputs.P_lo_bound = -peak_single_load * 100
    LDFinputs.Q_lo_bound = -peak_single_load * 10

    bigM = peak_load * 1e4

    return LDFinputs, bigM
end