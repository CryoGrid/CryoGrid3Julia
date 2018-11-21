#**************************************************************************************
#
#  $$$$$$\                                 $$$$$$\  $$$$$$$\  $$$$$$\ $$$$$$$\   $$$$$$\
# $$  __$$\                               $$  __$$\ $$  __$$\ \_$$  _|$$  __$$\ $$ ___$$\
# $$ /  \__| $$$$$$\  $$\   $$\  $$$$$$\  $$ /  \__|$$ |  $$ |  $$ |  $$ |  $$ |\_/   $$ |
# $$ |      $$  __$$\ $$ |  $$ |$$  __$$\ $$ |$$$$\ $$$$$$$  |  $$ |  $$ |  $$ |  $$$$$ /
# $$ |      $$ |  \__|$$ |  $$ |$$ /  $$ |$$ |\_$$ |$$  __$$<   $$ |  $$ |  $$ |  \___$$\
# $$ |  $$\ $$ |      $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |  $$ |  $$ |  $$ |  $$ |$$\   $$ |
# \$$$$$$  |$$ |      \$$$$$$$ |\$$$$$$  |\$$$$$$  |$$ |  $$ |$$$$$$\ $$$$$$$  |\$$$$$$  |
#  \______/ \__|       \____$$ | \______/  \______/ \__|  \__|\______|\_______/  \______/
#                     $$\   $$ |
#                     \$$$$$$  |
#                      \______/
#
#**************************************************************************************
# -------------------------------------------------------------------------
# CryoGRID3
# main script for running the model
#
# Developed by: S. Westermann and M. Langer 2015
#
# -------------------------------------------------------------------------
#module CG3
include("matlab.jl")
include("cryoGridTechnical.jl")
include("M1.jl")
include("cryoGridSEB.jl")
include("cryoGridSoil.jl")
include("cryoGridSnow.jl")
include("cryoGridInfiltrationUnfrozenSoil.jl")
include("cryoGridInitialize.jl")

using MAT
using cryoGridTechnical
using M1
using cryoGridSEB
using cryoGridSoil
using matlab
using cryoGridInitialize
using cryoGridInfiltrationUnfrozenSoil
using cryoGridSnow

#function CryoGrid3()
#workspace()
PARA = Dict()
PARA["soil"]=Dict()
PARA["snow"]=Dict()
PARA["constants"]=Dict()
PARA["location"]=Dict()
PARA["technical"]=Dict()
PARA["forcing"]=Dict()
PARA["surf"]=Dict()
PARA["modules"]=Dict()

#simple stratigraphy used to test water balance:
PARA["soil"]["layer_properties"]=[[0.0    0.4    0.25    0.00   1   0.75];
                                  [0.5    0.4    0.5     0.00   1   0.50];
                                  [10.0   0.25   0.75    0.00   1   0.25]]

# soil stratigraphy
# column 1: start depth of layer (first layer must start with 0) - each layer extends until the beginning of the next layer, the last layer
# extends until the end of the model domain
# column 2: volumetric water+ice content
# column 3: volumetric mineral content
# column 4: volumetric organic content
# column 5: code for soil type: 1: sand, 2: silt
# column 6: natural porosity - should be the same as 1-mineral-organic if no ground subsidence/thermokarst occurs


#------ model parameters --------------------------------------------------
PARA["soil"]["albedo"]=0.2;       # albedo snow-free surface
PARA["soil"]["albedoPond"]=0.07;  # albedo of water, used when the uppermost grod cell is 100# water due to modeled thermokarst development
PARA["soil"]["epsilon"]=0.97;     # emissvity snow-free surface
PARA["soil"]["z0"]=1e-3;          # roughness length [m] snow-free surface
PARA["soil"]["rs"]=50.0;          # surface resistance against evapotransiration [m^-1] snow-free surface
PARA["soil"]["Qgeo"]=0.05;        # geothermal heat flux [W/m2]
PARA["soil"]["kh_bedrock"]=3.0;   # thermal conductivity of the mineral soil fraction [W/mK]

# parameters related to hydrology scheme
PARA["soil"]["fieldCapacity"]=0.3;    #water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
PARA["soil"]["evaporationDepth"]=0.1; #depth to which evaporation occurs - place on grid cell boundaries
PARA["soil"]["rootDepth"]=0.2;        #depth affected by transpiration - place on grid cell boundaries
PARA["soil"]["wiltingPoint"]=0.2;     #point at which transpiration shuts off
PARA["soil"]["residualWC"]=0.05;      #water always remaining in the soil, not accessible to evaporation
PARA["soil"]["ratioET"]=0.5;          # 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc
PARA["soil"]["externalWaterFlux"]=0.0; # external water flux / drainage in [m/day]
PARA["soil"]["convectiveDomain"]=[];       # soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
PARA["soil"]["mobileWaterDomain"]=[];      # soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile
PARA["soil"]["waterTable"]=0.0;              # depth at which a water table will form [m] - above excess water is removed, below it pools up


# parameters related to snow
PARA["snow"]["max_albedo"]=0.85;      # albedo of fresh snow
PARA["snow"]["min_albedo"]=0.5;       # albedo of old snow
PARA["snow"]["epsilon"]=0.99;         # surface emissivity snow
PARA["snow"]["z0"]=5.0e-4;            # roughness length surface [m]
PARA["snow"]["rs"]=0.0;               # surface resistance -> should be 0 for snow
PARA["snow"]["rho_snow"]=200.0;       # density in [kg/m3]
PARA["snow"]["tau_1"]=86400.0;        # time constants of snow albedo change (according to ECMWF reanalysis) [sec]
PARA["snow"]["tau_a"]=0.008;          # [per day]
PARA["snow"]["tau_f"]=0.24;           # [per day]
PARA["snow"]["maxSnow"]=[];        # maximum snow depth that can be reached [m] - excess snow is removed in the model - if empty, no snow threshold
PARA["snow"]["extinction"]=25.0;      # light extinction coefficient of snow

PARA["location"]["altitude"]=20.0;    #used to generate pressure forcing based on barometric altitude formula, if pressure forcing is not given

PARA["technical"]["z"]=2.0;                       # height of input air temperature above ground in [m] - assumed constant even when snow depth increases
PARA["technical"]["SWEperCell"]=0.005;            # SWE per grid cell in [m] - determines size of snow grid cells
PARA["technical"]["maxSWE"]=1.0;                   # in [m] SWE
PARA["technical"]["arraySizeT"]=5002;
PARA["technical"]["starttime"]= Dates.datetime2julian(Dates.DateTime("1979.06.01 03:00:00", "Y.m.d H:M:S"));       # starttime of the simulation - if empty start from first value of time series
PARA["technical"]["endtime"] =  Dates.datetime2julian(Dates.DateTime("1983.06.01 00:00:00", "Y.m.d H:M:S"));       # endtime of the simulation - if empty end at last value of time series
PARA["technical"]["minTimestep"]=0.1 / 3600.0 / 24.0;   # smallest possible time step in [days] - here 0.1 seconds
PARA["technical"]["maxTimestep"]=300.0 / 3600.0 / 24.0;   # largest possible time step in [days] - here 300 seconds
PARA["technical"]["targetDeltaE"]=1.0e5;            # maximum energy change of a grid cell between time steps in [J/m3]
PARA["technical"]["outputTimestep"]=3.0 / 24.0;    # output time step in [days] - here three hours
PARA["technical"]["saveDate"]="01.06.";           # date of year when output file is written - no effect if "saveInterval" is empty
PARA["technical"]["saveInterval"]=1;             # interval [years] in which output files are written - if empty the entire time series is written - minimum is 1 year

#default grid used for publications and testing of water balance:
PARA["technical"]["subsurfaceGrid"] = [collect(0:0.02:2); collect(2.1:0.1:10); collect(10.2:0.2:20); collect(21:1:30); collect(35:5:50); collect(60:10:100); collect(200:100:1000)]; # the subsurface K-grid in [m]
#very simple grid used for testing of energy balance:

#initial temperature profile -> first column depth [m] -> second column temperature [degree C]
#default:
PARA["Tinitial"] = [[-5.0   10.0];
                     [0.0   10.0];
                     [5.0  -9.0];
                     [500.0  4.0];
                     [5000.0 10.0]]
#simple linear gradients

# important natural constants, given in SI units
PARA["constants"]["kappa"] = 0.4;                                     # von Kármán constant [-]
PARA["constants"]["sigma"] = 5.6704e-8;                               # Stefan-Boltzmann constant [ W / (m^2 K^4) ]
PARA["constants"]["g"] = 9.81;                                        # gravitational acceleration [m/s^2]
PARA["constants"]["p_0"] = 100500.0;                                    # normal pressure (sea level) [Pa=kg/(m s^2)]
#water
PARA["constants"]["rho_w"] = 1000.0;                                    # density of liquid water (and ice) [kg/m^3]
PARA["constants"]["c_w"] = 4200.0 * PARA["constants"]["rho_w"];         # volumetric heat capacity of water [ J / (m^3 K) ]
PARA["constants"]["k_w"] = 0.57;                                      # heat conductivity of water [ W/(mK) ] [Hillel(1982)]
#ice
PARA["constants"]["rho_i"] = 1000.0;                                    # density of ice, assumed to be equal to that of water [kg/m^3]
PARA["constants"]["c_i"] = 1900.0 * PARA["constants"]["rho_i"];         # volumetric heat capacity of ice [ J / (m^3 K) ]
PARA["constants"]["k_i"] = 2.2;                                       # heat conductivity of ice [ W/(mK) ] [Hillel(1982)]
#latent heat of water
PARA["constants"]["T_f"] = 273.15;                                    # freezing point of water / zero degree Celsius [K]
PARA["constants"]["L_sl"] = 334.0e3;                                    # specific latent heat of fusion of water [J/kg]            [AMS]
PARA["constants"]["L_lg"] = 2501.0e3;                                   # specific latent heat of vaporization of water [J/kg]      [AMS]
PARA["constants"]["L_sg"] = PARA["constants"]["L_sl"] + PARA["constants"]["L_lg"];# specific latent heat of sublimation of water [J/kg]
#air
PARA["constants"]["rho_a"] = 1.293;                                   # density of air [kg/m^3] @ 0°C, sea level
PARA["constants"]["c_a"] = 1005.7 * PARA["constants"]["rho_a"];             #c_a= 0.00125*10^6;#[J/m^3 K]   # volumetric heat capacity of dry air [J/(m^3 K)] @ 0°C, sea level, isobar
PARA["constants"]["k_a"] = 0.0243;                                    #ka=0.025; [Hillel(1982)]       # heat conductivity of air [ W/(mK)] @ 0 °C, sea level pressure
PARA["constants"]["R_a"] = 287.058;                                   # specific gas constant of air [ J/(kg K) ]
# organic
#PARA["constants"]["rho_o"] = 1; # n.a.
PARA["constants"]["c_o"] = 2.5e6; #[J/(K m^3)]                        # volumetric heat capacity of organic material [J/(K m^3)]
PARA["constants"]["k_o"] = 0.25;                                      # heat conductivity of organic material [ W/(mK) ] [Hillel(1982)]
# mineral
#PARA["constants"]["rho_m"] = 1; # n.a.
PARA["constants"]["c_m"] = 2.0e6; #[J/(K m^3)]                          # volumetric heat capacity of minearal material [J/(K m^3)]
PARA["constants"]["k_m"] = PARA["soil"]["kh_bedrock"];                  # heat conductivity of mineral material / bedrock [ W/(mK)] (specified above) #km=3.8 #mineral [Hillel(1982)]

#FORCING data mat-file
PARA["forcing"]["filename"]="CG3_CCLM_forcing_90_101.mat";  #must be in subfolder "forcing" and follow the conventions for CryoGrid 3 forcing files
PARA["forcing"]["rain_fraction"]=1.0;
PARA["forcing"]["snow_fraction"]=1.0;


PARA["modules"]["infiltration"] = true;  # true if infiltration into unfrozen ground occurs
PARA["modules"]["xice"] = false;

PARA["technical"]["run_number"] =  "test_001"

# ------make output directory (name depends on parameters) ----------------
mkpath(PARA["technical"]["run_number"])

#--------------------------------------------------------------------------
#-----------do not modify from here onwards--------------------------------
#--------------------------------------------------------------------------
FORCING, success = M1.load_forcing_from_file(PARA); # load FORCING mat-file

if success!==1
    return
end
clear!(:success)

PARA = M1.initializeParameters(PARA, FORCING); #set start time, etc.

#----------------create and initialize the grids --------------------------
GRID = cryoGridTechnical.makeGrids(PARA);              #create all grids
GRID = cryoGridSoil.createStratigraphy(PARA,GRID);     #interpolate input stratigraphy to the soil grid

#----- initializie soil thermal properties --------------------------------
GRID = cryoGridSoil.initializeSoilThermalProperties(GRID,PARA);

#------ initializie snow properties----------------------------------------
GRID = cryoGridInitialize.initializeSnow(GRID);

#---- initialize the surface energy balance struct ------------------------
SEB = cryoGridInitialize.initializeSEB();

#---- initialize temperature profile --------------------------------------
T = cryoGridInitialize.inititializeTemperatureProfile_simple(GRID, PARA, FORCING);

#---- modification for infiltration
wc=GRID["soil"]["cT_water"];
GRID["soil"]["E_lb"] = find(PARA["soil"]["evaporationDepth"].==GRID["soil"]["soilGrid"][:,1])-1;
GRID["soil"]["E_lb"] = GRID["soil"]["E_lb"][1];
GRID["soil"]["T_lb"] = find(PARA["soil"]["rootDepth"].==GRID["soil"]["soilGrid"][:,1])-1;
GRID["soil"]["T_lb"] = GRID["soil"]["T_lb"][1];

#---- preallocate temporary arrays for capacity and conductivity-----------
c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid = cryoGridInitialize.initializeConductivityCapacity(T,wc, GRID, PARA);

#---- energy and water balance initialization -----------------------------
BALANCE = cryoGridInitialize.initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA);

#__________________________________________________________________________
#-------- provide arrays for data storage ---------------------------------
t, TEMPORARY = cryoGridTechnical.generateTemporary(T, PARA);
OUT = cryoGridTechnical.generateOUT(GRID, PARA, TEMPORARY);

print("initialization successful\n");
outdict=Dict("FORCING" => FORCING, "PARA" => PARA,"GRID" => GRID)
matwrite(string(PARA["technical"]["run_number"],"/",PARA["technical"]["run_number"],"_settings.mat"), outdict)
#write(file, "test",test)
#close(file)
## ________________________________________________________________________
# Time Integration Routine                                                I
#                                                                         I
#_________________________________________________________________________I
timestep=1.0/(3600.0*24.0)*60.0;
tt=t
tic()

while t<PARA["technical"]["endtime"]


#------ interpolate forcing data to time t ----------------------------
FORCING = cryoGridTechnical.interpolateForcingData(t, FORCING);

#------determine the thermal properties of the model domains ----------
c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid = cryoGridInfiltrationUnfrozenSoil.getThermalPropertiesInfiltration(T, wc, c_cTgrid, k_cTgrid, k_Kgrid, lwc_cTgrid, GRID, PARA);
lwc = lwc_cTgrid[GRID["soil"]["cT_domain"]];

#------- water and energy balance calculations ------------------------
BALANCE = cryoGridTechnical.updateBALANCE(T, wc, c_cTgrid, lwc_cTgrid, BALANCE, GRID, PARA);

#------ surface energy balance module ---------------------------------
#set surface conditions (albedo, roughness length, etc.)
PARA, GRID = cryoGridSEB.surfaceCondition(GRID, PARA, T);

#calculate the surface energy balance
SEB, dwc_dt = cryoGridInfiltrationUnfrozenSoil.surfaceEnergyBalanceInfiltration(T, wc, FORCING, GRID, PARA, SEB);

#------ soil module  --------------------------------------------------
#calculate heat conduction
SEB = cryoGridSoil.heatConduction(T, k_Kgrid, GRID, PARA, SEB);

#------ sum up heat fluxes --------------------------------------------
SEB["dE_dt"] = SEB["dE_dt_cond"] + SEB["dE_dt_SEB"];

#------ determine optimal timestep [days] -----------------------------
# accounting for min and max timesteps specified, maximum energy change per grid cell and the CFT stability criterion
timestep = M1.timestep(t, c_cTgrid, k_cTgrid, GRID, PARA, SEB, TEMPORARY)

#------ update T array ------------------------------------------------
T = T + SEB["dE_dt"]./c_cTgrid./GRID["general"]["K_delta"].*timestep.*24.0.*3600.0;
T[GRID["air"]["cT_domain"]] = FORCING["i"]["Tair"];

#------- snow cover module --------------------------------------------
T, GRID, PARA, SEB, BALANCE = cryoGridSnow.CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_cTgrid, timestep, BALANCE);
GRID, T, BALANCE = cryoGridSnow.updateGRID_snow(T, GRID, PARA, BALANCE);


#------- infiltration module-------------------------------------------
#tic()
if PARA["modules"]["infiltration"]
    wc, GRID, BALANCE = cryoGridInfiltrationUnfrozenSoil.CryoGridInfiltration(T, wc, dwc_dt, timestep, GRID, PARA, FORCING, BALANCE);
end
#toc()

#------- update Lstar for next time step ------------------------------
SEB = cryoGridSEB.L_star(FORCING, PARA, SEB);

#------- next time step -----------------------------------------------
#print(string(timestep,"\n"))
TEMPORARY, OUT = M1.sum_up_output_store(GRID,PARA,SEB,TEMPORARY,OUT,T,t,timestep)


t=t+timestep;
end
toc()
#end
#end
