module cryoGridInitialize
    using MAT
    using matlab
    using cryoGridSoil
    using cryoGridSnow
    using cryoGridInfiltrationUnfrozenSoil

    function initializeSnow(GRID)
        GRID["snow"]["Snow_i"]=zeros(size(GRID["air"]["cT_domain"]));
        GRID["snow"]["Snow_w"]=zeros(size(GRID["air"]["cT_domain"]));
        GRID["snow"]["Snow_a"]=zeros(size(GRID["air"]["cT_domain"]));
        GRID["snow"]["SWEinitial"]=0.0;

        return GRID
    end

    function initializeSEB()
        size_L_star_smoothing=1;

        SEB=Dict();
        SEB["Qnet"]=0.0;
        SEB["Qh"]=0.0;
        SEB["Qe"]=0.0;
        SEB["Qg"]=0.0;
        SEB["Sout"]=0.0;
        SEB["Lout"]=0.0;
        SEB["newSnow"]=0.0;
        #SEB.sublim=0;
        #SEB.meltwater=0;
        SEB["L_star"]=-100000.0#+zeros(1,size_L_star_smoothing);
        SEB["u_star"]=10.0;


        SEB["Qsurf"] = 0.0;  # for EB checks
        return SEB
    end

    function inititializeTemperatureProfile_simple(GRID, PARA, FORCING)
        T = zeros(size(GRID["general"]["cT_grid"]));

        idx=indmin(abs.(FORCING["data"]["t_span"]-PARA["technical"]["starttime"]));
        T[GRID["air"]["cT_domain"]]   = FORCING["data"]["Tair"][idx];
        T[GRID["snow"]["cT_domain"]]  = FORCING["data"]["Tair"][idx];
        T[GRID["soil"]["cT_domain"]]  = matlab.interp1(PARA["Tinitial"][:,1], PARA["Tinitial"][:,2], GRID["general"]["cT_grid"][GRID["soil"]["cT_domain"]], "linear");
        return T
    end

    function initializeConductivityCapacity(T, wc, GRID, PARA)
        c_temp = zeros(size(GRID["general"]["cT_grid"]));
        k_temp = zeros(size(GRID["general"]["cT_grid"]));
        k_eff = zeros(size(GRID["general"]["K_grid"]));
        lwc_temp = zeros(size(GRID["general"]["cT_grid"]));


        #------- unused grid cells --------------------------------
        c_temp[GRID["air"]["cT_domain"]] = PARA["constants"]["c_a"];    #set some value e.g. air
        k_temp[GRID["air"]["cT_domain"]] = PARA["constants"]["k_a"];    #set some value e.g. air
        lwc_temp[GRID["air"]["cT_domain"]] = 0;

        #------- soil domain --------------------------------------
        c_temp[GRID["soil"]["cT_domain"]],
        k_temp[GRID["soil"]["cT_domain"]],
        lwc_temp[GRID["soil"]["cT_domain"]] = cryoGridSoil.readThermalParameters(T[GRID["soil"]["cT_domain"]], GRID, PARA);

        #adjust for the unfrozen part of the domain
        #JAN: this changes with infiltration scheme: frozen  (T<=0) remains unchanged,
        # #thawed (T>0) is calculated differently in dependence of wc
        c_temp[GRID["soil"]["cT_domain"]] = (T[GRID["soil"]["cT_domain"]].<=0).*c_temp[GRID["soil"]["cT_domain"]] + (T[GRID["soil"]["cT_domain"]].>0).* cryoGridInfiltrationUnfrozenSoil.capacityUnfrozen(wc,GRID,PARA);
        k_temp[GRID["soil"]["cT_domain"]] = (T[GRID["soil"]["cT_domain"]].<=0).*k_temp[GRID["soil"]["cT_domain"]] + (T[GRID["soil"]["cT_domain"]].>0).* cryoGridInfiltrationUnfrozenSoil.conductivityUnfrozen(wc,GRID,PARA);
        lwc_temp[GRID["soil"]["cT_domain"]] = (T[GRID["soil"]["cT_domain"]].<=0).*lwc_temp[GRID["soil"]["cT_domain"]] + (T[GRID["soil"]["cT_domain"]].>0).* wc;

        #------- snow domain --------------------------------------
        c_temp[GRID["snow"]["cT_domain"]] = cryoGridSnow.cap_snow(GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain"]],
        GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain"]],
        GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain"]],
        PARA);

        k_temp[GRID["snow"]["cT_domain"]] = cryoGridSnow.cond_snow(GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain"]],
        GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain"]],
        GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain"]]);

        lwc_temp[GRID["snow"]["cT_domain"]] = GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain"]]./GRID["general"]["K_delta"][GRID["snow"]["cT_domain"]];
        return c_temp, k_temp, k_eff, lwc_temp
    end

    function initializeBALANCE(T, wc, c_cTgrid, lwc_cTgrid, GRID, PARA)
        # ENERGY balance
        # energy content soil domain in [J/m^2] distinguished by sensible and latent part
        BALANCE=Dict();
        BALANCE["energy"]=Dict();
        BALANCE["water"]=Dict();

        E_soil_sens = (PARA["constants"]["c_w"].*lwc_cTgrid[GRID["soil"]["cT_domain"]]+
            PARA["constants"]["c_i"].*(wc-lwc_cTgrid[GRID["soil"]["cT_domain"]]) +
            PARA["constants"]["c_m"].*GRID["soil"]["cT_mineral"]+PARA["constants"]["c_o"].*GRID["soil"]["cT_organic"]).*
            T[GRID["soil"]["cT_domain"]].* GRID["general"]["K_delta"][GRID["soil"]["cT_domain"]];
        BALANCE["energy"]["E_soil_sens"]=sum(E_soil_sens[.!isnan.(E_soil_sens)])

        E_soil_lat = PARA["constants"]["rho_w"] .* PARA["constants"]["L_sl"] .* lwc_cTgrid[GRID["soil"]["cT_domain"]].*
            GRID["general"]["K_delta"][GRID["soil"]["cT_domain"]];
        BALANCE["energy"]["E_soil_lat"]=sum(E_soil_lat[.!isnan.(E_soil_lat)])

        BALANCE["energy"]["E_soil"] = BALANCE["energy"]["E_soil_sens"] + BALANCE["energy"]["E_soil_lat"];
        # energy content snow domain in [J/m^2] distinguished by sensible and latent part
        E_snow_sens = c_cTgrid[GRID["snow"]["cT_domain"]] .* T[GRID["snow"]["cT_domain"]] .* GRID["general"]["K_delta"][GRID["snow"]["cT_domain"]];
        BALANCE["energy"]["E_snow_sens"] = sum(E_snow_sens[.!isnan.(E_snow_sens)])

        E_snow_lat = PARA["constants"]["rho_w"] .* PARA["constants"]["L_sl"] .* lwc_cTgrid[GRID["snow"]["cT_domain"]].* GRID["general"]["K_delta"][GRID["snow"]["cT_domain"]];
        BALANCE["energy"]["E_snow_lat"] = sum(E_snow_lat[.!isnan.(E_snow_lat)])

        BALANCE["energy"]["E_snow"] = BALANCE["energy"]["E_snow_sens"] + BALANCE["energy"]["E_snow_lat"];
        # accumulated changes per output timestep
        BALANCE["energy"]["dE_soil_sens"] = 0.0;
        BALANCE["energy"]["dE_soil_lat"] = 0.0;
        BALANCE["energy"]["dE_soil"] = 0.0;
        BALANCE["energy"]["dE_snow_sens"] = 0.0;
        BALANCE["energy"]["dE_snow_lat"] = 0.0;
        BALANCE["energy"]["dE_snow"] = 0.0;

        # WATER balance
        # water content soil domain in [m]
        BALANCE["water"]["W_soil"] = sum(wc .* GRID["general"]["K_delta"][GRID["soil"]["cT_domain"]]);

        # water content snow domain in [m]
        BALANCE["water"]["W_snow"] = sum(GRID["snow"]["Snow_i"] + GRID["snow"]["Snow_w"]) + GRID["snow"]["SWEinitial"];
        # accumulated changes per output timestep
        # storage
        BALANCE["water"]["dW_soil"] = 0.0;
        BALANCE["water"]["dW_snow"] = 0.0;
        # precipitation
        BALANCE["water"]["dp_rain"]=0.0;
        BALANCE["water"]["dp_snow"]=0.0; # SWE
        # evapotranspiration and sublimation
        BALANCE["water"]["de"]=0.0;
        BALANCE["water"]["ds"]=0.0;
        # runoff
        BALANCE["water"]["dr_surface"]=0.0;
        BALANCE["water"]["dr_subsurface"]=0.0;
        BALANCE["water"]["dr_snowmelt"]=0.0;
        BALANCE["water"]["dr_excessSnow"]=0.0;
        BALANCE["water"]["dr_rain"]=0.0;  # this is only rain on frozen ground
        return BALANCE
    end

end
