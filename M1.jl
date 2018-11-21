module M1
    using MAT
    using cryoGridTechnical

    function load_forcing_from_file(PARA)

        var=matread(string("forcing/",PARA["forcing"]["filename"]))
        FORCING=var["FORCING"]

        #convert to julian
        FORCING["data"]["t_span"]=FORCING["data"]["t_span"] - 1 + Dates.datetime2julian(Dates.DateTime("0000.01.01 00:00:00", "Y.m.d H:M:S"))

        FORCING["data"]["rainfall"]=FORCING["data"]["rainfall"]*PARA["forcing"]["rain_fraction"];
        FORCING["data"]["snowfall"]=FORCING["data"]["snowfall"]*PARA["forcing"]["snow_fraction"];

        if std(FORCING["data"]["t_span"][2:end,1]-FORCING["data"]["t_span"][1:end-1,1])>1e-8
            print("timestamp of forcing data is not in regular intervals -> check, fix and restart")
            success=0;
        end

        #here, consistency checks, RH->q calculation, set threhsolds for wind, etc"] could be placed

        FORCING["data"]["wind"][FORCING["data"]["wind"].<0.5]=0.5; #set min wind speed to 0.5 m/sec to avoid breakdown of turbulence

        #set pressure to mean pressure at corresponding altitude (international
        #altitude formula) if now provided by the forcing time series
        if !haskey(FORCING["data"], "p")
            FORCING["data"]["p"]=FORCING["data"]["Tair"]*0.0 + 1013.25*100.0*(1.0-0.0065/288.15*PARA["location"]["altitude"])^5.255
        end

        success=1;

        #initialize
        FORCING["i"]=Dict()
        FORCING["i"]["snowfall"]=0.0;
        FORCING["i"]["rainfall"]=0.0;
        FORCING["i"]["Lin"]=0.0;
        FORCING["i"]["Sin"]=0.0;
        FORCING["i"]["Tair"]=0.0;
        FORCING["i"]["wind"]=0.0;
        FORCING["i"]["RH"]=0.0;
        FORCING["i"]["q"]=0.0;
        FORCING["i"]["p"]=0.0;

        return FORCING, success
    end

    function initializeParameters(PARA, FORCING)

        PARA["surf"]["albedo"]=PARA["soil"]["albedo"];
        PARA["surf"]["epsilon"]=PARA["soil"]["epsilon"];
        PARA["surf"]["z0"]=PARA["soil"]["z0"];
        PARA["surf"]["rs"]=PARA["soil"]["rs"];
        PARA["snow"]["albedo"]=PARA["snow"]["max_albedo"]; # sets the initial albedo

        if isempty(PARA["technical"]["starttime"])
            PARA["technical"]["starttime"]=FORCING["data"]["t_span"][1,1];
        end

        if isempty(PARA["technical"]["endtime"])
            PARA["technical"]["endtime"]=FORCING["data"]["t_span"][end-1,1]; #be on the save side end-1
        end
        return PARA
    end

    function timestep(t, c_cTgrid, k_cTgrid, GRID, PARA, SEB, TEMPORARY)
        technical_maxTimestep=PARA["technical"]["maxTimestep"]::Float64
        technical_minTimestep=PARA["technical"]["minTimestep"]::Float64
        technical_targetDeltaE=PARA["technical"]["targetDeltaE"]::Float64
        general_K_delta=GRID["general"]["K_delta"]::Array{Float64,1}
        soil_cT_domain=GRID["soil"]["cT_domain"]::Array{Bool,1}
        snow_cT_domain=GRID["snow"]["cT_domain"]::Array{Bool,1}
        dE_dt=SEB["dE_dt"]::Array{Float64,1}

        soilsnow_domain=soil_cT_domain .| snow_cT_domain;

        timestep = minimum([maximum([minimum([0.5*minimum(general_K_delta[soilsnow_domain].^2.0.*c_cTgrid[soilsnow_domain]./k_cTgrid[soilsnow_domain])./(24.0.*3600.0),
                                            technical_targetDeltaE.* minimum(abs.(general_K_delta./dE_dt))./(24.0.*3600.0),
                                            technical_maxTimestep]),
                                            technical_minTimestep]),
                                            TEMPORARY["outputTime"]-t]);
        return timestep
    end

    function sum_up_output_store(GRID,PARA, SEB, TEMPORARY, OUT, T, t, timestep)
        TEMPORARY["timestep_sum"]=TEMPORARY["timestep_sum"]+(timestep*24.0*3600.0)*timestep;
        TEMPORARY["T_sum"]=TEMPORARY["T_sum"]+T.*timestep;
        TEMPORARY["Qe_sum"]=TEMPORARY["Qe_sum"]+SEB["Qe"].*timestep;
        TEMPORARY["Qh_sum"]=TEMPORARY["Qh_sum"]+SEB["Qh"].*timestep;
        TEMPORARY["Qnet_sum"]=TEMPORARY["Qnet_sum"]+SEB["Qnet"].*timestep;
        TEMPORARY["Qg_sum"]=TEMPORARY["Qg_sum"]+SEB["Qg"].*timestep;

        #----store in output table --------------------------------------------
        if  t==TEMPORARY["outputTime"]

            TEMPORARY["counter"] = TEMPORARY["counter"]+1;

            #average over timespan:
            TEMPORARY["dt_out"]=t-TEMPORARY["t_last"];
            TEMPORARY["t_last"]=t;

            TEMPORARY["T_out"]=TEMPORARY["T_sum"]./TEMPORARY["dt_out"];


            TEMPORARY["Qh"]=TEMPORARY["Qh_sum"]./TEMPORARY["dt_out"];
            TEMPORARY["Qe"]=TEMPORARY["Qe_sum"]./TEMPORARY["dt_out"];
            TEMPORARY["Qnet"]=TEMPORARY["Qnet_sum"]./TEMPORARY["dt_out"];
            TEMPORARY["Qg"]=TEMPORARY["Qg_sum"]./TEMPORARY["dt_out"];

            TEMPORARY["timestep_out"]=TEMPORARY["timestep_sum"]./TEMPORARY["dt_out"];

            #store new values in OUT struct -----------------------------------
            nan_mask_air=GRID["air"]["cT_domain"]*NaN+1.0
            OUT["cryoGrid3"][:,TEMPORARY["counter"]] = TEMPORARY["T_out"].*nan_mask_air;


            #reset sum variables
            TEMPORARY["T_sum"]=zeros(size(T));

            TEMPORARY["Qh_sum"]=0.0;
            TEMPORARY["Qe_sum"]=0.0;
            TEMPORARY["Qnet_sum"]=0.0;
            TEMPORARY["Qg_sum"]=0.0;

            TEMPORARY["timestep_sum"]=0.0;


            print(string(now(),":  at ",Dates.julian2datetime(t), ",  Average timestep: ",  string(TEMPORARY["timestep_out"]), " seconds","\n"))

            TEMPORARY["outputTime"]=round((TEMPORARY["outputTime"]+PARA["technical"]["outputTimestep"])./PARA["technical"]["outputTimestep"]).*PARA["technical"]["outputTimestep"];
            #write output files
            if  round((t-TEMPORARY["saveTime"]).*48)==0

                outdict=Dict("OUT" => OUT)
                matwrite(string(PARA["technical"]["run_number"],"/",PARA["technical"]["run_number"],"_output",Dates.year(Dates.julian2datetime(t)),".mat"), outdict)
                TEMPORARY["saveTime"]=Dates.datetime2julian(Dates.DateTime(Dates.year(Dates.julian2datetime(t))+1,
                                                                            Dates.month(Dates.julian2datetime(t)),
                                                                            Dates.day(Dates.julian2datetime(t)),
                                                                            Dates.hour(Dates.julian2datetime(t)),
                                                                            Dates.minute(Dates.julian2datetime(t)),0));

                OUT = cryoGridTechnical.generateOUT(GRID, PARA, TEMPORARY);
                TEMPORARY["counter"]=0;
            end

        end

        return TEMPORARY, OUT
    end

end
