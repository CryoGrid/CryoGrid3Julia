module cryoGridSnow
    using MAT
    using cryoGridTechnical

    function cond_snow(snow_i, snow_w, snow_a)
        total=snow_w + snow_i + snow_a;
        #cond_snow=conductivity2(snow_w./total, snow_i./total, 0, 0, 2);
        cond_snow = 2.2.*((snow_w+snow_i)./total).^1.88; # Yen (1981)
#        if length(cond_snow)>1
#            cond_snow[isnan.(cond_snow)]=0.3;
#        else
#            if isnan(cond_snow)
#                cond_snow=NaN
#            end
#        end
        return cond_snow
    end

    function cap_snow(Snow_i, Snow_w, Snow_a, PARA)
        c_i = PARA["constants"]["c_i"]::Float64; #1.9*10^6;#[J/m�K]
        c_w = PARA["constants"]["c_w"]::Float64; #4.2*10^6; #[J/m�K]

        cap_snow = (Snow_i.* c_i + Snow_w.* c_w)./ (Snow_i + Snow_w + Snow_a);  # JAN: why not also add Snow_a.*c_a
#        if length(cap_snow)>1
#            cap_snow[isnan.(cap_snow)]=0.3;
#        else
#            if isnan(cap_snow)
#                cap_snow=NaN
#            end
#        end
        return cap_snow
    end

    function CryoGridSnow(T, GRID, FORCING, SEB, PARA, c_temp, timestep, BALANCE)

        if !isempty(GRID["snow"]["cT_domain_ub"]) #snow cover already exitis


            #----------calculate snow surface albedo ------------------------------

            if maximum(T[GRID["snow"]["cT_domain"]])>=0.0    # melting conditions
                PARA["snow"]["albedo"] = PARA["snow"]["min_albedo"] + (PARA["snow"]["albedo"] - PARA["snow"]["min_albedo"]) .*
                                        exp(-PARA["snow"]["tau_f"] .* timestep.*24.0.*3600.0 ./ PARA["snow"]["tau_1"]);
            else
                PARA["snow"]["albedo"]=max(PARA["snow"]["albedo"]-PARA["snow"]["tau_a"].*timestep.*24.0.*3600.0./PARA["snow"]["tau_1"], PARA["snow"]["min_albedo"]);
            end

            #--------SEB.sublimation-----------------------------------------------
            delta_SnowSub=SEB["Qe"].*timestep.*24.0.*3600.0./PARA["constants"]["L_sg"]./PARA["constants"]["rho_i"]

            GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] = GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] - delta_SnowSub;#- SEB["Qe"].*timestep.*24.*3600./(L+L_lv)./1000;

            nonAirFractionUppermostGridCell = (GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]]+GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]])./
            (GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]]+GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]]+GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]]);

            GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]] = GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]] - (delta_SnowSub./nonAirFractionUppermostGridCell -  delta_SnowSub); #-  SEB["Qe"].*timestep.*24.*3600./(L+L_lv)./1000);


            BALANCE["water"]["ds"] = BALANCE["water"]["ds"] - delta_SnowSub*1000.0; # sublimation in [mm]

            #---------- melt and infiltration -------------------------------------
            if maximum(T[GRID["snow"]["cT_domain"]])>0.0 || FORCING["i"]["rainfall"]>0.0 || sum(GRID["snow"]["Snow_w"])>0.0  #cases when melt or infiltration occur

                T[GRID["snow"]["cT_domain"]],
                GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain"]],
                GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain"]],
                GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain"]],
                newMelt =
                cryoGridSnow.snowMelt(T[GRID["snow"]["cT_domain"]],
                                        GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain"]],
                                        GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain"]],
                                        GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain"]],
                                        FORCING["i"]["rainfall"].*timestep./1000.0,
                                        c_temp[GRID["snow"]["cT_domain"]],
                                        PARA);

                #account for meltwater in water balance
                BALANCE["water"]["dr_snowmelt"] = BALANCE["water"]["dr_snowmelt"] + (-newMelt.*1000.0);    # in [mm]
            end

            #-------- add the new snow to the upper most snow cell in the case of a exisiting snow cover -------------------
            if isempty(PARA["snow"]["maxSnow"])
                deltaSnow_i = max(0.0, FORCING["i"]["snowfall"].*timestep./1000.0);
            else
                snowHeight = abs( GRID["general"]["K_grid"][GRID["snow"]["cT_domain_ub"]] - GRID["general"]["K_grid"][GRID["snow"]["cT_domain_lb"]+1] );
                if snowHeight>PARA["snow"]["maxSnow"]
                    print(" excess snow occurs ");
                end
                deltaSnow_i = max(0.0,
                              min(FORCING["i"]["snowfall"].*timestep./1000.0,
                                 (PARA["snow"]["maxSnow"] - snowHeight ) .* PARA["snow"]["rho_snow"] ./ PARA["constants"]["rho_w"])); #ensures that no more than maxSnow can accumulate
                #account for excess snow in water balance
                BALANCE["water"]["dr_excessSnow"] = BALANCE["water"]["dr_excessSnow"] -(FORCING["i"]["snowfall"].*timestep - deltaSnow_i*1000.0); #defined as negative when snow is removed; in [mm]
            end

            GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] = GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] + deltaSnow_i;

            GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]] = GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]] + (deltaSnow_i./(PARA["snow"]["rho_snow"]./1000.0) - deltaSnow_i);


        else #no snow cover

            #---------- add the new snow into initial SWE variable in case of no snow cover------------------

            GRID["snow"]["SWEinitial"] = GRID["snow"]["SWEinitial"] + FORCING["i"]["snowfall"].*timestep./1000.0 - GRID["snow"]["SWEinitial"].*0.1.*timestep;
            # account for decrease in SWEinitial in water balance
            BALANCE["water"]["dr_snowmelt"] = BALANCE["water"]["dr_snowmelt"] - GRID["snow"]["SWEinitial"].*0.1.*timestep*1000.0; #SWEinitial decreasing counted as surface runoff

            #----- add the rainfall as runoff in case of no infiltration into frozen ground
            if !PARA["modules"]["infiltration"] || (PARA["modules"]["infiltration"] && T[GRID["soil"]["cT_domain_ub"]]<=0.0)#no infiltration scheme or uppermost soil cell frozen
                BALANCE["water"]["dr_rain"] = BALANCE["water"]["dr_rain"] - FORCING["i"]["rainfall"].*timestep;
            end


        end

        #--------- update albedo after fresh fallen snow --------------------------
        # determine time of last snowfall for albedo calculation
        SEB["newSnow"] =  SEB["newSnow"]-SEB["newSnow"].*0.1.*timestep + FORCING["i"]["snowfall"].*timestep./1000.0;
        if SEB["newSnow"]>= PARA["technical"]["SWEperCell"]/2.0
            PARA["snow"]["albedo"]=PARA["snow"]["max_albedo"];
            SEB["newSnow"]=0.0;
        end

        return T, GRID, PARA, SEB, BALANCE
    end

    function snowMelt(T::Array{Float64,1}, snow_i::Array{Float64,1}, snow_w::Array{Float64,1}, snow_a::Array{Float64,1}, water_flux::Float64, c_temp::Array{Float64,1}, PARA)

        runoff=0.0;

        #------------melt the snow for cells with T>0------------
        poreSpace=(snow_w+snow_a)./(snow_i+snow_w+snow_a);
        #poreSpace[isnan.(poreSpace)]=300.0/1000.0;
        T, snow_i, snow_w, snow_a = cryoGridSnow.melt(T, snow_i, snow_w, snow_a, poreSpace, c_temp, PARA);

        pS=(snow_w+snow_a)./(snow_i+snow_w+snow_a);
        #pS[isnan.(pS)]=300.0/1000.0;
        pS_idx=pS.>0.0

        #-----------calculate how much water (in m) a snow cell can hold-----
        maxWater = cryoGridSnow.maxLiqWater(T[pS_idx], snow_i[pS_idx], snow_w[pS_idx], snow_a[pS_idx], poreSpace[pS_idx], c_temp[pS_idx], PARA);

        if water_flux>0.0 || (!isempty(maxWater) && minimum(maxWater)<0.0)   #infiltration occurs

            #------------infiltrate from top to bottom--------------------
            snow_w[pS_idx], snow_a[pS_idx], water_flux = cryoGridSnow.infiltrateTop2Bottom(snow_i[pS_idx], snow_w[pS_idx], snow_a[pS_idx], poreSpace[pS_idx], maxWater, water_flux);

            #------------infiltrate bottom to top----------------------
            if water_flux>0.0
                snow_w[pS_idx], snow_a[pS_idx], runoff = cryoGridSnow.infiltrateBottom2Top(snow_i[pS_idx], snow_w[pS_idx], snow_a[pS_idx], water_flux);
            end
        end

        #----------shift energy from water to T due to refreezing---------
        T, snow_i, snow_w = cryoGridSnow.refreeze(T, snow_i, snow_w, snow_a, c_temp, PARA);
        return T, snow_i, snow_w, snow_a, runoff
    end

    function melt(T::Array{Float64,1}, snow_i::Array{Float64,1}, snow_w::Array{Float64,1}, snow_a::Array{Float64,1}, poreSpace::Array{Float64,1}, c_temp::Array{Float64,1}, PARA)
        L=PARA["constants"]["L_sl"].*PARA["constants"]["rho_w"]; #3.34e8;

        pot_SWE=zeros(size(T))
        delta_SWE=zeros(size(T))

        @inbounds @fastmath for i=1:length(T)
            pot_SWE[i]=(T[i].>0.0).*T[i].*c_temp[i].*(snow_i[i]+snow_w[i]+snow_a[i])./L;
            T[i]=(T[i].<=0.0).*T[i];
            delta_SWE[i]=pot_SWE[i];
        end
        #if sum(pot_SWE.>snow_i)!=0.0# in one cell more energy  than needed to melt entire cell
        @inbounds @fastmath for i=1:length(T)-1
            delta_SWE[i]=min(snow_i[i], pot_SWE[i]);
            SWEres=pot_SWE[i] - delta_SWE[i];
            if (snow_i[i+1]+snow_w[i+1]+snow_a[i+1])!=0.0
                T[i+1]=T[i+1] + SWEres.*L./(c_temp[i+1].*(snow_i[i+1]+snow_w[i+1]+snow_a[i+1]));
            end
            pot_SWE[i+1]=pot_SWE[i+1] + (T[i+1]>0.0).*T[i+1].*c_temp[i+1].*(snow_i[i+1]+snow_w[i+1]+snow_a[i+1])./L;
            T[i+1]=(T[i+1]<=0.0).*T[i+1];
        end
        delta_SWE[end]=min(snow_i[end], pot_SWE[end]);
        #end
        #T(find(isnan(T(:,1))==1),1)=0;
        #delta_SWE
        # delta_SWE = min([snow_i  double(T>0).*T.*c_temp.*(snow_i+snow_w+snow_a)./L]');   #Energy conserving since sensible heat term is only excess energy from SEB+conduction
        # delta_SWE=delta_SWE';

        snow_i=snow_i - delta_SWE; #melting only changes SWE, not density theta_s
        snow_w=snow_w + delta_SWE;
#        if minimum(snow_i_out)<0.0
#            print(minimum(snow_i), minimum(snow_i_out) ,"\n")
#            ede
#        end
        #snow_a=(poreSpace.*snow_i + poreSpace.*snow_w - snow_w)./(1-poreSpace);  #pore space stays stays constant
        return T, snow_i, snow_w, snow_a
    end

    function maxLiqWater(T::Array{Float64,1}, snow_i::Array{Float64,1}, snow_w::Array{Float64,1}, snow_a::Array{Float64,1}, poreSpace::Array{Float64,1}, c_temp::Array{Float64,1}, PARA)
        L=PARA["constants"]["L_sl"].*PARA["constants"]["rho_w"]; #3.34e8;
        maxLiqWater=zeros(size(T))
        @inbounds @fastmath for i=1:length(T)
            waterHoldingPot=0.05.* (poreSpace[i].*snow_i[i])./(1.0-poreSpace[i]);   #in m; 5# of the pore space can be filled by water
            maxLiqWater_ = (waterHoldingPot[i] - snow_w[i] - T[i].*c_temp[i].*(snow_i[i]+snow_w[i]+snow_a[i])./L);  #in m
            maxLiqWater[i] = min(snow_a[i], maxLiqWater_);
        end
        #maxLiqWater = maxLiqWater'; # negative if snow_w>waterHoldingPot
        return maxLiqWater
    end

    function infiltrateBottom2Top(snow_i::Array{Float64,1}, snow_w::Array{Float64,1}, snow_a::Array{Float64,1}, waterFlux::Float64)
        j=size(snow_w,1);
        @inbounds @fastmath while waterFlux>0.0 && j>=1
            delta_SWE = ((waterFlux - snow_a[j]).>0.0).*snow_a[j] + ((waterFlux - snow_a[j]).<=0.0).*waterFlux;
            snow_w[j] = snow_w[j] + delta_SWE;
            snow_a[j] = snow_a[j] - delta_SWE;
            waterFlux=waterFlux - delta_SWE;
            j=j-1;
        end
        return snow_w, snow_a, waterFlux
    end

    function infiltrateTop2Bottom(snow_i::Array{Float64,1}, snow_w::Array{Float64,1}, snow_a::Array{Float64,1}, poreSpace::Array{Float64,1}, maxLiqWater::Array{Float64,1}, waterFlux::Float64)
        @inbounds @fastmath for j=1:size(snow_w,1)
            delta_SWE=((waterFlux - maxLiqWater[j]).>0.0).*maxLiqWater[j] + ((waterFlux - maxLiqWater[j]).<=0.0).*waterFlux;
            snow_w[j] = snow_w[j] + delta_SWE;
            if delta_SWE>=0.0
                snow_a[j] = snow_a[j] - delta_SWE;
            else
                snow_a[j]=(poreSpace[j].*snow_i[j] + poreSpace[j].*snow_w[j] - snow_w[j])./(1.0-poreSpace[j]);
            end
            waterFlux=((waterFlux - maxLiqWater[j])>0.0).*(waterFlux - maxLiqWater[j]);
        end
        return snow_w, snow_a, waterFlux
    end

    function refreeze(T::Array{Float64,1}, snow_i::Array{Float64,1}, snow_w::Array{Float64,1}, snow_a::Array{Float64,1}, c_temp::Array{Float64,1}, PARA)
        L=PARA["constants"]["L_sl"].*PARA["constants"]["rho_w"];#3.34e8;
        @inbounds @fastmath for i=1:length(T)
            delta_SWE = min(snow_w[i],  -1.0.*(T[i].<=0.0).*T[i].*c_temp[i].*(snow_i[i]+snow_w[i]+snow_a[i])./L);
            #delta_SWE=delta_SWE';
            snow_i[i]=snow_i[i] + delta_SWE;
            snow_w[i]=snow_w[i] - delta_SWE;
            T[i] = T[i] + delta_SWE.*L ./ ((snow_i[i]+snow_w[i]+snow_a[i]).*c_temp[i]);
            if isnan(T[i])
                T[i]=0.0
            end
        end
        #T[isnan.(T)]=0.0;
        return T, snow_i, snow_w
    end

    function updateGRID_snow(T, GRID, PARA, BALANCE)
        snowCellSize=GRID["snow"]["snowCellSize"];
        check_change=false;
        check_minor_change=false;

        if isempty(GRID["snow"]["cT_domain_lb"]) #no snow exists

            if GRID["snow"]["SWEinitial"]>=0.5*PARA["technical"]["SWEperCell"]   #create and initialize first cell of snow
                #------ modify snow and air grid -----------------------------
                GRID["snow"]["cT_domain"][GRID["air"]["cT_domain_lb"]]=true;
                GRID["snow"]["K_domain"][GRID["air"]["K_domain_lb"]]=true;
                GRID["snow"]["cT_domain_lb"], GRID["snow"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["cT_domain"]);
                GRID["snow"]["K_domain_lb"], GRID["snow"]["K_domain_ub"] =   cryoGridTechnical.LayerIndex(GRID["snow"]["K_domain"]);

                GRID["air"]["cT_domain"][GRID["air"]["cT_domain_lb"]]=false;
                GRID["air"]["K_domain"][GRID["air"]["K_domain_lb"]]=false;
                GRID["air"]["cT_domain_lb"], GRID["air"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["cT_domain"]);
                GRID["air"]["K_domain_lb"], GRID["air"]["K_domain_ub"] =   cryoGridTechnical.LayerIndex(GRID["air"]["K_domain"]);

                # ------- update SWE grid -------------------------------------
                GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] = GRID["snow"]["SWEinitial"];
                GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]] = 0.0;
                GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]] = (GRID["snow"]["SWEinitial"]./(PARA["snow"]["rho_snow"]./1000.0) - GRID["snow"]["SWEinitial"]);
                GRID["snow"]["SWEinitial"]=0.0;

                # -------- update K grid -------------------------------------
                GRID["general"]["K_grid"][GRID["snow"]["K_domain_ub"]]=-1.0.*( GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] + GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]] + GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]] );
                T[GRID["snow"]["cT_domain_ub"]]=T[GRID["air"]["cT_domain_lb"]];

                check_change=true;
            end

        else   #snow exists

            GRID["general"]["K_grid"][GRID["snow"]["cT_domain_ub"]] = GRID["general"]["K_grid"][GRID["snow"]["cT_domain_ub"]+1] -
            (GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] + GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]] + GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]]); #updates the position of the uppermost snow grid cell
            check_minor_change=true;

            if GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]]>=1.5.*PARA["technical"]["SWEperCell"]  #create new grid cell      # JAN: why 1.5 ???

                #------ modify snow and air grid -----------------------------
                GRID["snow"]["cT_domain"][GRID["air"]["cT_domain_lb"]]=true;
                GRID["snow"]["K_domain"][GRID["air"]["K_domain_lb"]]=true;
                GRID["snow"]["cT_domain_lb"], GRID["snow"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["cT_domain"]);
                GRID["snow"]["K_domain_lb"], GRID["snow"]["K_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["K_domain"]);

                GRID["air"]["cT_domain"][GRID["air"]["cT_domain_lb"]]=false;
                GRID["air"]["K_domain"][GRID["air"]["K_domain_lb"]]=false;
                GRID["air"]["cT_domain_lb"], GRID["air"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["cT_domain"]);
                GRID["air"]["K_domain_lb"], GRID["air"]["K_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["K_domain"]);

                # ------- update SWE grid
                GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]]=1.0./3.0.*GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]+1];
                GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]]=1.0./3.0.*GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]+1];
                GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]]=1./3.*GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]+1];
                GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]+1]=GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]+1] - GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]];
                GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]+1]=GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]+1] - GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]];
                GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]+1]=GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]+1] - GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_ub"]];
                T[GRID["snow"]["cT_domain_ub"]]=T[GRID["snow"]["cT_domain_ub"]+1];

                check_change=true;
            end

            if minimum(GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain"]])<=0.5.*PARA["technical"]["SWEperCell"] #avoid looping when unnecessary
                for i=GRID["snow"]["cT_domain_ub"]:GRID["snow"]["cT_domain_lb"]-1  #check all snow cells except for the lowermost one for too small ice and water contents - merge with lower cell

                    if GRID["snow"]["Snow_i"][i]<=0.5.*PARA["technical"]["SWEperCell"]

                        #------ modify snow and air grid -----------------------------
                        GRID["snow"]["cT_domain"][GRID["snow"]["cT_domain_ub"]]=false;
                        GRID["snow"]["K_domain"][GRID["snow"]["K_domain_ub"]]=false;
                        GRID["snow"]["cT_domain_lb"], GRID["snow"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["cT_domain"]);
                        GRID["snow"]["K_domain_lb"], GRID["snow"]["K_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["K_domain"]);

                        GRID["air"]["cT_domain"][GRID["air"]["cT_domain_lb"]+1]=true;
                        GRID["air"]["K_domain"][GRID["air"]["K_domain_lb"]+1]=true;
                        GRID["air"]["cT_domain_lb"], GRID["air"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["cT_domain"]);
                        GRID["air"]["K_domain_lb"], GRID["air"]["K_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["K_domain"]);

                        #-------- rearrange SWE and T grids --------------------------
                        GRID["snow"]["Snow_i"][i+1]=GRID["snow"]["Snow_i"][i+1]+GRID["snow"]["Snow_i"][i];
                        GRID["snow"]["Snow_w"][i+1]=GRID["snow"]["Snow_w"][i+1]+GRID["snow"]["Snow_w"][i];
                        GRID["snow"]["Snow_a"][i+1]=GRID["snow"]["Snow_a"][i+1]+GRID["snow"]["Snow_a"][i];
                        GRID["snow"]["Snow_i"][2:i]=GRID["snow"]["Snow_i"][1:i-1];
                        GRID["snow"]["Snow_w"][2:i]=GRID["snow"]["Snow_w"][1:i-1];
                        GRID["snow"]["Snow_a"][2:i]=GRID["snow"]["Snow_a"][1:i-1];
                        T[i+1]=(T[i+1]+T[i])/2.0;
                        T[2:i]=T[1:i-1];
                    end

                end
                check_change=true;
            end

            if GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_lb"]]<=0.5.*PARA["technical"]["SWEperCell"] && sum(GRID["snow"]["cT_domain"])>=2 #lowermost grid cell has too little snow, but there still is 2 or more snow cells - merge with upper cell

                #------ modify snow and air grid ----------------------------------
                GRID["snow"]["cT_domain"][GRID["snow"]["cT_domain_ub"]]=false;
                GRID["snow"]["K_domain"][GRID["snow"]["cT_domain_ub"]]=false;
                GRID["snow"]["cT_domain_lb"], GRID["snow"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["cT_domain"]);
                GRID["snow"]["K_domain_lb"], GRID["snow"]["K_domain_ub"] =   cryoGridTechnical.LayerIndex(GRID["snow"]["K_domain"]);

                GRID["air"]["cT_domain"][GRID["air"]["cT_domain_lb"]+1]=true;
                GRID["air"]["K_domain"][GRID["air"]["K_domain_lb"]+1]=true;
                GRID["air"]["cT_domain_lb"], GRID["air"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["cT_domain"]);
                GRID["air"]["K_domain_lb"], GRID["air"]["K_domain_ub"] =   cryoGridTechnical.LayerIndex(GRID["air"]["K_domain"]);

                #-------- rearrange SWE and T grids --------------------------
                GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_lb"]]=GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_lb"]]+GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_lb"]-1];
                GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_lb"]]=GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_lb"]]+GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_lb"]-1];
                GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_lb"]]=GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_lb"]]+GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain_lb"]-1];
                GRID["snow"]["Snow_i"][2:GRID["snow"]["cT_domain_lb"]-1]=GRID["snow"]["Snow_i"][1:GRID["snow"]["cT_domain_lb"]-2];
                GRID["snow"]["Snow_w"][2:GRID["snow"]["cT_domain_lb"]-1]=GRID["snow"]["Snow_w"][1:GRID["snow"]["cT_domain_lb"]-2];
                GRID["snow"]["Snow_a"][2:GRID["snow"]["cT_domain_lb"]-1]=GRID["snow"]["Snow_a"][1:GRID["snow"]["cT_domain_lb"]-2];

                T[GRID["snow"]["cT_domain_lb"]]=(T[GRID["snow"]["cT_domain_lb"]]+T[GRID["snow"]["cT_domain_lb"]-1])/2.0;
                T[2:GRID["snow"]["cT_domain_lb"]-1]=T[1:GRID["snow"]["cT_domain_lb"]-2];

                check_change=true;

            end


            if check_change #update grid spacings
                GRID["general"]["K_grid"][GRID["snow"]["K_domain"]] = GRID["general"]["K_grid"][GRID["snow"]["cT_domain_lb"]+1] -
                                                                      flipdim(cumsum(flipdim(GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain"]] +
                                                                      GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain"]] +
                                                                      GRID["snow"]["Snow_a"][GRID["snow"]["cT_domain"]],1)),1);

                #GRID["general"]["K_grid"][GRID["air"]["K_domain"]]=collect(GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]]+(-snowCellSize)*(GRID["air"]["cT_domain_lb"]-1):
                #snowCellSize:GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]+1]-snowCellSize);

                GRID["general"]["K_grid"][GRID["air"]["K_domain"]]=flipdim(collect(range(GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]+1]-snowCellSize, -GRID["snow"]["snowCellSize"],GRID["air"]["cT_domain_lb"])),1)
            end


            if (GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]]<=0.5.*PARA["technical"]["SWEperCell"] && sum(GRID["snow"]["cT_domain"])<2)  #remove last grid cell if snow threshold is reached

                # for water balance: add snow of last grid cell to runoff
                BALANCE["water"]["dr_snowmelt"] = BALANCE["water"]["dr_snowmelt"] - ( GRID["snow"]["Snow_i"][GRID["snow"]["cT_domain_ub"]] + GRID["snow"]["Snow_w"][GRID["snow"]["cT_domain_ub"]] ).*1000.0;

                #------ modify snow and air grid ----------------------------------
                GRID["snow"]["cT_domain"][GRID["snow"]["cT_domain_ub"]]=false;
                GRID["snow"]["K_domain"][GRID["snow"]["cT_domain_ub"]]=false;
                GRID["snow"]["cT_domain_lb"], GRID["snow"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["snow"]["cT_domain"]);
                GRID["snow"]["K_domain_lb"], GRID["snow"]["K_domain_ub"] =   cryoGridTechnical.LayerIndex(GRID["snow"]["K_domain"]);

                GRID["air"]["cT_domain"][GRID["air"]["cT_domain_lb"]+1]=true;
                GRID["air"]["K_domain"][GRID["air"]["K_domain_lb"]+1]=true;
                GRID["air"]["cT_domain_lb"], GRID["air"]["cT_domain_ub"] = cryoGridTechnical.LayerIndex(GRID["air"]["cT_domain"]);
                GRID["air"]["K_domain_lb"], GRID["air"]["K_domain_ub"] =   cryoGridTechnical.LayerIndex(GRID["air"]["K_domain"]);

                GRID["snow"]["Snow_i"][GRID["air"]["cT_domain_lb"]]=false;
                GRID["snow"]["Snow_w"][GRID["air"]["cT_domain_lb"]]=false;
                GRID["snow"]["Snow_a"][GRID["air"]["cT_domain_lb"]]=false;
                T[GRID["air"]["cT_domain_lb"]]=0.0;

                check_change=true;
            end

        end

        if check_change || check_minor_change
            # update grid spacings
            GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]] = GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]+1]-snowCellSize;
            GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]-1] = GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]+1]-2.0.*snowCellSize;
            GRID["general"]["cT_grid"] = ( GRID["general"]["K_grid"][1:end-1] + GRID["general"]["K_grid"][2:end])./2.0; #grid on which capacity and temperature information lives (midpoints of grid cells)
            GRID["general"]["cT_delta"] = (- GRID["general"]["cT_grid"][1:end-1] + GRID["general"]["cT_grid"][2:end]);
            GRID["general"]["K_delta"] = (- GRID["general"]["K_grid"][1:end-1] + GRID["general"]["K_grid"][2:end]);
        end

        return GRID, T, BALANCE
    end

end
