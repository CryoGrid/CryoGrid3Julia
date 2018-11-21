module cryoGridInfiltrationUnfrozenSoil
using MAT
using cryoGridSoil
using cryoGridSnow
using cryoGridSEB

         function conductivityUnfrozen(wc::Array{Float64,1},GRID,PARA)
                  ka = PARA["constants"]["k_a"]::Float64; #0.025;       #air [Hillel(1982)]
                  kw = PARA["constants"]["k_w"]::Float64; #0.57;        #water [Hillel(1982)]
                  ko = PARA["constants"]["k_o"]::Float64; #0.25;        #organic [Hillel(1982)]
                  km = PARA["constants"]["k_m"]::Float64; #soil.kh_bedrock;     #mineral
                  soil_cT_mineral=GRID["soil"]["cT_mineral"]::Array{Float64,1};
                  soil_cT_organic=GRID["soil"]["cT_organic"]::Array{Float64,1};

                  k_temp=zeros(size(wc));
                  @inbounds @fastmath for i=1:length(wc)
                           air=1.0-wc[i]-soil_cT_mineral[i]-soil_cT_organic[i];

                           k_temp[i] = (soil_cT_mineral[i]+soil_cT_organic[i].>1.0e-6) .*
                           (wc[i].* kw.^0.5 + soil_cT_mineral[i].* km.^0.5 + soil_cT_organic[i].* ko.^0.5 + air.* ka.^0.5).^2.0 +
                           (soil_cT_mineral[i]+soil_cT_organic[i].<=1.0e-6) .* kw; # assume pure water for cells which consist partly of air and water
                  end
                  return k_temp
         end

         function capacityUnfrozen(wc::Array{Float64,1},GRID,PARA)
                  c_w = PARA["constants"]["c_w"]::Float64; # 4.2*10^6; #[J/m�K]
                  c_o = PARA["constants"]["c_o"]::Float64; # 2.5*10^6; #[J/m�K]
                  c_m = PARA["constants"]["c_m"]::Float64; # 2*10^6; #[J/m�K]
                  soil_cT_mineral=GRID["soil"]["cT_mineral"]::Array{Float64,1};
                  soil_cT_organic=GRID["soil"]["cT_organic"]::Array{Float64,1};

                  c_temp=zeros(size(wc));
                  @inbounds @fastmath for i=1:length(wc)
                           c_temp[i] = (soil_cT_mineral[i]+soil_cT_organic[i].>1.0e-6) .* (soil_cT_mineral[i].*c_m + soil_cT_organic[i].*c_o + wc[i].*c_w) +
                           (soil_cT_mineral[i]+soil_cT_organic[i].<=1.0e-6) .* c_w; # assume pure water for cells which consist partly of air and water
                  end
                  return c_temp
         end

         function getThermalPropertiesInfiltration(T::Array{Float64,1}, wc::Array{Float64,1}, c_cTgrid::Array{Float64,1}, k_cTgrid::Array{Float64,1}, k_Kgrid::Array{Float64,1}, lwc_cTgrid::Array{Float64,1}, GRID, PARA)
                  c_a=PARA["constants"]["c_a"]::Float64;
                  k_a=PARA["constants"]["k_a"]::Float64;

                  c_w = PARA["constants"]["c_w"]::Float64; # 4.2*10^6; #[J/m�K]
                  c_o = PARA["constants"]["c_o"]::Float64; # 2.5*10^6; #[J/m�K]
                  c_m = PARA["constants"]["c_m"]::Float64; # 2*10^6; #[J/m�K]
                  ka = PARA["constants"]["k_a"]::Float64; #0.025;       #air [Hillel(1982)]
                  kw = PARA["constants"]["k_w"]::Float64; #0.57;        #water [Hillel(1982)]
                  ko = PARA["constants"]["k_o"]::Float64; #0.25;        #organic [Hillel(1982)]
                  km = PARA["constants"]["k_m"]::Float64; #soil.kh_bedrock;     #mineral
                  soil_cT_mineral=GRID["soil"]["cT_mineral"]::Array{Float64,1};
                  soil_cT_organic=GRID["soil"]["cT_organic"]::Array{Float64,1};

                  air_cT_domain=GRID["air"]["cT_domain"]::Array{Bool,1};
                  air_K_domain_lb=GRID["air"]["K_domain_lb"]::Int64;
                  soil_cT_domain=GRID["soil"]["cT_domain"]::Array{Bool,1};
                  soil_cT_domain_ub=GRID["soil"]["cT_domain_ub"]::Int64
                  snow_cT_domain=GRID["snow"]["cT_domain"]::Array{Bool,1};
                  snow_Snow_i=GRID["snow"]["Snow_i"]::Array{Float64,1};
                  snow_Snow_w=GRID["snow"]["Snow_w"]::Array{Float64,1};
                  snow_Snow_a=GRID["snow"]["Snow_a"]::Array{Float64,1};
                  general_K_delta=GRID["general"]["K_delta"]::Array{Float64,1};
                  general_cT_delta=GRID["general"]["cT_delta"]::Array{Float64,1};

                  cT_frozen = GRID["soil"]["cT_frozen"]::Array{Float64,1};
                  cT_thawed = GRID["soil"]["cT_thawed"]::Array{Float64,1};
                  capacity = GRID["soil"]["capacity"]::Array{Float64,2};
                  K_frozen = GRID["soil"]["K_frozen"]::Array{Float64,1};
                  K_thawed = GRID["soil"]["K_thawed"]::Array{Float64,1};
                  arraySizeT = PARA["technical"]["arraySizeT"]::Int64
                  conductivity = GRID["soil"]["conductivity"]::Array{Float64,2};
                  liquidWaterContent = GRID["soil"]["liquidWaterContent"]::Array{Float64,2}; #added by JAN for liquid water content

                  c_temp=zeros(size(c_cTgrid));
                  k_temp=zeros(size(k_cTgrid));
                  lwc_temp=zeros(size(lwc_cTgrid));
                  k_eff=zeros(size(k_cTgrid));

                  @inbounds @fastmath for i=1:length(T)
                           #c_old = c_temp;
                           #------- unused grid cells --------------------------------------------
                           if air_cT_domain[i]==true
                                    c_temp[i] = c_a;
                                    k_temp[i] = k_a;
                                    lwc_temp[i] = 0.0;
                           end

                           #------- soil domain --------------------------------------------------
                           if soil_cT_domain[i]==true
                                    #c_temp_, k_temp_, lwc_temp_ = cryoGridSoil.readThermalParameters(T[i], GRID, PARA);
                                    #c_temp[i]=c_temp_[1];
                                    #k_temp[i]=k_temp_[1];
                                    #lwc_temp[i]=lwc_temp_[1];
                                    j=i-soil_cT_domain_ub+1
                                    #a=round(Int64,(T[i]-cT_frozen[j])./(cT_thawed[j]-cT_frozen[j])*(arraySizeT-2)+1); #T and c information live on same grid
                                    #if a<1
                                    #         a=1
                                    #elseif a>arraySizeT
                                    #         a=arraySizeT
                                    #end

                                    #c_temp_=capacity[j,a];
                                    #lwc_temp_=liquidWaterContent[j,a]; #added by JAN for liquid water content
                                    #k_temp_=conductivity[j,a];

                                    #adjust for the unfrozen part of the domain
                                    #JAN: this changes with infiltration scheme: frozen  (T<=0) remains unchanged,
                                    #thawed (T>0) is calculated differently in dependence of wc
                                    #c_temp_unfrozen = (soil_cT_mineral[j]+soil_cT_organic[j].>1.0e-6).*
                                    #                  (soil_cT_mineral[j].*c_m + soil_cT_organic[j].*c_o + wc[j].*c_w)+
                                    #                  (soil_cT_mineral[j]+soil_cT_organic[j].<=1.0e-6).*c_w; # assume pure water for cells which consist partly of air and water


                                    #c_temp[i] = (T[i].<=0.0).*c_temp_ + (T[i].>0.0).* c_temp_unfrozen;


                                    #air=1.0-wc[j]-soil_cT_mineral[j]-soil_cT_organic[j];
                                    #k_temp_unfrozen = (soil_cT_mineral[j]+soil_cT_organic[j].>1.0e-6) .*
                                    #(wc[j].* kw.^0.5 + soil_cT_mineral[j].* km.^0.5 + soil_cT_organic[j].* ko.^0.5 + air.* ka.^0.5).^2.0 +
                                    #(soil_cT_mineral[j]+soil_cT_organic[j].<=1.0e-6) .* kw; # assume pure water for cells which consist partly of air and water

                                    #k_temp[i] = (T[i].<=0.0).*k_temp_ + (T[i].>0.0).* k_temp_unfrozen;
                                    #lwc_temp[i] = (T[i].<=0.0).*lwc_temp_ + (T[i].>0.0).* wc[i-soil_cT_domain_ub+1];
                                    if T[i]<=0.0
                                             a=round(Int64,(T[i]-cT_frozen[j])./(cT_thawed[j]-cT_frozen[j])*(arraySizeT-2)+1); #T and c information live on same grid
                                             if a<1
                                                      a=1
                                             elseif a>arraySizeT
                                                      a=arraySizeT
                                             end

                                             c_temp_=capacity[j,a];
                                             lwc_temp_=liquidWaterContent[j,a]; #added by JAN for liquid water content
                                             k_temp_=conductivity[j,a];

                                             c_temp[i] = c_temp_::Float64;
                                             k_temp[i] = k_temp_::Float64;
                                             lwc_temp[i] = lwc_temp_::Float64;
                                    else
                                             c_temp_unfrozen = (soil_cT_mineral[j]+soil_cT_organic[j].>1.0e-6).*
                                             (soil_cT_mineral[j].*c_m + soil_cT_organic[j].*c_o + wc[j].*c_w)+
                                             (soil_cT_mineral[j]+soil_cT_organic[j].<=1.0e-6).*c_w; # assume pure water for cells which consist partly of air and water

                                             air=1.0-wc[j]-soil_cT_mineral[j]-soil_cT_organic[j];
                                             k_temp_unfrozen = (soil_cT_mineral[j]+soil_cT_organic[j].>1.0e-6) .*
                                             (wc[j].* kw.^0.5 + soil_cT_mineral[j].* km.^0.5 + soil_cT_organic[j].* ko.^0.5 + air.* ka.^0.5).^2.0 +
                                             (soil_cT_mineral[j]+soil_cT_organic[j].<=1.0e-6) .* kw; # assume pure water for cells which consist partly of air and water

                                             c_temp[i]=c_temp_unfrozen::Float64;
                                             k_temp[i] = k_temp_unfrozen::Float64;
                                             lwc_temp[i] = wc[i-soil_cT_domain_ub+1]::Float64;
                                    end
                           end
                           #------- snow domain --------------------------------------------------
                           if snow_cT_domain[i]==true
                                    c_temp[i] = cryoGridSnow.cap_snow(snow_Snow_i[i],snow_Snow_w[i],snow_Snow_a[i],PARA)::Float64;
                                    k_temp[i] = cryoGridSnow.cond_snow(snow_Snow_i[i],snow_Snow_w[i],snow_Snow_a[i])::Float64;
                                    lwc_temp[i] = snow_Snow_w[i]./general_K_delta[i]::Float64;
                           end
                           #------- interpolate conductivity to K-grid ---------------------------
                           if (i>=2 && i<=length(k_eff)-1)
                                    k_eff_ = general_K_delta[i-1]./(2.0.*general_cT_delta[i]) .* (1.0./k_temp[i-1]).^2.0 +
                                             general_K_delta[i]./(2.0.*general_cT_delta[i]) .* (1.0./k_temp[i]).^2.0;

                                    k_eff[i] = k_eff_.^(-0.5)::Float64;
                           end
                  end
                  k_eff[1]     = k_temp[1]::Float64;
                  k_eff[end]   = k_temp[end]::Float64;

                  #------ correct upper most value below air-domain ---------------------
                  k_eff[air_K_domain_lb+1] = k_temp[air_K_domain_lb+1]::Float64;

                  return c_temp, k_temp, k_eff, lwc_temp
         end

         function surfaceEnergyBalanceInfiltration(T::Array{Float64,1}, wc::Array{Float64,1}, FORCING, GRID, PARA, SEB)
                  Lstar=mean(SEB["L_star"])::Float64;

                  sigma=PARA["constants"]["sigma"]::Float64; #5.67e-8; #Stefan-Boltzmann const.
                  L=PARA["constants"]["L_lg"].*PARA["constants"]["rho_w"]::Float64;  #2.8*10^6.*1000;   #check this, this seems to be for ice? JAN: yes, should be smaller!
                  z=PARA["technical"]["z"];

                  surf_z0=PARA["surf"]["z0"]::Float64;
                  surf_albedo=PARA["surf"]["albedo"]::Float64;
                  surf_epsilon=PARA["surf"]["epsilon"]::Float64;
                  surf_rs=PARA["surf"]["rs"]::Float64;

                  i_wind=FORCING["i"]["wind"]::Float64;
                  i_Tair=FORCING["i"]["Tair"]::Float64;
                  i_p=FORCING["i"]["p"]::Float64;
                  i_q=FORCING["i"]["q"]::Float64;
                  i_Sin=FORCING["i"]["Sin"]::Float64;
                  i_Lin=FORCING["i"]["Lin"]::Float64;

                  air_cT_domain_lb=GRID["air"]["cT_domain_lb"]::Int64;
                  general_cT_grid=GRID["general"]["cT_grid"]::Array{Float64,1};
                  snow_cT_domain_ub=GRID["snow"]["cT_domain_ub"];
                  snow_cT_domain_lb=GRID["snow"]["cT_domain_lb"];
                  general_K_grid=GRID["general"]["K_grid"]::Array{Float64,1};
                  soil_cT_domain_ub=GRID["soil"]["cT_domain_ub"]::Int64;
                  general_K_delta=GRID["general"]["K_delta"]::Array{Float64,1};


                  Qh=real(cryoGridSEB.Q_h(i_wind, z, surf_z0, i_Tair, T[air_cT_domain_lb+1], Lstar, i_p, i_q, PARA));

                  dwc_dt=zeros(size(wc));

                  #______here SW radiation is calculated_____________________________________
                  dE_dt=zeros(size(general_cT_grid));
                  Qsolar=zeros(size(general_cT_grid));

                  dE_dt[air_cT_domain_lb+1]=(1.0-surf_albedo)*i_Sin;
                  #------ snow surface (solid state green house effect) ---------------------
                  if !isempty(snow_cT_domain_ub)
                           beta=PARA["snow"]["extinction"]::Float64;
                           Qsolar[snow_cT_domain_ub:snow_cT_domain_lb+1] = dE_dt[snow_cT_domain_ub].*
                                    exp.(-beta.*(general_K_grid[snow_cT_domain_ub:snow_cT_domain_lb+1]-general_K_grid[snow_cT_domain_ub]));

                           dE_dt[snow_cT_domain_ub:snow_cT_domain_lb] = -Qsolar[snow_cT_domain_ub+1:snow_cT_domain_lb+1]+
                                    Qsolar[snow_cT_domain_ub:snow_cT_domain_lb];
                           #put the rest to cell below snow
                           dE_dt[snow_cT_domain_lb+1] = Qsolar[snow_cT_domain_lb+1];
                  end

                  # JAN : here a modification for water bodies would be needed (extinction),
                  # but this would probably make no difference for summer due to mixing,
                  # maybe different in winter/spring

                  #__________________________________________________________________________
                  Sout = surf_albedo*i_Sin;
                  Lout = surf_epsilon*sigma*(T[air_cT_domain_lb+1]+273.15).^4.0 + (1.0-surf_epsilon).*i_Lin;
                  Qnet = i_Sin-Sout + i_Lin - Lout ;

                  #calculate ET
                  if PARA["modules"]["infiltration"]

                           #snow cover or uppermost grid cell frozen --> no ET ; JAN: this includes the case of a frozen water body
                           if !isempty(snow_cT_domain_ub) || T[soil_cT_domain_ub]<=0.0
                                    Qe=real(cryoGridSEB.Q_eq(i_wind, z, surf_z0, i_q, i_Tair, T[air_cT_domain_lb+1], Lstar, surf_rs, i_p, PARA));
                                    # unfrozen water body at surface
                           elseif GRID["lake"]["unfrozenWaterSurface"]
                                    Qe=real(cryoGridSEB.Q_eq(i_wind, z, surf_z0, i_q, i_Tair, T[air_cT_domain_lb+1], Lstar, surf_rs, i_p, PARA));
                                    dwc_dt[1]=-Qe./L; #in m water per sec, this can be evaporation or condensation

                                    # JAN: this is the "default" case of an unfrozen soil surface
                           else
                                    Qe_pot=real(cryoGridSEB.Q_eq(i_wind, z, surf_z0, i_q, i_Tair, T[air_cT_domain_lb+1], Lstar, 0.0, i_p, PARA));  #potential ET
                                    if Qe_pot>0.0
                                             fraction_T = getET_fraction(T[soil_cT_domain_ub:soil_cT_domain_ub+GRID["soil"]["T_lb"]-1], wc[1:GRID["soil"]["T_lb"]], PARA["soil"]["fieldCapacity"], PARA["soil"]["wiltingPoint"]);
                                             fraction_E = getET_fraction(T[soil_cT_domain_ub:soil_cT_domain_ub+GRID["soil"]["E_lb"]-1], wc[1:GRID["soil"]["E_lb"]], PARA["soil"]["fieldCapacity"], PARA["soil"]["residualWC"]);
                                             fraction_ET = fraction_T.*PARA["soil"]["ratioET"];
                                             fraction_ET[1:GRID["soil"]["E_lb"]] = fraction_ET[1:GRID["soil"]["E_lb"]] + fraction_E.*(1.0-PARA["soil"]["ratioET"]);

                                             Qe=sum(fraction_ET.*general_K_delta[soil_cT_domain_ub:soil_cT_domain_ub+GRID["soil"]["T_lb"]-1])./sum(general_K_delta[soil_cT_domain_ub:soil_cT_domain_ub+GRID["soil"]["T_lb"]-1]).*Qe_pot;
                                             fraction_ET=fraction_ET.*general_K_delta[soil_cT_domain_ub:soil_cT_domain_ub+GRID["soil"]["T_lb"]-1]./sum(fraction_ET.*general_K_delta[soil_cT_domain_ub:soil_cT_domain_ub+GRID["soil"]["T_lb"]-1]);
                                             # sum(fraction_ET) is always 1
                                             dwc_dt[1:GRID["soil"]["T_lb"]]=-Qe./L.*fraction_ET;    #in m water per sec
                                    else  #condensation
                                             Qe=Qe_pot;
                                             dwc_dt[1]=-Qe./L; #in m water per sec, put everything in uppermost grid cell
                                    end
                           end
                  else # this is identical to case with snow cover or frozen ground
                           Qe=real(cryoGridSEB.Q_eq(i_wind, z, surf_z0, i_q, i_Tair, T[air_cT_domain_lb+1], Lstar, surf_rs, i_p, PARA));
                  end
                  #ground heat flux
                  Qg = Qnet-Qh-Qe;

                  #surface heat flux (into upper cell, ground heat flux regards also other
                  #grid cells, should be identical if no snow cover and no evapotranspiration
                  #occur
                  dE_dt[air_cT_domain_lb+1] = dE_dt[air_cT_domain_lb+1]+
                  surf_epsilon.*i_Lin-
                  surf_epsilon.*sigma.*(T[air_cT_domain_lb+1]+273.15).^4.0- Qh - Qe;  # Qe positive: cooling of soil => evaporation/subl. => loss of SWE

                  # fluxes are in [ W / m^2 ]
                  SEB["Qsurf"] = dE_dt[air_cT_domain_lb+1];

                  # if abs( SEB["Qsurf"]-Qg ) > 1e-6
                  #     warning ( ' Qsurf != Qg ' );
                  # end

                  SEB["dE_dt_SEB"] = dE_dt;
                  SEB["Qnet"] = Qnet;
                  SEB["Qh"] = Qh;
                  SEB["Qe"] = Qe;
                  SEB["Qg"] = Qg;
                  SEB["Sout"] = Sout;
                  SEB["Lout"] = Lout;
                  return SEB, dwc_dt
         end

         function getET_fraction(T::Array{Float64,1}, wc::Array{Float64,1}, start_reduction, zero_reached)
                  #fraction=double(T>0).*min(1, max(0, (wc-zero_reached)./(start_reduction-zero_reached)));
                  fraction=(T.>0.0).*((wc.>=start_reduction) + (wc.<start_reduction).*0.25.*(1.0-cos.(pi.*wc./start_reduction)).^2.0);
                  return fraction
         end

         function CryoGridInfiltration(T, wc, dwc_dt, timestep, GRID, PARA, FORCING, BALANCE)

                  # possible meltwater contribution from xice
                  if !PARA["modules"]["xice"]
                           meltwaterGroundIce = 0.0;
                  else
                           meltwaterGroundIce = GRID["lake"]["residualWater"];
                           GRID["lake"]["residualWater"] = 0.0;
                  end

                  # external flux
                  external_flux_rate = PARA["soil"]["externalWaterFlux"];             # in m/day
                  BALANCE["water"]["dr_subsurface"] = BALANCE["water"]["dr_subsurface"] + external_flux_rate.*timestep.*1000.0;    #in mm

                  if isempty(GRID["snow"]["cT_domain_ub"]) && T[GRID["soil"]["cT_domain_ub"]]>0.0   #no snow cover and uppermost grid cell unfrozen

                           ### step 1: infiltrate rain and meltwater and external flux through bucket scheme
                           # changes due to evapotranspiration and condensation
                           dwc_dt=dwc_dt.*timestep.*24.0.*3600.0;   #now in m water per grid cell
                           BALANCE["water"]["de"] = BALANCE["water"]["de"] + sum(dwc_dt)*1000.0; # in mm accumulated over soil column

                           # changes due to rainfall
                           dwc_dt[1]=dwc_dt[1]+FORCING["i"]["rainfall"]./1000.0.*timestep;

                           # changes due to meltwater from excess ice
                           dwc_dt[1]=dwc_dt[1]+meltwaterGroundIce;

                           # routing of water
                           wc, surface_runoff = bucketScheme(T, wc, dwc_dt, GRID, PARA, external_flux_rate.*timestep);

                           # consistency check
                           if sum( wc.<0.0 )!=0.0
                                    print("warning negative water content occured");
                                    #here one could correct the water balance
                           end

                           # remove water above water table in case of ponding, e.g. through rain (independent of xice module)
                           if GRID["soil"]["cT_mineral"][1]+GRID["soil"]["cT_organic"][1]<1e-6 &&
                                    GRID["general"]["K_grid"][GRID["soil"]["cT_domain_ub"]]<PARA["soil"]["waterTable"]

                                    cellSize = GRID["general"]["K_delta"][GRID["soil"]["cT_domain_ub"]];
                                    actualWater = wc[1]*cellSize;
                                    h = GRID["general"]["K_grid"][GRID["soil"]["K_domain_ub"]+1]-PARA["soil"]["waterTable"];
                                    if h<0.0
                                             print("warning h<0. too much water above water table!")
                                    end

                                    if actualWater>h
                                             print("infiltration - removing excess water from upper cell");
                                             wc[1]=h./cellSize;
                                             surface_runoff = surface_runoff + actualWater-h;
                                    end

                           end

                           ### step 2: update GRID including reomval of excess water above water table and ponding below water table
                           wc, GRID, surface_runoff = updateGRID_infiltration(wc, GRID, PARA, surface_runoff);

                           # store remaining surface runoff
                           BALANCE["water"]["dr_surface"] = BALANCE["water"]["dr_surface"] - surface_runoff*1000.0; # in [mm]
                  end

                  # step 3: LUT update
                  #         JAN:recalculate lookup tables when water content of freezing grid cells
                  #         has changed (infiltrated cells can freeze --> LUT is updated)
                  if sum((wc.!=GRID["soil"]["cT_water"]) .& (T[GRID["soil"]["cT_domain"]].<=0.0))>0
                           print("infiltration - reinitializing LUT - freezing of infiltrated cell(s)");
                           GRID["soil"]["cT_water"] = wc;
                           GRID = cryoGridSoil.initializeSoilThermalProperties(GRID, PARA);
                  end
                  return wc, GRID, BALANCE
         end

         function bucketScheme(T::Array{Float64,1}, wc::Array{Float64,1}, dwc_dt::Array{Float64,1}, GRID, PARA, external_flux::Float64)

                  T_soil=T[GRID["soil"]["cT_domain"]]::Array{Float64,1};
                  K_delta=GRID["general"]["K_delta"][GRID["soil"]["cT_domain"]]::Array{Float64,1};  #in m
                  porosity=1.0-GRID["soil"]["cT_mineral"]::Array{Float64,1}-GRID["soil"]["cT_organic"]::Array{Float64,1};

                  T_soil=T[GRID["soil"]["cT_domain"]]::Array{Float64,1};  #in percent  JAN: why not use GRID["soil"]["cT_natPor"] here?

                  # A=sum(K_delta(1:30).*wc(1:30)+dwc_dt(1:30))
                  max_water=0.0::Float64
                  actual_water=0.0::Float64
                  soil_fieldCapacity=PARA["soil"]["fieldCapacity"]::Float64;

                  i=1::Int64;
                  i_max=70::Int64;  # maximum infiltration depth, must be defined somehow before
                  @inbounds @fastmath while T_soil[i]>0.0 && i<=i_max
                           max_water=K_delta[i].*soil_fieldCapacity::Float64;  #maximum amount of water (in m) that a grid cell can hold
                           actual_water=wc[i].*K_delta[i]+dwc_dt[i];       # should be dwc (already multiplied with timestep)
                           dwc_dt[i+1]=dwc_dt[i+1] + max(0.0, actual_water-max_water);  #when excess water, move it to next grid cell
                           wc[i]=min(max_water, actual_water)./K_delta[i];

                           i=i+1;
                  end

                  excess_water=dwc_dt[i]+external_flux::Float64; #add external flux

                  # B=sum(K_delta(1:30).*wc(1:30))+excess_water
                  i=i-1;

                  @inbounds @fastmath while i>=1 && excess_water>0.0
                           max_water=K_delta[i].*porosity[i]::Float64;
                           actual_water=wc[i].*K_delta[i]+excess_water;
                           wc[i]=min(actual_water, max_water)./K_delta[i];
                           excess_water=max(0.0, actual_water-wc[i].*K_delta[i]);

                           i=i-1;
                  end

                  surface_runoff=excess_water::Float64;

                  # C=sum(K_delta(1:30).*wc(1:30))+surface_runoff
                  return wc, surface_runoff
         end

         function updateGRID_infiltration(wc, GRID, PARA, surface_runoff)

                  ### step 2: GRID update
                  ### TODO: add a function updateGRID_infiltration

                  soilGRIDsizeOld = sum(GRID["soil"]["cT_domain"]);

                  ### step 2a) remove cells filled with air (e.g. due to evaporation
                  ### of uppermost grid cell )
                  @inbounds @fastmath while GRID["soil"]["cT_mineral"][1]+GRID["soil"]["cT_organic"][1]+wc[1]<=0.0
                           print("infiltration - update GRID - removing air cell")

                           # adjust air and soil domains and boundaries
                           GRID["air"]["cT_domain"][GRID["soil"]["cT_domain_ub"]]=true;
                           GRID["air"]["K_domain"][GRID["soil"]["K_domain_ub"]]=true;
                           GRID["air"]["cT_domain_lb"]=GRID["air"]["cT_domain_lb"]+1;
                           GRID["air"]["K_domain_lb"]=GRID["air"]["K_domain_lb"]+1;
                           GRID["soil"]["cT_domain"][GRID["soil"]["cT_domain_ub"]]=false;
                           GRID["soil"]["K_domain"][GRID["soil"]["K_domain_ub"]]=false;
                           GRID["soil"]["cT_domain_ub"]=GRID["soil"]["cT_domain_ub"]+1;
                           GRID["soil"]["K_domain_ub"]=GRID["soil"]["K_domain_ub"]+1;
                           deleteat!(GRID["soil"]["soilGrid"],1);

                           deleteat!(wc,1);

                           deleteat!(GRID["soil"]["cT_organic"],1);
                           deleteat!(GRID["soil"]["cT_natPor"],1);
                           deleteat!(GRID["soil"]["cT_mineral"],1);
                           deleteat!(GRID["soil"]["cT_soilType"],1);
                           # K fields are not used currently
                           #                 GRID["soil"]["K_water"](1)=[];
                           #                 GRID["soil"]["K_organic"](1)=[];
                           #                 GRID["soil"]["K_mineral"](1)=[];
                           #                 GRID["soil"]["K_soilType"](1)=[];
                           #deleteat!(GRID["soil"]["excessGroundIce"],1);

                  end

                  ### step 2b) ponding of surface runoff below water table
                  @inbounds @fastmath while surface_runoff>1e-6 &&                                # not >0 as sometimes numerical errors occur during calculation of surface_runoff
                           GRID["general"]["K_grid"][GRID["soil"]["K_domain_ub"]]>PARA["soil"]["waterTable"] #&& ...
                           #wc(1)>=1   # this prevents a bug for very small
                           #surface_runoff when upper cell not filled // but this
                           #does not allow ponding on top of actual soil
                           print("infiltration - update GRID - ponding of water below water table")

                           h = GRID["general"]["K_grid"][GRID["soil"]["K_domain_ub"]]-PARA["soil"]["waterTable"];    # this is guruanteed to be >0


                           # create new water cell / change GRID domains
                           GRID["soil"]["cT_domain"][GRID["air"]["cT_domain_lb"]]=true;
                           GRID["soil"]["K_domain"][GRID["air"]["K_domain_lb"]]=true;
                           GRID["soil"]["cT_domain_ub"]=GRID["soil"]["cT_domain_ub"]-1;
                           GRID["soil"]["K_domain_ub"]=GRID["soil"]["K_domain_ub"]-1;
                           GRID["air"]["cT_domain"][GRID["air"]["cT_domain_lb"]]=false;
                           GRID["air"]["K_domain"][GRID["air"]["K_domain_lb"]]=false;
                           GRID["air"]["cT_domain_lb"]=GRID["air"]["cT_domain_lb"]-1;
                           GRID["air"]["K_domain_lb"]=GRID["air"]["K_domain_lb"]-1;

                           # fill new water cell
                           #cellSize = GRID["general"]["K_delta"](GRID["soil"]["cT_domain_ub"]);
                           cellSize = PARA["technical"]["waterCellSize"];
                           waterAdded = min(surface_runoff, cellSize, h); # add water until water table is reached or surface_runoff "empty"
                           wc = collect([waterAdded./cellSize ; wc]);
                           surface_runoff = surface_runoff - waterAdded;

                           # update remaining soil fields with exception of cT_water
                           GRID["soil"]["cT_organic"]  = collect([ 0.0 ; GRID["soil"]["cT_organic"] ]);
                           GRID["soil"]["cT_natPor"]   = collect([ GRID["soil"]["cT_natPor"][1]; GRID["soil"]["cT_natPor"] ]);    # take natPor of cell below
                           GRID["soil"]["cT_mineral"]  = collect([ 0.0 ; GRID["soil"]["cT_mineral"] ]);
                           GRID["soil"]["cT_soilType"] = collect([ 1.0; GRID["soil"]["cT_soilType"]]);                        # assume sand as soil type for water cell
                           # K fields are not used currently
                           #GRID["soil"]["K_water"] = [ wc(1); GRID["soil"]["K_water"] ];
                           #GRID["soil"]["K_organic"] = [ 0 ; GRID["soil"]["K_organic"] ];
                           #GRID["soil"]["K_mineral"] = [ 0 ; GRID["soil"]["K_mineral"] ];
                           #GRID["soil"]["K_soilType"] = [ GRID["soil"]["K_soilType"](1); GRID["soil"]["K_soilType"]];
                           #GRID["soil"]["excessGroundIce"] = collect([ 0.0 ; GRID["soil"]["excessGroundIce"] ]);

                           # update GRID spacings
                           GRID["general"]["K_grid"][GRID["soil"]["cT_domain_ub"]] = GRID["general"]["K_grid"][GRID["soil"]["cT_domain_ub"]+1]-cellSize;
                           GRID["general"]["K_grid"][GRID["air"]["cT_domain"]] = collect(GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]]+(-GRID["snow"]["snowCellSize"])*(GRID["air"]["cT_domain_lb"]-1.0):
                                                                                          GRID["snow"]["snowCellSize"]:
                                                                                          GRID["general"]["K_grid"][GRID["air"]["cT_domain_lb"]]);
                           GRID["general"]["cT_grid"] = (GRID["general"]["K_grid"][1:end-1]+ GRID["general"]["K_grid"][2:end])/2.0; #grid on which capacity and temperature information lives (midpoints of grid cells)
                           GRID["general"]["cT_delta"] = (- GRID["general"]["cT_grid"][1:end-1] + GRID["general"]["cT_grid"][2:end]);
                           GRID["general"]["K_delta"] = (- GRID["general"]["K_grid"][1:end-1] + GRID["general"]["K_grid"][2:end]);
                           GRID["soil"]["soilGrid"] = collect([ GRID["general"]["K_grid"][GRID["soil"]["cT_domain_ub"]] ; GRID["soil"]["soilGrid"] ]);
                  end

                  ### step 2c)  check if soil/air domains changed --> LUT update
                  soilGRIDsizeNew = sum(GRID["soil"]["cT_domain"]);
                  if soilGRIDsizeOld!=soilGRIDsizeNew
                           print("infiltration - reinitializing LUT - soil/air domains changed");
                           GRID["soil"]["cT_water"] = wc;
                           GRID = cryoGridSoil.initializeSoilThermalProperties(GRID, PARA);
                  end

                  return wc, GRID, surface_runoff
         end

end
