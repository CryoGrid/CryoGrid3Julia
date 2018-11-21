module cryoGridSoil
    using matlab
    function createStratigraphy(PARA, GRID)

        i=1;
        soilParam=Array{Float64}(size(PARA["soil"]["layer_properties"],1)*2,size(PARA["soil"]["layer_properties"],2))
        soilParam[1,:]=[-100 PARA["soil"]["layer_properties"][i, 2:6]' ];
        i=2;
        ii=2;
        while i<=size(PARA["soil"]["layer_properties"],1)
            soilParam[ii,:]=[PARA["soil"]["layer_properties"][i, 1]-0.01 PARA["soil"]["layer_properties"][i-1, 2:6]'];
            soilParam[ii+1,:]=[PARA["soil"]["layer_properties"][i, 1]+0.01 PARA["soil"]["layer_properties"][i, 2:6]'];
            i=i+1;
            ii=ii+2;
        end
        soilParam[end,:]=[1e6 PARA["soil"]["layer_properties"][i-1, 2:6]'];


        #interpolate to grid
        GRID["soil"]["cT_water"] = matlab.interp1(soilParam[:,1],soilParam[:,2],GRID["general"]["cT_grid"][GRID["soil"]["cT_domain"]],"linear");
        GRID["soil"]["cT_mineral"] = matlab.interp1(soilParam[:,1],soilParam[:,3],GRID["general"]["cT_grid"][GRID["soil"]["cT_domain"]],"linear");
        GRID["soil"]["cT_organic"] = matlab.interp1(soilParam[:,1],soilParam[:,4],GRID["general"]["cT_grid"][GRID["soil"]["cT_domain"]],"linear");
        GRID["soil"]["cT_soilType"] = matlab.interp1(soilParam[:,1],soilParam[:,5],GRID["general"]["cT_grid"][GRID["soil"]["cT_domain"]],"nearest");
        GRID["soil"]["cT_natPor"] = matlab.interp1(soilParam[:,1],soilParam[:,6],GRID["general"]["cT_grid"][GRID["soil"]["cT_domain"]],"linear");

        GRID["soil"]["K_water"] = matlab.interp1(soilParam[:,1],soilParam[:,2],GRID["general"]["K_grid"][GRID["soil"]["K_domain"]],"linear");
        GRID["soil"]["K_mineral"] = matlab.interp1(soilParam[:,1],soilParam[:,3],GRID["general"]["K_grid"][GRID["soil"]["K_domain"]],"linear");
        GRID["soil"]["K_organic"] = matlab.interp1(soilParam[:,1],soilParam[:,4],GRID["general"]["K_grid"][GRID["soil"]["K_domain"]],"linear");
        GRID["soil"]["K_soilType"] = matlab.interp1(soilParam[:,1],soilParam[:,5],GRID["general"]["K_grid"][GRID["soil"]["K_domain"]],"nearest");
        return GRID
    end

    #---------------------------------------------------
    # function initialize
    #creates matrices for heat capacity and conductivity
    function initializeSoilThermalProperties(GRID, PARA)

        cT_water = GRID["soil"]["cT_water"]::Array{Float64,1};
        cT_mineral = GRID["soil"]["cT_mineral"]::Array{Float64,1};
        cT_organic = GRID["soil"]["cT_organic"]::Array{Float64,1};
        cT_soilType = GRID["soil"]["cT_soilType"]::Array{Float64,1};
        K_water = GRID["soil"]["K_water"]::Array{Float64,1};
        K_mineral = GRID["soil"]["K_mineral"]::Array{Float64,1};
        K_organic = GRID["soil"]["K_organic"]::Array{Float64,1};
        K_soilType = GRID["soil"]["K_soilType"]::Array{Float64,1};
        arraySize = PARA["technical"]["arraySizeT"]::Int64;
        cT_grid = GRID["general"]["cT_grid"][GRID["soil"]["cT_domain"]]::Array{Float64,1};
        kh_bedrock = PARA["soil"]["kh_bedrock"]::Float64;

        c_w = PARA["constants"]["c_w"]::Float64; #4.2*10^6; #[J/m�K]
        c_o = PARA["constants"]["c_o"]::Float64; #2.5*10^6; #[J/m�K]
        c_m = PARA["constants"]["c_m"]::Float64; #2*10^6; #[J/m�K]
        c_a = PARA["constants"]["c_a"]::Float64; #0.00125*10^6;#[J/m�K]
        c_i = PARA["constants"]["c_i"]::Float64; #1.9*10^6;#[J/m�K]

        #density of water
        rho_w = PARA["constants"]["rho_w"]::Float64; #1000; #[kg/m�]
        #rho_i=900;
        #latent heat of freezing
        L_si = PARA["constants"]["L_sl"]::Float64; #334000; # [J/kg]

        # JAN: modification to assume pure water for mixed air/water cells
        cT_water[cT_mineral+cT_organic.<=1e-6]=1.0;

        #------- capacity part ----------------------------------------------------
        waterMin=0.0;
        water=cT_water;
        mineral=cT_mineral;
        organic=cT_organic;
        a=cT_soilType;
        a_length=length(a);
        deltaT=0.001*ones(a_length);

        cT_thawed=zeros(a_length);
        #cT_frozen=-15*ones(size(a,1),1);

        # JAN: Why the additional loop from -30 to -1??? --> determine cT_frozen
        #preallocate variable
        dchdT=ones(a_length,length(-30.0:0.01:-1.0));
        c_h2o_check=ones(a_length,length(-30.0:0.01:-1.0));
        ch_check=mineral*c_m+organic*c_o+waterMin.*c_w+(water-waterMin)*c_i::Float64;
        j=1;
        @inbounds @fastmath for i= -30.0:0.01:-1.0
            dchdT[:,j]=L_si*rho_w*(cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, i+deltaT/2.0, PARA)-cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, i-deltaT/2.0, PARA))/deltaT[1]
            c_h2o_check[:,j] = 1.0.*(dchdT[:,j] .< 0.05.*ch_check);
            j=j+1;

        end


        #preallocate variables
        cT_frozen=Array{Float64}(a_length)
        cT_frozen=-30.0+(squeeze(sum(c_h2o_check',1)',2)-1.0).*0.01;


        c_h2o=ones(a_length,arraySize);
        water_c=ones(a_length,arraySize);
        ch=ones(a_length,arraySize);

        water_c[:,1] = cryoGridSoil.freezeC(water,1-mineral-organic, a, cT_frozen, PARA);
        ch[:,1]      = mineral * c_m + organic * c_o +  water_c[:,1] * (c_w-c_i) + water * c_i;
        #here the derivative of the freeze curve dwc / dt is computed
        c_h2o[:,1]   = L_si*rho_w* (cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, cT_frozen+deltaT/2.0, PARA)-cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, cT_frozen-deltaT/2.0, PARA))/deltaT[1];
        @inbounds @fastmath for i=1:arraySize-2
            T_tab=(cT_frozen-cT_thawed)*(arraySize-2-i)/(arraySize-2)
            water_c[:,i+1] = cT_thawed+cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, T_tab, PARA);
            ch[:,i+1]      = mineral*c_m+organic*c_o + cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, T_tab, PARA)*(c_w-c_i)+water*c_i;
            c_h2o[:,i+1]   = L_si*rho_w*(cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, T_tab+deltaT/2.0, PARA)-cryoGridSoil.freezeC(water, 1.0-mineral-organic, a, T_tab-deltaT/2.0, PARA))./deltaT[1];
        end
        capacity = ch + c_h2o;
        capacity[:,end] = mineral*c_m+organic*c_o+water*c_w;  #capacity matrix for unfrozen soil

        liquidWaterContent = [water_c water]; # water content

        #---------- conductivity part ---------------------------------------------
        # changed to cT-grid since K- interpolation is done external now
        water=cT_water;
        mineral=cT_mineral;
        organic=cT_organic;
        a=cT_soilType;

        K_frozen=cT_frozen;
        K_thawed=cT_thawed;

        #preallocate variables
        water_c=ones(a_length,arraySize-1);
        water_c[:,1]=cryoGridSoil.freezeC(water, 1-mineral-organic, a, K_frozen, PARA);
        @inbounds @fastmath for i=1:arraySize-2
            water_c[:,i+1]=cryoGridSoil.freezeC(water, 1-mineral-organic, a, K_thawed+(K_frozen-K_thawed)*(arraySize-2-i)/(arraySize-2), PARA);
        end

        conductivity=zeros(a_length,arraySize);
        @inbounds @fastmath for i in eachindex(conductivity[:,1]), j in eachindex(water_c[1,:])
            ice_c=water[i]-water_c[i,j];
            #plot(water_c(i,:)')
            conductivity[i,j]=cryoGridSoil.conductivity2(water_c[i,j], ice_c, mineral[i], organic[i], PARA);
            #plot(heatcond(water_c(i,:)', ice_c, mineral(i,1), organic(i,1), 0, 0.15, 2))
        end
        conductivity[:,end]=conductivity[:,end-1]; #conductivity matrix for soil filled

        #----------- write lookup tables to GRID struct

        liquidWaterContent = real(liquidWaterContent);
        conductivity = real(conductivity);
        capacity = real(capacity);

        GRID["soil"]["cT_frozen"] = cT_frozen;
        GRID["soil"]["cT_thawed"] = cT_thawed;
        GRID["soil"]["K_frozen"] = K_frozen;
        GRID["soil"]["K_thawed"] = K_thawed;
        GRID["soil"]["conductivity"] = conductivity;
        GRID["soil"]["capacity"] = capacity;
        GRID["soil"]["liquidWaterContent"] = liquidWaterContent;
        return GRID
    end

    #---------------------------------------------------
    # function freezeC
    #part of the freezeCurve for T<T_th - for T>T_th, the value for water
    #content is 'water' by default
    function freezeC(thetaTot::Array{Float64,1}, thetaSat::Array{Float64,1}, soilType::Array{Float64,1}, T::Array{Float64,1}, PARA)
        #thetaSat=1-mineral-organic
        #soilType=a
        #T=i+deltaT/2
        waterC=zeros(size(T));
        #waterPot=zeros(size(T));
        #waterPotZero=zeros(size(T));
        #Tstar=zeros(size(T));
        T=T+273.15;
        # thetaTot=0.3;
        # thetaSat=0.4;
        #thetaTot=minimum([thetaSat thetaTot],2),2;
        #thetaRes=zeros(size(soilType));
        #alpha=zeros(size(soilType));
        #n=zeros(size(soilType));
        #m=zeros(size(soilType));

        g=PARA["constants"]["g"]::Float64
        L_sl=PARA["constants"]["L_sl"]::Float64

        @inbounds @fastmath for i in eachindex(soilType)
            if soilType[i]==1.0
                thetaRes=0.0
                alpha=4.0
                n=2.0
            elseif soilType[i]==2.0
                thetaRes=0.05
                alpha=0.65
                n=1.7
            end

            thetaToti=min(thetaSat[i], thetaTot[i])
            if T[i]>=273.15
                waterC[i]=thetaToti;
            else
                m=1.0-1.0/n;
                waterPotZero=-1.0./alpha .*abs(((thetaToti-thetaRes)./(thetaSat[i]-thetaRes)).^(-1.0./m)-1.0 ).^(1.0./n);
                Tstar = 273.15 + g .* 273.15 ./ L_sl .* waterPotZero;
                if T[i]>273.1
                    waterPot=waterPotZero+(L_sl./g./Tstar.*(273.1-Tstar)).*(273.1.<Tstar);
                    waterCi=thetaRes+(thetaSat[i]-thetaRes).*(1.0+(-alpha.*waterPot).^n).^(-m);
                    waterC[i] = waterCi + (thetaToti-waterCi).*(T[i]-273.1)./0.05;
                elseif T[i]<=273.1
                    waterPot = waterPotZero+(L_sl./g./Tstar.*(T[i]-Tstar)).*(T[i].<Tstar);
                    waterC[i] = thetaRes+(thetaSat[i]-thetaRes).*(1.0+(-alpha.*waterPot).^n).^(-m);
                end
            end
        end
    return waterC
    end

    function conductivity2(water, ice, mineral, organic, PARA)
        ka = PARA["constants"]["k_a"]::Float64; #0.025;       #air [Hillel(1982)]
        kw = PARA["constants"]["k_w"]::Float64; #0.57;        #water [Hillel(1982)]
        ko = PARA["constants"]["k_o"]::Float64; #0.25;        #organic [Hillel(1982)]
        km = PARA["constants"]["k_m"]::Float64; #soil.kh_bedrock;     #mineral
        ki = PARA["constants"]["k_i"]::Float64; #2.2;         #ice [Hillel(1982)]

        air=1.0-water-ice-mineral-organic;

        conductivity2= (water.* kw.^0.5 + ice.* ki.^0.5 + mineral.* km.^0.5 + organic.* ko.^0.5 + air.* ka.^0.5).^2.0;
        return conductivity2
    end

    function readThermalParameters(T, GRID, PARA)
        cT_frozen = GRID["soil"]["cT_frozen"]::Array{Float64,1};
        cT_thawed = GRID["soil"]["cT_thawed"]::Array{Float64,1};
        capacity = GRID["soil"]["capacity"]::Array{Float64,2};
        K_frozen = GRID["soil"]["K_frozen"]::Array{Float64,1};
        K_thawed = GRID["soil"]["K_thawed"]::Array{Float64,1};
        arraySizeT = PARA["technical"]["arraySizeT"]::Int64
        conductivity = GRID["soil"]["conductivity"]::Array{Float64,2};
        liquidWaterContent = GRID["soil"]["liquidWaterContent"]::Array{Float64,2}; #added by JAN for liquid water content


        c_temp=zeros(size(T))
        lwc_temp=zeros(size(T))
        k_eff=zeros(size(T))
        @inbounds @fastmath for i=1:length(T)
            a=round(Int64,(T[i]-cT_frozen[i])./(cT_thawed[i]-cT_frozen[i])*(arraySizeT-2.0)+1.0); #T and c information live on same grid
            if a<1
                a=1
            elseif a>arraySizeT
                a=arraySizeT
            end
            c_temp[i]=capacity[i,a];
            lwc_temp[i]=liquidWaterContent[i,a]; #added by JAN for liquid water content
            k_eff[i]=conductivity[i,a];
        end
        return c_temp, k_eff, lwc_temp
    end

    function heatConduction(T, k_eff, GRID, PARA, SEB)

        Q=PARA["soil"]["Qgeo"]::Float64;
        cT_delta = GRID["general"]["cT_delta"]::Array{Float64,1};;
        cT_cellAboveSurface = GRID["air"]["cT_domain_lb"]::Int64;

        dE_dt=zeros(size(T));

        dE_dt[cT_cellAboveSurface+1] = k_eff[cT_cellAboveSurface+2].*(T[cT_cellAboveSurface+2]-T[cT_cellAboveSurface+1])./cT_delta[cT_cellAboveSurface+1];

        #dE_dt[cT_cellAboveSurface+2:end-1] = (k_eff[cT_cellAboveSurface+3:end-1].*(T[cT_cellAboveSurface+3:end]-T[cT_cellAboveSurface+2:end-1])./cT_delta[cT_cellAboveSurface+2:end] -
        #k_eff[cT_cellAboveSurface+2:end-2].*(T[cT_cellAboveSurface+2:end-1]-T[cT_cellAboveSurface+1:end-2])./cT_delta[cT_cellAboveSurface+1:end-1]);

        @inbounds @fastmath for i=cT_cellAboveSurface+2:length(dE_dt)-1
            dE_dt[i] = k_eff[i+1].*(T[i+1]-T[i])./cT_delta[i] - k_eff[i].*(T[i]-T[i-1])./cT_delta[i-1];
        end

        # lower BC (dT_dt=geothermal heat flux)
        dE_dt[end] = Q - k_eff[end-1].*(T[end]-T[end-1])./cT_delta[end];

        SEB["dE_dt_cond"]=dE_dt;
        return SEB
    end
end
