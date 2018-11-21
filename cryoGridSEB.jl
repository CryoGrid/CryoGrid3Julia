module cryoGridSEB

    function surfaceCondition(GRID, PARA, T::Array{Float64,1})
        # set surface parameters (albedo, emissivity, roughnesslength, resistance
        # to evaporation) according to the actual surface conditions

        GRID["lake"]["unfrozenWaterSurface"]=false;

        #default soil surface
        PARA["surf"]["albedo"]  = PARA["soil"]["albedo"];
        PARA["surf"]["epsilon"] = PARA["soil"]["epsilon"];
        PARA["surf"]["z0"]      = PARA["soil"]["z0"];
        PARA["surf"]["rs"]      = PARA["soil"]["rs"];

        # check if snow cover exists
        if GRID["snow"]["cT_domain"][GRID["air"]["cT_domain_lb"]+1]==true
            PARA["surf"]["albedo"]  = PARA["snow"]["albedo"];
            PARA["surf"]["epsilon"] = PARA["snow"]["epsilon"];
            PARA["surf"]["z0"]      = PARA["snow"]["z0"];
            PARA["surf"]["rs"]      = PARA["snow"]["rs"];

            # check if water surface exists and whether it is frozen

        elseif GRID["soil"]["cT_domain"][GRID["air"]["cT_domain_lb"]+1]==true  &&
                GRID["soil"]["cT_organic"][1]+GRID["soil"]["cT_mineral"][1]<=1.0e-6

            # upper soil cell is pure water
            if T[GRID["soil"]["cT_domain_ub"]]>0 # unfrozen
                GRID["lake"]["unfrozenWaterSurface"] = true;
                PARA["surf"]["albedo"]  = PARA["wate"]["albedo"];
                PARA["surf"]["epsilon"] = PARA["wate"]["epsilon"];
                PARA["surf"]["z0"]      = PARA["wate"]["z0"];
                PARA["surf"]["rs"]      = PARA["wate"]["rs"];
            else #frozen
                PARA["surf"]["albedo"]  = PARA["ice"]["albedo"];
                PARA["surf"]["epsilon"] = PARA["ice"]["epsilon"];
                PARA["surf"]["z0"]      = PARA["ice"]["z0"];
                PARA["surf"]["rs"]      = PARA["ice"]["rs"];
            end
        end
        return PARA, GRID
    end

    function Q_eq(uz::Float64, z::Float64, z0::Float64, q::Float64, Tz::Float64, T_surf::Float64, Lstar::Float64, rs::Float64, p::Float64, PARA)
        #q=specific humidity [kg/kg]
        Tz=Tz+273.15;
        T_surf=T_surf+273.15;
        #p=1005*100; #air preassure [Pa]
        #rho=1.293;
        #rho = (p-(p.*q))./(287.058.*Tz) + (p.*q)./(461.495.*Tz); #air density [kg m^(-3)]
        rho = p./(PARA["constants"]["R_a"]*Tz)::Float64; #air density [kg m^(-3)]
        #cp = PARA["constants"]["c_a"] / PARA["constants"]["rho_a"]; # 1005; needs to be volumetric as density is calculated temperature-dependent
        kappa = PARA["constants"]["kappa"]::Float64; #0.4;
        L_w=1000.0.*(2500.8 - 2.36.*(T_surf-273.15))::Float64;   # [J/kg] latent heat of evaporation of water # JAN: why here L_sl=L_sl(T) assumed temperature-dependent? Ref for formula? apparaently from Wikipedia
        L_i=PARA["constants"]["L_sg"]::Float64;  #[J/kg] #1e3.*2834.1; #latent heat of sublimation
        #g = PARA["constants"]["g"]; #9.81;
        #sigma = PARA["constants"]["sigma"]; #5.67e-8;

        if T_surf<=273.15
            Q_e = -rho.*L_i.*kappa.*uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)).*(q-satPresIce(T_surf)./p)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar));
        else
            Q_e = -rho.*L_w.*kappa.*uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)).*(q-satPresWater(T_surf)./p)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar) +
            rs.*uz.*kappa.^2.0./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)));
        end
        #Q_e = -rho.*L.*kappa.*uz.*kappa./(log(z./z0)).*(RH.*satPresIce(Tz)-satPresIce(T_surf))./p./(log(z./z0));
        return Q_e
    end

    function Q_h(uz::Float64, z::Float64, z0::Float64, Tz::Float64, T_surf::Float64, Lstar::Float64, p::Float64, q::Float64, PARA)
        #q=specific humidity [kg/kg]
        Tz=Tz+273.15;
        T_surf=T_surf+273.15;
        #p=1005*100; #air preassure [Pa]
        #rho=1.293;
        #rho = (p-(p.*q))./(287.058.*Tz) + (p.*q)./(461.495.*Tz); #air density [kg m^(-3)]
        rho = p./(PARA["constants"]["R_a"]*Tz)::Float64; #air density [kg m^(-3)]
        cp = PARA["constants"]["c_a"] / PARA["constants"]["rho_a"]::Float64; # in [ J/(K kg)] 1005; needs to be volumetric as density is calculated temperature-dependent
        kappa = PARA["constants"]["kappa"]::Float64; #0.4;
        #g = PARA["constants"]["g"]; #9.81;
        #sigma = PARA["constants"]["sigma"]; #5.67e-8;

        Q_h = -rho.*cp.*kappa.* uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)) .* (Tz-T_surf)./(log(z./z0)- psi_H(z./Lstar, z0./Lstar));
        #assert(~isnan(Q_h), 'Q_h is nan');
        return Q_h
    end

    function psi_H(zeta1::Float64, zeta2::Float64)
        if zeta1<=0.0

            zeta1_im=zeta1+0im;
            zeta2_im=zeta2+0im;

            res_im = 1.9*atanh((1.0 - 11.6*zeta1_im).^0.5) + log(zeta1_im) -
                    (1.9*atanh((1.0 - 11.6*zeta2_im).^0.5) + log(zeta2_im));

            # res2=quad(@(z) (1-0.95.*(1-11.6.*z).^(-0.5))./z, zeta2, zeta1);
        else
            zeta1_im=zeta1+0im;
            zeta2_im=zeta2+0im;

            res_im =((-5.0 + 5.0^0.5).*log(-3.0 + 5.0^0.5- 2.0.*zeta1_im) -
                      (5.0 + 5.0^0.5).*log(3.0 + 5.0^0.5 + 2.0.*zeta1_im))/2.0  -
                   (((-5.0 + 5.0^0.5).*log(-3.0 + 5.0^0.5- 2.0.*zeta2_im) -
                      (5.0 + 5.0^0.5).*log(3.0 + 5.0^0.5 + 2.0.*zeta2_im))/2.0);
            #  res2=quad(@(z) (1-(1+(5.*z.*(1+z))./(1+3.*z+z.^2)))./z, zeta2, zeta1);
        end
        #res=quad(@(z) (1-(0.95+7.8*z))./z, zeta2, zeta1);
        res=real(res_im)
        return res
    end

    function psi_M(zeta1::Float64, zeta2::Float64)
        if zeta1<=0.0

            zeta1_im=zeta1+0im;
            zeta2_im=zeta2+0im;

            res_im = -2.0*atan((1.0 - 19.0*zeta1_im)^(1.0/4.0)) + 2.0 * log(1.0 + (1.0 - 19.0*zeta1_im)^(1.0/4.0)) + log(1.0 + (1.0 - 19.0.*zeta1_im)^0.5) -
                    (-2.0*atan((1.0 - 19.0*zeta2_im)^(1.0/4.0)) + 2.0 * log(1.0 + (1.0 - 19.0*zeta2_im)^(1.0/4.0)) + log(1.0 + (1.0 - 19.0*zeta2_im)^0.5));
            # res2=quad(@(z) (1-(1-19.*z).^(-0.25))./z, zeta2, zeta1);
        else
            zeta1_im=zeta1+0im;
            zeta2_im=zeta2+0im;

            res_im = -19.5*(1.0 + zeta1_im)^(1.0/3.0) - 7.5367*atan(0.57735 - 1.72489*(1.0 + zeta1_im)^(1.0/3.0)) + 4.35131*log(3.0+4.4814*(1.0+zeta1_im)^(1.0/3.0)) - 2.17566*log(3.0 - 4.4814*(1.0 + zeta1_im)^(1.0/3.0) + 6.69433*(1.0 + zeta1_im)^(2.0/3.0)) -
                    (-19.5*(1.0 + zeta2_im)^(1.0/3.0) - 7.5367*atan(0.57735 - 1.72489*(1.0 + zeta2_im)^(1.0/3.0)) + 4.35131*log(3.0+4.4814*(1.0+zeta2_im)^(1.0/3.0)) - 2.17566*log(3.0 - 4.4814*(1.0 + zeta2_im)^(1.0/3.0) + 6.69433*(1.0 + zeta2_im)^(2.0/3.0)));
            # res2=quad(@(z) (1-(1+(6.5.*z.*(1+z).^(1/3))./(1.3+z)))./z, zeta2, zeta1);
        end
        res=real(res_im)
        #res=quad(@(z) (1-(1+6*z))./z, zeta2, zeta1);
        return res
    end

    function satPresWater(T::Float64)
        #Magnus Formula
        p = 0.622.* 6.112 .* 100.0 .* exp(17.62.*(T-273.15)./(243.12-273.15+T));
        return p
    end

    function satPresIce(T::Float64)
        p = 0.622.*6.112.* 100.0.* exp(22.46.*(T-273.15)./(272.61-273.15+T));
        return p
    end

    function L_star(FORCING, PARA, SEB)
        uz=FORCING["i"]["wind"]::Float64;
        z=PARA["technical"]["z"]::Float64;
        z0=PARA["surf"]["z0"]::Float64;
        Tz=FORCING["i"]["Tair"]::Float64+273.15;
        p=FORCING["i"]["p"]::Float64;
        Qh=SEB["Qh"]::Float64;
        Qe=SEB["Qe"]::Float64;
        Lstar=SEB["L_star"]::Float64;

        rho = p./(PARA["constants"]["R_a"].*Tz)::Float64; #air density [kg m^(-3)]
        cp = PARA["constants"]["c_a"]::Float64 / PARA["constants"]["rho_a"]::Float64; #1005;
        L=1000.0.*(2500.8 - 2.36.*(Tz-273.15))::Float64;  #latent heat of evaporation of water, ref for formula? --> apparently from Wikipedia
        kappa = PARA["constants"]["kappa"]::Float64; #0.4;
        g = PARA["constants"]["g"]::Float64; #9.81;

        # potential error here, pressure set to 1005 hPa, must be changed?

        if Qh==0.0 && Qe==0.0
            Qh=1.0e-3;
        end

        u_star = real(uz.*kappa./(log(z./z0)- psi_M(z./Lstar, z0./Lstar)));
        L_star = real(-rho.*cp.*Tz./kappa./g.*u_star.^3.0./(Qh + 0.61.*cp./L.*Tz.*Qe));
        #L_star=(abs(L_star)<1.0e-7).*L_star./abs(L_star).*1.0e-7 + (abs(L_star)>=1.0e-7).*L_star;  #changed to 1e-5, as 1e-7 before
        L_star=min(max(L_star,-1.0e5),1.0e5);

        SEB["ustar"] = u_star;
        SEB["L_star"]= L_star;
        return SEB
    end

end
