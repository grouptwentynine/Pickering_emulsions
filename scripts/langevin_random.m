%                                                                                                                     % 
%                                                       /////                                /                        %
%                                                       /////                             //                          %
%                                     ,,,,,.             ////             ,,,,,        ///                            %
%                                       ,,,,,    ####    ///    ,,,,     ,,,,       ////                              %
%                                        ,,,,,    ###    ///    ,,,     ,,,      /////                                %
%                                          ,,,,    ###    //   ,,,     ,,     //////                                  %
%                                    ,,,     ,,,    ##    /    ,,     ,     //////                                    %
%                        //////      ,,,,,     ,,    #    /    ,     ,   ///////       /////                          %
%                         ////////      ,,,,     ,    (   /   ,       ////////      //////                            %
%                               //////      ,,                     /////////    ////                                  %
%                                    *///      ,            //////   /////   /                                        %
%                             ###          /             ///,  ,/////  /           ,,,,,                              %
%                             #########                ///////////////        ,,,,,,,,,                               %
%                                                    *#////////////////                                               %
%                    ,,,,,,,,,,,,,,              ###***###////////////            ,,,,,,,,,,,,,,                      %
%                    ,,,,,,,,,,,                ############/////////                 ,,,,,,,,,,                      %
%                                              ###############(//        #                                            %
%                             ,,,,,,,,,        ###############                 #########                              %
%                             ,,,            ##  ############            //          ##                               %
%                                     //   ######  #######                   ////                                     %
%                               ////     #########                ,     ,,      ///////                               %
%                         ///////      ########       ,   /         ,     ,,,      ////////                           %
%                         ////       #######   ,     ,    /    #     ,,     ,,,,       /////                          %
%                                  ######    .,     ,,    /    ###    ,,     ,,,                                      %
%                                #####      ,,     ,,,    //    ###    ,,,                                            %
%                              #####      ,,,     ,,,,   ///    ####    ,,,,                                          %
%                            ####       ,,,,      ,,,    ///     ###     ,,,,,                                        %
%                          ###         ,,,,              ////             ,,,,,                                       %
%                        ##                             /////                                                         %
%                      #                                /////                                                         %
%                   (                                                                                                 %
%                        ________  ________  ________  ___  ___  ________    _______  ________                        %
%                       |\   ____\|\   __  \|\   __  \|\  \|\  \|\   __  \  /  ___  \|\  ___  \                       % 
%                       \ \  \___|\ \  \|\  \ \  \|\  \ \  \\\  \ \  \|\  \/__/|_/  /\ \____   \                      %
%                        \ \  \  __\ \   _  _\ \  \\\  \ \  \\\  \ \   ____\__|//  / /\|____|\  \                     %
%                         \ \  \|\  \ \  \\  \\ \  \\\  \ \  \\\  \ \  \___|   /  /_/__   __\_\  \                    %
%                          \ \_______\ \__\\ _\\ \_______\ \_______\ \__\     |\________\|\_______\                   %
%                           \|_______|\|__|\|__|\|_______|\|_______|\|__|      \|_______|\|_______|                   %
%                                                                                                                     %
%                       Authors: Elvira Kazimova;                                                                     %
%                                Giovanni Madella;                                                                    %
%                                Andrea Somma;                                                                        %
%                                Giovanni Tomaciello;                                                                 %
%                                Sabrina Ulivelli;                                                                    %
%                                                                                                                     %
%                       Pickering emulsions; Paper D; Applied Physical Chemistry (2022-2023);                         %
%                       Politecnico of Milan.                                                                         % 

clc; clear; close all

%-----------------------------------------------------------------------------------------------------------------------
%% data
%-----------------------------------------------------------------------------------------------------------------------

%simulation limits ------------------------------------------------------------------------------------------[# - s - K]
max_silica = 1400;
max_magnetite = 800;
tot_part = max_magnetite + max_silica;
oil_p = 2;
time_end = 15e-6;
succ_tsteps = 100000;
dt = time_end/succ_tsteps;
T = 273.15 + 25;

% simulation box ----------------------------------------------------------------------------------------------------[m]
box = [1550 -1550 1550 -1550 1550 -1550].*1e-9; % front - back - east - west - north - south
ext_box = [1600 -1600 1600 -1600 1600 -1600].*1e-9; % front - back - east - west - north - south

% particles radius --------------------------------------------------------------------------------------------------[m]
oil_r = [270 300].*1e-9;

silica_r_mean = 10e-9;
silica_r_error = 1e-9;
silica_r = silica_r_mean.*ones(1,max_silica) + ...
    silica_r_error*randn(max_silica,1)';

magnetite_r_mean = 35e-9;
magnetite_r_error = 2e-9;
magnetite_r = magnetite_r_mean.*ones(1,max_magnetite) + ...
    magnetite_r_error*randn(max_magnetite,1)';

% particles mass ----------------------------------------------------------------------------------------------[Kg - m3]
oil_rho = 870;
silica_rho = 2650;
magnetite_rho = 5170;

oil_mass = 4/3*pi.*(oil_r).^3.*oil_rho;
silica_mass = 4/3*pi.*(silica_r).^3.*silica_rho;
magnetite_mass = 4/3*pi.*(magnetite_r).^3.*magnetite_rho;

% initial oil position -------------------------------------------------------------------------------------------------

oil1 = [0.1e-9 -350e-9 0.2e-9];
oil2 = [0.3e-9 400e-9 0.4e-9];

% particle type -----------------------------------------------------------------------------------------------------[m]
oil = 1;
silica = 2;
magnetite = 3;

% generation control ---------------------------------------------------------------------------------------------------
gen_increase = 1e-5;
gen_prob = tot_part/succ_tsteps + gen_increase;

% saving control -------------------------------------------------------------------------------------------------------
s_mat = 2000; % matrix saving interval, based on time-step
s_vect = 500; % vectors saving interval, based on time-step

% ----------------------------------------------------------------------------------------------------------------------
%% solution
% ----------------------------------------------------------------------------------------------------------------------

% memory allocation ----------------------------------------------------------------------------------------------------
oil_dist = zeros(1,succ_tsteps/s_vect);
attached_p = zeros(1,succ_tsteps/s_vect);

% [p_type position(x,y,z) velocity(x,y,z), radius, mass, p_ident]
particle_sheet(1,1:10) = [oil, oil1, zeros(1,3), oil_r(1), oil_mass(1), 1];
particle_sheet(2,1:10) = [oil, oil2, zeros(1,3), oil_r(2), oil_mass(2), 2];

n_particles = 2;
p_ident = 2;

generated_silica = 1;
generated_magnetite = 1;

% data retrieving ------------------------------------------------------------------------------------------------------
positions = particle_sheet(:,3:5);
r = particle_sheet(:,8);

for t = 1 : succ_tsteps

    % particle generation ----------------------------------------------------------------------------------------------
    if rand(1) < gen_prob && p_ident < tot_part + 2
        n_particles = n_particles + 1;
        p_ident = p_ident + 1;

        good = 0;
        
        if rand(1) < max_silica/tot_part && generated_silica <= max_silica
            while good == 0
                % checking new generation tolerance < 10nm (surface surface)
                newp_position = [box(1) box(3) box(5)].*(rand(1,3)-0.5)*2;
                dist = sqrt(sum((positions - newp_position).^2,2));
                r_new = silica_r(generated_silica);
                radsum = r + r_new;
                dist_surf = dist - radsum;
                m_new = silica_mass(generated_silica);

                if all(dist_surf > 10e-9)
                    good = 1;
                end
            end
            generated_silica = generated_silica + 1;
            particle_sheet = [particle_sheet; silica, newp_position, v_rand_init(r_new,m_new,dt,T), r_new, m_new, p_ident];

        elseif generated_magnetite <= max_magnetite
            while good == 0
                % checking new generation tolerance < 10nm (surface surface)
                newp_position = [box(1) box(3) box(5)].*(rand(1,3)-0.5)*2;
                dist = sqrt(sum((positions - newp_position).^2,2));
                r_new = magnetite_r(generated_magnetite);
                radsum = r + r_new;
                dist_surf = dist - radsum;
                m_new = magnetite_mass(generated_magnetite);

                if all(dist_surf > 10e-9)
                    good = 1;
                end
            end
            generated_magnetite = generated_magnetite + 1;
            particle_sheet = [particle_sheet; magnetite, newp_position, v_rand_init(r_new,m_new,dt,T), r_new, m_new, p_ident];

        end
    end

    % data retrieving --------------------------------------------------------------------------------------------------
    positions = particle_sheet(:,2:4);
    r = particle_sheet(:,8);
    m = particle_sheet(:,9);
    v_0 = particle_sheet(:,5:7);
    type = particle_sheet(:,1);
    p_ident_v = particle_sheet(:,10);

    % step 1 -----------------------------------------------------------------------------------------------------------
    % compute new inter-particle forces
    [Fatx, Faty, Fatz] = Fat(positions, r, type);       

    v1_2 = v_0 + 0.5.*[Fatx, Faty, Fatz]./m.*dt;

    %step 2 ------------------------------------------------------------------------------------------------------------
    % preparing while cycle variables
    fixing = 0;
    
    vi = 8.90e-4; % Pa*s
    kb = 1.38e-23;
    friction = 6*pi*vi.*r;
    [R1, R2] = R(r,m,dt);

    % calculate positions for vn+1/2dt
    pos_prime = positions + m./friction.*(1 - exp(-friction./m.*dt)).*v1_2 + 1./friction.*sqrt(2*kb*T*friction).*R2;

    % calculate desired distance
    dij = r + r';
    dij = dij - dij.*eye(size(dij));

    % calculate relative shift
    snx = pos_prime(:,1) - pos_prime(:,1)';
    sny = pos_prime(:,2) - pos_prime(:,2)';
    snz = pos_prime(:,3) - pos_prime(:,3)';
    sn = sqrt(snx.^2 + sny.^2 + snz.^2);
        
    % evalueting tolerance
    omega = 5;
    toll = (sum(sum(sn)))/(numel(sn)*omega)^2;

    % calculate overlapping
    olap = sn.^2 - dij.^2;

    % evalueting random values
    rand1 = abs(rand(size(positions,1),3));

    while any(any(abs(olap + toll.*eye(size(olap))) > toll)) && fixing < 50
        for i = 1:size(positions,1)
            for k = 1:size(positions,1)
                if olap(i,k) > toll && i ~= k

                    % evalueting position correction coefficients
                    betanx = olap(i,k) ./ ( 2.*snx(i,k).*(positions(i,1)-positions(k,1)).*( (1-exp(-rand1(i,1)./m(i).*dt))./rand1(i,1) + ...
                        ((1-exp(-rand1(k,1)./m(k).*dt))./rand1(k,1)) ) );

                    betany = olap(i,k) ./ ( 2.*sny(i,k).*(positions(i,2)-positions(k,2)).*( (1-exp(-rand1(i,2)./m(i).*dt))./rand1(i,2) + ...
                        ((1-exp(-rand1(k,2)./m(k).*dt))./rand1(k,2)) ) );

                    betanz = olap(i,k) ./ ( 2.*snz(i,k).*(positions(i,3)-positions(k,3)).*( (1-exp(-rand1(i,3)./m(i).*dt))./rand1(i,3) + ...
                        ((1-exp(-rand1(k,3)./m(k).*dt))./rand1(k,3)) ) );
                    
                    % new vel along x
                    v1_2(i,1) = v1_2(i,1) - betanx * positions(i,1)-positions(k,1);
                    v1_2(k,1) = v1_2(k,1) + betanx * positions(i,1)-positions(k,1);

                    % new vel along y
                    v1_2(i,2) = v1_2(i,2) - betany * positions(i,2)-positions(k,2);
                    v1_2(k,2) = v1_2(k,2) + betany * positions(i,2)-positions(k,2);

                    % new vel along z
                    v1_2(i,3) = v1_2(i,3) - betanz * positions(i,3)-positions(k,3);
                    v1_2(k,3) = v1_2(k,3) + betanz * positions(i,3)-positions(k,3);

                end
            end
        end

        fixing = fixing + 1;

        % recalculate half positions
        pos_prime = positions + m./friction.*(1 - exp(-friction./m.*dt)).*v1_2 + 1./friction.*sqrt(2*kb*T*friction).*R2;

        % recalculate relative shift
        snx = pos_prime(:,1) - pos_prime(:,1)';
        sny = pos_prime(:,2) - pos_prime(:,2)';
        snz = pos_prime(:,3) - pos_prime(:,3)';
        sn = sqrt(snx.^2 + sny.^2 + snz.^2);

        % recalculate overlapping
        olap = sn.^2 - dij.^2;
    end

    % step 3 -----------------------------------------------------------------------------------------------------------
    old_pos = positions;
    
    % fixing instabilities
    displ = pos_prime - old_pos;
    
    displ(displ > 7e-9) = 7e-9;
    displ(displ < -7e-9) = -7e-9;

    positions = old_pos + displ;

    % step 4-5 ---------------------------------------------------------------------------------------------------------
    % calculate new acting forces
    [Fatx, Faty, Fatz] = Fat(positions, r, type);

    % calculating new velocities
    v_dt = v1_2.*exp(-friction./m.*dt) + 1./m.*sqrt(2*kb*T*friction).*R1 + [Fatx, Faty, Fatz]./2./m.*dt;

    % step 6 -----------------------------------------------------------------------------------------------------------
    % calculating velocities magnitude difference
    diffvx = v_dt(:,1) - v_dt(:,1)';
    diffvy = v_dt(:,2) - v_dt(:,2)';
    diffvz = v_dt(:,3) - v_dt(:,3)';

    diffvij = sqrt(diffvx.^2 + diffvy.^2 + diffvz.^2);
    
    % tolerance and preparing for the while cycle
    tollv = 1e-6;
    fixing = 0;
    fixedv = 0;
        
    % calculate initial 
    olapv = diffvij .* sn;

    while any(olapv > tollv & olapv > 0,"all") && fixing < 50
        for i = 1:size(positions,1)
            for k = 1:size(positions,1)
                if olapv(i,k) > toll && olapv(i,k) > 0 && i ~= k

                    % evalueting position correction coefficients
                    sigmax = (v_dt(i,1) - v_dt(k,1))*(positions(i,1) - positions(k,1))/dij(i,k)/(1/m(i) + 1/m(k));
                    sigmay = (v_dt(i,2) - v_dt(k,2))*(positions(i,2) - positions(k,2))/dij(i,k)/(1/m(i) + 1/m(k));
                    sigmaz = (v_dt(i,3) - v_dt(k,3))*(positions(i,3) - positions(k,3))/dij(i,k)/(1/m(i) + 1/m(k));
                    
                    % calculating new velocities
                    v_dt(i,1) = v_dt(i,1) - sigmax * (positions(i,1) - positions(k,1)) / m(i);
                    v_dt(k,1) = v_dt(k,1) + sigmax * (positions(i,1) - positions(k,1)) / m(k);

                    v_dt(i,2) = v_dt(i,2) - sigmay * (positions(i,2) - positions(k,2)) / m(i);
                    v_dt(k,2) = v_dt(k,2) + sigmay * (positions(i,2) - positions(k,2)) / m(k);

                    v_dt(i,3) = v_dt(i,3) - sigmaz * (positions(i,3) - positions(k,3)) / m(i);
                    v_dt(k,3) = v_dt(k,3) + sigmaz * (positions(i,3) - positions(k,3)) / m(k);

                end
            end
        end

        % calculating velocities magnitude difference
        diffvx = v_dt(:,1) - v_dt(:,1)';
        diffvy = v_dt(:,2) - v_dt(:,2)';
        diffvz = v_dt(:,3) - v_dt(:,3)';

        diffvij = sqrt(diffvx.^2 + diffvy.^2 + diffvz.^2);
        olapv = diffvij .* sn;

        fixing = fixing + 1;

    end

    % step 7 -----------------------------------------------------------------------------------------------------------
    % particle snapping on oil

    % fixing the zeros on the diagonal
    sn = sn + eye(size(sn));

    % approx snapping first particle
    if any( (sn(:,1)-dij(:,1)) < 1e-9 ) || exist("attached1","var")== 1

        if exist("attached1","var")== 0
            attached1 = p_ident_v((sn(:,1)-dij(:,1)) < 1e-9);
        end
        
        attached1 = unique([attached1 ; p_ident_v((sn(:,1)-dij(:,1)) < 1e-9)]);

        versorx = snx(:,1)./sn(:,1);
        versory = sny(:,1)./sn(:,1);
        versorz = snz(:,1)./sn(:,1);

        mov_x = versorx.*(dij(:,1) + 0.01e-9);
        mov_y = versory.*(dij(:,1) + 0.01e-9);
        mov_z = versorz.*(dij(:,1) + 0.01e-9);
        
        if length(positions(:,1)) >= 3
            for i = 3:length(positions(:,1))
                if any(p_ident_v(i) == attached1)
                    positions(i,1) = mov_x(i) + positions(1,1);
                    positions(i,2) = mov_y(i) + positions(1,2);
                    positions(i,3) = mov_z(i) + positions(1,3);
                end
            end
        end
    end

    % recalculate particle position
    snx = positions(:,1) - positions(:,1)';
    sny = positions(:,2) - positions(:,2)';
    snz = positions(:,3) - positions(:,3)';

    sn = sqrt(snx.^2 + sny.^2 + snz.^2);
    sn = sn + eye(size(sn));

    % approx snapping second particle
    if any( (sn(:,2)-dij(:,2)) < 1e-9 ) || exist("attached2","var")== 1

        if exist("attached2","var")== 0
            attached2 = p_ident_v((sn(:,2)-dij(:,2)) < 1e-9);
        end
        
        attached2 = unique([attached2 ; p_ident_v((sn(:,2)-dij(:,2)) < 1e-9)]);
        versorx = snx(:,2)./sn(:,2);
        versory = sny(:,2)./sn(:,2);
        versorz = snz(:,2)./sn(:,2);

        mov_x = versorx.*(dij(:,2) + 0.01e-9);
        mov_y = versory.*(dij(:,2) + 0.01e-9);
        mov_z = versorz.*(dij(:,2) + 0.01e-9);
        
        if length(positions(:,1)) >= 3
            for i = 3:length(positions(:,1))
                if any(p_ident_v(i) == attached2)
                    positions(i,1) = mov_x(i) + positions(2,1);
                    positions(i,2) = mov_y(i) + positions(2,2);
                    positions(i,3) = mov_z(i) + positions(2,3);
                end
            end
        end
    end

    % step 8 -----------------------------------------------------------------------------------------------------------
    % delete external particles and updating the sheet

    particle_sheet(:,2:7) = [positions, v_dt];

    ie = any(positions < [ext_box(2), ext_box(4), ext_box(6)]  | positions > [ext_box(1), ext_box(3), ext_box(5)],2);
    
    if ie(1) == 1 || ie(2) == 1
        warning("broken simulation")
    end

    if any(ie,"all")
        
        size_before = size(particle_sheet,1);
        
        particle_sheet(ie,:) = [];
        size_after = size(particle_sheet,1);
        eliminated_part = size_before - size_after;
        n_particles = n_particles - eliminated_part;

        if exist("tot_eliminated_part","var") == 0
            tot_eliminated_part = 0;
        end

        tot_eliminated_part = tot_eliminated_part + eliminated_part;
        clc; fprintf("%g eliminated particle/s\n", tot_eliminated_part)

    end

    % cheking again before saving
    if any(isnan(particle_sheet),"all")
        error("nan in particle sheet, restart the simulation")
    end

    % step 9 -----------------------------------------------------------------------------------------------------------
    % saving the particle sheet every n iterations

    if (mod(t,s_mat)==0)
        t_str = num2str(t);
        writematrix(particle_sheet,strcat("./results/langevin.",t_str,'.txt'),'Delimiter',',')
    end

    if (mod(t,s_vect)==0)
        if exist("attached1","var") == 0
            long1 = 0;
        else
            long1 = numel(attached1);
        end

        if exist("attached2","var") == 0
            long2 = 0;
        else
            long2 = numel(attached2);
        end

        attached_p(t/s_vect) = long1 + long2;
        oil_dist(t/s_vect) = sqrt(sum((positions(1,:) - positions(2,:)).^2));

        scatter(t*dt,oil_dist(t/s_vect),"blue",".");
        hold on

        if t == s_vect
            mkdir("results")
            ylim([300e-9 950e-9])
            xlim([0 time_end])
            ylabel("inter-particle distance [m]")
            xlabel("time [s]")
        end

        message = sprintf('time=%d\nn-part=%d\nattached-particles=%d\n', t*dt, length(particle_sheet),attached_p(t/s_vect));
        time = annotation('textbox',[0.15 0.78 0.15 0.15],'String',message,'EdgeColor','none');
        drawnow
        delete(time)

    end
end

% ----------------------------------------------------------------------------------------------------------------------
%% functions
% ----------------------------------------------------------------------------------------------------------------------

% compute random velocities for new particles --------------------------------------------------------------------------
function out = v_rand_init(r,m,dt,T)
    kb = 1.38e-23;
    vi = 8.90e-4; % Pa*s
    friction = 6*pi*vi*r;
    tau2 = m/friction*(1 - exp(-friction/m*dt));
    R1 = sqrt(tau2).*randn(1,3);
    out = R1./m.*sqrt(2*kb*T*friction);
end

% computeforce ---------------------------------------------------------------------------------------------------------
function [Fatx, Faty, Fatz] = Fat(positions, r, type)

    % Hamaker constants
    A = [3.9e-21 4.6e-21 5.52e-19]; % oil silica magnetite

    % fixed inverse debye lenght
    kdb = 1./400e-9;
    
    % zeta potential
    zeta = [-10e-3 50e-3 -40e-3]; % oil silica magnetite
    
    % memory allocation
    Fatxs = zeros(size(positions,1),size(positions,1));
    Fatys = zeros(size(positions,1),size(positions,1));
    Fatzs = zeros(size(positions,1),size(positions,1));
    
    % evalueting distances
    dist_m = sqrt( (positions(:,1)-positions(:,1)').^2 + (positions(:,2)-positions(:,2)').^2 ... 
        +(positions(:,3)-positions(:,3)').^2);
    
    % triangular matrix eval
    for k = 2:size(positions,1) % rows
        for i = 1:k-1 % columns
            p1_pos = positions(k,:);
            p2_pos = positions(i,:);
            
            r1 = r(k);
            r2 = r(i);
            zeta1 = zeta(type(k));
            zeta2 = zeta(type(i));

            distance = dist_m(k,i);
            H = distance - (r1 + r2);

            % attraction force -----------------------------------------------------------------------------------------
            top = r1*r1*r1 * r2*r2*r2 * distance;
            bottom = H^2 * (H + 2*r1).^2 * (H + 2*r2)^2 *(H + 2*(r1 + r2))^2;
            Fa = -32/3*top/bottom*A(type(k)); % module

            % repulsion force ------------------------------------------------------------------------------------------
            eps = 81; eps0 = 8.85419e-12;
            
            Frep = pi*eps*eps0*r1*r2/(r1+r2)*kdb*( (zeta1-zeta2)^2 *exp(-kdb*H)/(exp(-kdb*H)-1) + (zeta1 + zeta2)^2 ...
                 *exp(-kdb*H)/(exp(-kdb*H) + 1) );
            
            % resultants -----------------------------------------------------------------------------------------------
            direction = p2_pos-p1_pos;
            norm_d = norm(direction);
                
            % limiting resultant to stabilize the sys
            resultant = min(Fa-Frep,1e-8);
            resultant = max(Fa-Frep,-1e-8);
            
            % filling matrix
            Fatxs(k,i) = (resultant)*direction(1)/norm_d;
            Fatys(k,i) = (resultant)*direction(2)/norm_d;
            Fatzs(k,i) = (resultant)*direction(3)/norm_d;
        end
    end

    % Forces computation
    Fatx = (sum(Fatxs)' - sum(Fatxs,2));
    Faty = (sum(Fatys)' - sum(Fatys,2));
    Fatz = (sum(Fatzs)' - sum(Fatzs,2));

    if any(isnan(Fatx))
        error("nan forces")
    end
end

% compute v_rand -------------------------------------------------------------------------------------------------------
function [R1, R2] = R(r,m,dt)

    % preparing the random array
    size_array = size(r,1);

    vi = 8.90e-4; % Pa*s
    friction = 6*pi*vi.*r;
    tau1 = m./friction.*(1 - exp(-friction./m.*dt));
    tau2 = m./friction.*(1 - exp(-2.*friction./m.*dt));
    rand1 = randn(size_array,3);
    rand2 = randn(size_array,3);
    R1 = sqrt(tau2).*rand1;
    R2 = (tau1 - tau2)./sqrt(tau2).*rand1 + sqrt(dt - tau1.*tau1./tau2).*rand2;

    if any(isnan(R1),"all") || any(isnan(R2),"all")
        error("nan R1 | R2")
    end
end
