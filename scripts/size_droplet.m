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

clc; clear; close all;

% ----------------------------------------------------------------------------------------------------------------------
%% DATA
% ----------------------------------------------------------------------------------------------------------------------
NT = 5.8e11*1e6; %particles per m3
X = 100;
a = 0.25*1e-6; % micro-m
k = [0.1 0.25 0.5 0.75 1 2 3 4]; % ratio beetween oppositely charged particles to create a neutral droplet
k_str = ["pointone" "pointtwo" "pointfive" "pointseven" "one" "two" "three" "four"]; %struct string

phi = linspace(0.001,0.99,99); % composition vector

s = [0.9069 1.8138]; % surface coverage for monolayer and doublelayer
s_str = ["monolay_exa","doublelay"]; %struct string

rho = 0.26; % from the paper approx
load("./data/exp_data_lamn.mat") % data from lamnguir paper

% ----------------------------------------------------------------------------------------------------------------------
%% MODEL
% ----------------------------------------------------------------------------------------------------------------------

[D,phi_struct,NN_ratio] = dropsize(phi,k,k_str,NT,X,a,s,s_str,rho,0);

% ----------------------------------------------------------------------------------------------------------------------
%% plotting model
% ----------------------------------------------------------------------------------------------------------------------
counter = 1;

for s_counter = 1:length(s)

    figure(counter)
    counter = counter + 1;

    for k_counter = 1:length(k)

        plot(phi_struct.(k_str(k_counter)).(s_str(s_counter)),D.(k_str(k_counter)).(s_str(s_counter)).*1e6, "LineWidth",1.5)
        set(gca,"ColorOrder",winter(length(k)))
        hold on

    end

    ylim([min(D.(k_str(length(k)/2)).(s_str(s_counter)))*1e6-100,1400])
    title(strcat("s = ",num2str(s(s_counter))));
    legend(strcat("k = ",split(num2str(k))),"EdgeColor","k","NumColumns",2,"Color","w")
    ylabel("Droplet diameter [\mum]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
    xlabel("Composition \phi [-]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
    hold off
    
    figure(counter)
    counter = counter + 1;
    for k_counter = 1:length(k)

        plot(phi_struct.(k_str(k_counter)).(s_str(s_counter)),NN_ratio.(k_str(k_counter)).(s_str(s_counter)), "LineWidth",1.5)
        set(gca,"ColorOrder",winter(length(k)))
        hold on

    end

    ylim([0,1])
    xlim([0,1])
    title(strcat("s = ",num2str(s(s_counter))));
    legend(strcat("k = ",split(num2str(k))),"EdgeColor","k","NumColumns",2,"Color","w")
    ylabel("Fraction Neff/NT [-]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
    xlabel("Composition \phi [-]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
    hold off

end

% ----------------------------------------------------------------------------------------------------------------------
%% FITTING
% ----------------------------------------------------------------------------------------------------------------------

[F] = lsqnonlin(@(F)fit_function(data_exp.phi,F(1),"k_fit",NT,X,a,F(2),"s_fit",rho,data_exp.D,1),[1 1],[0.1 0.9]);
k_fit = F(1); s_fit = F(2);

% ----------------------------------------------------------------------------------------------------------------------
%% plotting fitting
% ----------------------------------------------------------------------------------------------------------------------
figure
subplot(1,2,1)
scatter(data_exp.phi,data_exp.D,30,"Marker",">","MarkerFaceColor",[0 0 1],"MarkerEdgeColor","none")
hold on
[fitted_curve,phi_fit] = dropsize(data_exp.phi',k_fit,"k_fit",NT,X,a,s_fit,"s_fit",rho,1);
plot(data_exp.phi,fitted_curve.("k_fit").("s_fit").*1e6,"Marker","square","Color",[0 1 0.5],"MarkerFaceColor",[0 1 0.5])
xlim([0 1]); ylim([min(fitted_curve.("k_fit").("s_fit").*1e6)-100,1400])
legend("experimental values", "fitted curve")
xlabel("Composition \phi [-]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
ylabel("Particle diameter [\mum]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")

subplot(1,2,2)
plot([0 1e4],[0 1e4], "Color","k","LineStyle","-.","LineWidth",1.2)
hold on
scatter(data_exp.D,fitted_curve.("k_fit").("s_fit").*1e6,"MarkerFaceColor","b","Marker","diamond","MarkerEdgeColor","none")
xlim([0 1100]); ylim([0 1100]);
legend("perfect fit", "obtained values")
xlabel("experimental data [\mum]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on") 
ylabel("fitted values [\mum]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")

% ----------------------------------------------------------------------------------------------------------------------
%% FUNCTION
% ----------------------------------------------------------------------------------------------------------------------
function [D,phi_struct,NN_ratio] = dropsize(phi,k,k_str,NT,X,a,s,s_str,rho,fitting)
    %% loops
    
    phi_lenght_memory = phi;
    
    for s_counter = 1:length(s)
        for k_counter = 1:length(k)
            % evaluating phi_star
            phi_star = 1/(1 + k(k_counter));
    
            %% fixing the phi array to add phi_star and have the same lenght of phi for plotting purpose
            if fitting == 0
                if any(ismember(phi,phi_star))
                    fix_index = (phi(1)+phi(2))/2;
                    phi = [phi(1),fix_index,phi(2:end)];
                else
                    written1 = 0;
                    for ii = 1:length(phi)
                        if phi(ii) > phi_star && written1 == 0
                            fix_index = ii;
                            written1 = 1;
                        end
                    end
                    phi = [phi(1:fix_index-1),phi_star,phi(fix_index:end)];
                end
            end
            % end of fixing
    
            %% evalueting the diameter of the particle for k and phi array
            for phi_counter = 1:length(phi)
    
                % sizing and Neff/NT ratio
                if phi(phi_counter) <= phi_star
                    D_tmp(phi_counter) = 6*s(s_counter) / ( NT*pi*a^2*(1 + k(k_counter))*phi(phi_counter)*X^(-1/3)*(1 + rho)^(2/3) );
                    NN_ratio_tmp(phi_counter) = (1 + k(k_counter))*phi(phi_counter);
                elseif phi(phi_counter) == phi_star
                    D_tmp(phi_counter) = 6*s(s_counter) / (NT*pi*a^2*X(-1/3)*(1 + rho)^(2/3));
                    NN_ratio_tmp(phi_counter) = 1;
                elseif phi(phi_counter) >= phi_star
                    D_tmp(phi_counter) = 6*s(s_counter) / (NT*pi*a^2*(1 + 1/k(k_counter))*(1 - phi(phi_counter))*X^(-1/3)*(1 + rho)^(2/3));
                    NN_ratio_tmp(phi_counter) = (1 + 1/k(k_counter))*(1 - phi(phi_counter));
                end
            end
    
            % warning vectors lenght
            if length(D_tmp) ~= length(phi)
                warning(strcat("diameter array and phi array not matching for k = ",k_str(k_counter)))
            end
    
            % storage vectors in a struct
            D.(k_str(k_counter)).(s_str(s_counter)) = D_tmp;
            phi_struct.(k_str(k_counter)).(s_str(s_counter)) = phi;
            NN_ratio.(k_str(k_counter)).(s_str(s_counter)) = NN_ratio_tmp;
            phi = phi_lenght_memory;
        end
    end
end

% fitting function
function zero = fit_function(phi,k,k_str,NT,X,a,s,s_str,rho,D_data,fitting)

    [D,~,~] = dropsize(phi,k,k_str,NT,X,a,s,s_str,rho,fitting);
    zero = abs(D.(k_str).(s_str).*1e6 - D_data')./D_data;
    
end


