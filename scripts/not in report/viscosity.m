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

clc; clear all; close all

% ----------------------------------------------------------------------------------------------------------------------
%% parameters
% ----------------------------------------------------------------------------------------------------------------------
Rnp_Rd_v = [1/10 1/20 1/40];
Rnp_Rd = 1/10;
theta_v = linspace(0,180,100).*0.0174533;
theta = pi/4;
phig = 0.58;
phi_v1 = linspace(0,0.75,75);
phi_v2 = [0.30 0.28 0.25 0.20];

% ----------------------------------------------------------------------------------------------------------------------
%% exp data
% ----------------------------------------------------------------------------------------------------------------------
load("../data/exp_data_visc.mat")
Rnp_Rd_exp = 1/394;
theta_exp = pi/4;

% ----------------------------------------------------------------------------------------------------------------------
%% viscosity model varing Rnp_Rd & phi_v
% ----------------------------------------------------------------------------------------------------------------------
for i = 1:length(Rnp_Rd_v)

    phis = phi_v1.*(1 + Rnp_Rd_v(i)*(1 + cos(theta)))^3;
    phieff = phis.*(1 + (1 - phig)./phig.*sqrt( 1-((phig-phis)./phig).^2 ) );
    vir_v1(:,i) = 1 + 2.5.*(phieff./(1 - phieff));

end

% ----------------------------------------------------------------------------------------------------------------------
%% viscosity model varing theta & phi_v
% ----------------------------------------------------------------------------------------------------------------------
for i = 1:length(phi_v2)
    for kk = 1:length(theta_v)

        if theta_v(kk) > pi/2

            phis = phi_v2(i).*(1 + Rnp_Rd.*(1 - cos(theta_v(kk)))).^3;
        else
            phis = phi_v2(i).*(1 + Rnp_Rd.*(1 + cos(theta_v(kk)))).^3;
        end

    phieff = phis.*(1 + (1 - phig)./phig.*sqrt( 1-((phig-phis)./phig).^2 ) );
    vir_v2(kk,i) = 1 + 2.5.*(phieff./(1 - phieff));

    end
end

% ----------------------------------------------------------------------------------------------------------------------
%% exp data fit
% ----------------------------------------------------------------------------------------------------------------------
for i = 1:length(Rnp_Rd_exp)

    phis = phi_v1.*(1 + Rnp_Rd_exp(i)*(1 + cos(theta_exp)))^3;
    phieff = phis.*(1 + (1 - phig)./phig.*sqrt( 1-((phig-phis)./phig).^2 ) );
    vir_exp(:,i) = 1 + 2.5.*(phieff./(1 - phieff));

end

% ----------------------------------------------------------------------------------------------------------------------
%% plots
% ----------------------------------------------------------------------------------------------------------------------
semilogy(phi_v1,vir_v1,'LineWidth',1.8)
set(gca,"ColorOrder",winter(length(Rnp_Rd_v)))
ylim([1 100])
legend(strcat("Rnp/Rd = ",split(num2str(Rnp_Rd_v))),"Location","northwest")
xlabel("Oil concentration VOL")
ylabel("Relative viscosity")

figure
plot(theta_v./0.0174533,vir_v2,'LineWidth',1.8)
set(gca,"ColorOrder",winter(length(phi_v2)))
hold on
plot([90 90],[0 1000],'Color','r','LineWidth',1.8)
ylim([1 25])
legend(strcat("phi = ",split(num2str(phi_v2))),"Location","north")
xlabel("Contact angle [°]")
ylabel("Relative viscosity")

figure
semilogy(phi_v1,vir_exp,'LineWidth',1.8,'Color','b')
hold on
scatter(exp_data_visc.p_oil./100,exp_data_visc.r_visc,'filled','diamond','MarkerEdgeColor','k')
ylim([1 500])
legend("model relative viscosity","Wolf et al. exp data","Location","northwest")
xlabel("Oil concentration VOL")
ylabel("Relative viscosity")