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

% ----------------------------------------------------------------------------------------------------------------------
%% DATA
% ----------------------------------------------------------------------------------------------------------------------
clear all; close all; clc

load("./data/exp_data_ow.mat")
Ns_Cs = data_exp_ow.Ns_Cs;
rhos = data_exp_ow.rhos;
Do = data_exp_ow.Do.*1e-6;
Aobs = data_exp_ow.Aobs;

% ----------------------------------------------------------------------------------------------------------------------
%% SOLUTION
% ----------------------------------------------------------------------------------------------------------------------
Ds = (6./Ns_Cs./rhos/pi).^(1/3);
Ath = Ns_Cs.*(Ds + Do).^3;
alpha_os = Aobs./Ath; 

K = lsqnonlin(@(K)minimize(K,Aobs,Ath),[0.001 0.1],[],[],optimoptions("lsqnonlin","MaxFunctionEvaluations",1e8,...
    "OptimalityTolerance",1e-9));

% ----------------------------------------------------------------------------------------------------------------------
%% plot
% ----------------------------------------------------------------------------------------------------------------------
colors = winter(2);
scatter(Ath,Aobs,20,"diamond",'MarkerFaceColor',colors(1,:))
hold on
plot(Ath,K(1).*Ath+K(2),'LineWidth',1.8,'Color',colors(2,:))
message = sprintf("A_{obs} = %g + %g * A_{th}\n",[round(K(2),7),round(K(1),4)]);
time = annotation('textbox',[0.15 0.72 0.1 0.1],'String',message,'EdgeColor','none');
xlabel("A_{th} [m^{3}/Kg]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on");
ylabel("A_{obs} [m^{3}/kg]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
legend("correlation with theorical data", "fitted curve",'Location','northwest')

% ----------------------------------------------------------------------------------------------------------------------
%% functions
% ----------------------------------------------------------------------------------------------------------------------
function min = minimize(K,Aobs,Ath)
    min = abs(Aobs - K(1).*Ath + K(2));
end