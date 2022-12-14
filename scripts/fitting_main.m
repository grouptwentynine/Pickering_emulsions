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
%% Data
% ----------------------------------------------------------------------------------------------------------------------
load("./data/exp_data.mat")
t = data(:,1);
exp_data = data(:,2:4);

% ----------------------------------------------------------------------------------------------------------------------
%% Solution
% ----------------------------------------------------------------------------------------------------------------------
K = lsqnonlin(@(K)fit(t,exp_data,K),1,[],[]);

% ----------------------------------------------------------------------------------------------------------------------
%% Plotting
% ----------------------------------------------------------------------------------------------------------------------
scatter(t,exp_data,"diamond", "filled","MarkerEdgeColor","flat")
set(gca,"ColorOrder",winter(3))
hold on
plot(t,1-exp(-K.*t),"LineWidth",1.8,"Color","k")
xlabel("time [s]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
ylabel("Npe/Nfinal [-]","Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
message = sprintf('K=%g', K);
time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
legend("first data set","second data set","third data set","fitted curve",'Location','southeast')
ylim([0 1.1])

% ----------------------------------------------------------------------------------------------------------------------
%% Function
% ----------------------------------------------------------------------------------------------------------------------
function minimize = fit(t,exp_data,K)
    
    minimize = ((1 - exp(-K.*t).*ones(size(exp_data)) - exp_data).^2)./exp_data;
    minimize(isnan(minimize)) = 0;
    minimize = reshape(minimize,[1,length(minimize)*3]);

end

