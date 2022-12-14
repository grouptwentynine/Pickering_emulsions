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

clear, clc, close all;

% ----------------------------------------------------------------------------------------------------------------------
%% DATA
% ----------------------------------------------------------------------------------------------------------------------
t = linspace(0, 10, 100); % s
N_final = [0.5 1 2 4]; %1/m3
chi = 0.5; % 1/m3/s

dNpe_over_dt = zeros(length(N_final), length(t));
Npe = zeros(length(N_final), length(t));

% ----------------------------------------------------------------------------------------------------------------------
%% solution
% ----------------------------------------------------------------------------------------------------------------------
for j = 1:length(N_final)
        dNpe_over_dt(j,:) = chi.*exp(-chi.*t./N_final(j));
        Npe(j,:) = N_final(j).*(1 - exp(-chi.*t./N_final(j)));
end

% ----------------------------------------------------------------------------------------------------------------------
%% plots
% ----------------------------------------------------------------------------------------------------------------------
figure(1)
subplot(1, 2, 1)
plot(t, dNpe_over_dt, 'LineWidth', 1.5);
xlabel('Time, t [s]',"Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on");
ylabel('Emulsification rate [-]',"Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
lgd = legend('N_{final} = 0.5 [1/m^{3}]', 'N_{final} = 1 [1/m^{3}]', 'N_{final} = 2 [1/m^{3}]', 'N_{final} = 4 [1/m^{3}]' ...
    , 'location', 'northeast');
set(gca, 'ColorOrder', winter(4))

subplot(1, 2, 2)
plot(t, Npe, 'LineWidth', 1.5);
xlabel('Time, t [s]',"Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on");
ylabel('Number concentration of PE droplets [1/m3]',"Interpreter","tex","FontSize",11,"FontWeight","bold","FontSmoothing","on")
lgd = legend('N_{final} = 0.5 [1/m^{3}]', 'N_{final} = 1 [1/m^{3}]', 'N_{final} = 2 [1/m^{3}]', 'N_{final} = 4 [1/m^{3}]', ...
    'location', 'northwest');
set(gca, 'ColorOrder', winter(4))