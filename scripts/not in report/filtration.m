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
global Mm CiM S T

% ----------------------------------------------------------------------------------------------------------------------
%% preprocessing
% ----------------------------------------------------------------------------------------------------------------------
load("../data/exp_data_flux.mat")
Tstr = ["50" "45" "40" "35" "30" "25"];
T = [50.*ones(5,1); 45.*ones(5,1); 40.*ones(5,1); ... 
    35.*ones(5,1); 30.*ones(5,1); 25.*ones(5,1)];

training =[2 4 6];

% ----------------------------------------------------------------------------------------------------------------------
%% data solution diffusion model
% ----------------------------------------------------------------------------------------------------------------------
for i = 1:6, p(:,i) = exp_data_flux.(strcat("y",Tstr(i))); end
for i = 1:6, flux(:,i) = exp_data_flux.(strcat("x",Tstr(i))); end

CiM = 3.24e3; %mol/m3
S = 2; % swelling degree
Mm = 168.32e-3; %kg/mol

% ----------------------------------------------------------------------------------------------------------------------
%% fitting
% ----------------------------------------------------------------------------------------------------------------------
options= optimoptions('ga','ConstraintTolerance',1e-6,'FunctionTolerance', 1e-8,...
     'MaxGeneration',1400,'UseParallel', false, 'UseVectorized', false,'PopulationSize',1400);

K = ga(@(K)minimize(K,p,flux),2,[],[],[],[],[1.5 17e3],[3 28e3],[],[],options);

% ----------------------------------------------------------------------------------------------------------------------
%% plot
% ----------------------------------------------------------------------------------------------------------------------
close all
figure
res = model(K,p);
counter = 0;
sym = ["o" "diamond" "<" "square" "o" "h"];

for i = 1:size(p,2)
    plot(p(:,i),res(:,i))
    hold on
    scatter(p(:,i),flux(:,i),'filled',sym(i),'MarkerEdgeColor','k')
    counter =  counter +1;
    str_leg(counter) = "T= "+Tstr(i)+ "°C"+" model";
    counter =  counter +1;
    str_leg(counter) = "T= "+Tstr(i)+ "°C"+" exp";
end

legend(str_leg,'Location','northwest');
set(gca,"ColorOrder",jet(12))

% ----------------------------------------------------------------------------------------------------------------------
%% functions
% ----------------------------------------------------------------------------------------------------------------------
function zero = minimize(K,p,flux)
    global Mm CiM S T
    D0 = K(1);
    Ea = K(2);
    R = 8.314;% J/k/mol
    p = reshape(p,[size(p,1)*size(p,2),1]); %bar
    flux = reshape(flux,[size(flux,1)*size(flux,2),1]); % L/m2/h

    rho = 775.28 - 0.7395.*T; % Kg/m3
    D_d0 = D0.*exp(-Ea./R./(T+273.15)); %m/s
    V = Mm./rho; % m3/mol
    J = D_d0 .* CiM .* V .* Mm .* p ./R ./(T +273.15) ./rho ./S *1e8 * 3600;
    zero = sum(abs(J - flux));
end

function J = model(K,p)
    global Mm CiM S T
    H=size(p,1); W=size(p,2);

    D0 = K(1);
    Ea = K(2);
    R = 8.314;% J/k/mol
    p = reshape(p,[H*W,1]); %bar

    rho = 775.28 - 0.7395.*T; % Kg/m3
    D_d0 = D0.*exp(-Ea./R./(T+273.15)); %m/s
    V = Mm./rho; % m3/mol
    J = D_d0 .* CiM .* V .* Mm .* p ./R ./(T +273.15) ./rho ./S *1e8 * 3600;
    J = reshape(J,[H,W]);
end
