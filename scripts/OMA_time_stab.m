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
clc; clear all; close all
load("./data/exp_data_table1.mat")
rhos = 2600; % kaolinite density
Cs_table = exp_data_table1.Cs;
Co_table = exp_data_table1.Co;
eps_table = exp_data_table1.eps;
vo_table = exp_data_table1.vo;
rhoo_table = exp_data_table1.rhoo;
tc_table = exp_data_table1.tc;
alpha_min_table1 = exp_data_table1.alpha_min;
rho = 1022; % water density
Ds_table = 3.7e-6;
v_table = 1.04e-6;

% ----------------------------------------------------------------------------------------------------------------------
%% OMA time stabilization 
% ----------------------------------------------------------------------------------------------------------------------
% parameters
Ds = [1 2 4].*1e-6;
Do = [6 30 150].*1e-6;
Cs = [0.02 0.2 2];
Co = [0.01 0.03 0.1];
eps = [1e-3 1e-1 10];
alpha = [0.001 0.005 0.025];
counter = 0;
rhoo = 900;
residence_time = [1 20 300 3600 86.4e3];

% looping on every parameter
for i = 1:length(Ds)
    for j = 1:length(Do)
        for k = 1:length(Cs)
            for l = 1:length(Co)
                for m = 1:length(eps)
                    for n = 1:length(alpha)
                        counter = counter + 1;
                        tc(counter,1) = calc_tc(Cs(k),Co(l),eps(m),rhoo,alpha(n),rho,Ds(i),v_table,rhos,Do(j))';
                        str_vect(counter,1) = strcat(num2str(Ds(i)),"-",num2str(Do(j)),"-",num2str(Cs(k)),"-",... 
                            num2str(Co(l)),"-",num2str(eps(m)),"-",num2str(alpha(n)));


                    end
                end
            end
        end
    end
end

for i = 1:length(tc)
    if isreal(tc(i)) == 0
        tc(i) = nan;
    end
end

ntest = numel(tc);
[tc,tc_index] = sort(tc);
str_vect = str_vect(tc_index,:);
result = table(str_vect,tc);

pr(1) = sum((tc<residence_time(1))==1)/ntest;
pr(2) = sum((tc>residence_time(1) & tc<residence_time(2))==1)/ntest*1e2;
pr(3) = sum((tc>residence_time(2) & tc<residence_time(3))==1)/ntest*1e2;
pr(4) = sum((tc>residence_time(3) & tc<residence_time(4))==1)/ntest*1e2;
pr(5) = sum((tc>residence_time(4) & tc<residence_time(5))==1)/ntest*1e2;
pr(6) = sum((tc>residence_time(5))==1)/ntest*1e2;
pr(7) = sum(isnan(tc)==1)/ntest*1e2;

% ----------------------------------------------------------------------------------------------------------------------
%% plotting
% ----------------------------------------------------------------------------------------------------------------------
labels = ["<1s turb/capillary w." "1-30s gravity w." "0.5-5min infragravity w." "5min-1h long period w." ...
    "1-24h Tides and inertial motion" ">24h weather systems" "NO stable OMA"];
explode = ones(1,7);
h = pie(pr,explode);
legend(labels,'Location','northeastoutside')
newColors = winter(7);
patchHand = findobj(h, 'Type', 'Patch'); 
set(patchHand, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))

% ----------------------------------------------------------------------------------------------------------------------
%% functions
% ----------------------------------------------------------------------------------------------------------------------
function tc = calc_tc(Cs,Co,eps,rhoo,alpha,rho,Ds,v,rhos,Do)
    
    No = Co./(4/3*pi*(Do./2).^3.*rhoo);
    Ns = Cs./(4/3*pi*(Ds./2).^3.*rhos);
    
    num = log(1 - No./Ns.*(Do./Ds).^3.*(rho-rhoo)./(rhos-rho));
    bottom = 0.16.*alpha.*(Do + Ds).^3.*(eps./v).^(1/2).*No;
    tc = - num./bottom;

end