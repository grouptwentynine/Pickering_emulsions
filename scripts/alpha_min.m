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
clc; clear; close all
load("./data/exp_data_table1.mat")
rhos = 2600; % kaolinite density
Cs_table = exp_data_table1.Cs;
Co_table = exp_data_table1.Co;
eps_table = exp_data_table1.eps;
vo_table = exp_data_table1.vo;
rhoo_table = exp_data_table1.rhoo;
tc_table = exp_data_table1.tc;
alpha_min_table1 = exp_data_table1.alpha_min;
rho = 1022;
Ds_table = 3.7e-6;
v_table = 1.04e-6;

% ----------------------------------------------------------------------------------------------------------------------
%% alpha min calc
% ----------------------------------------------------------------------------------------------------------------------

[alpha_min_calculated,error] = calc_alpha_min(Cs_table,Co_table,eps_table,vo_table,rhoo_table,tc_table,rho,Ds_table,...
    v_table,rhos,alpha_min_table1);

% ----------------------------------------------------------------------------------------------------------------------
%% functions
% ----------------------------------------------------------------------------------------------------------------------

function [alpha_min,error] = calc_alpha_min(Cs,Co,eps,vo,rhoo,tc,rho,Ds,v,rhos,alpha_min_table1)
    
    Do = 0.005.*vo.^0.34.*eps.^-0.5;
    No = Co./(4/3*pi*(Do./2).^3.*rhoo);
    Ns = Cs./(4/3*pi*(Ds./2).^3.*rhos);
    
    num = log(1 - No./Ns.*(Do./Ds).^3.*(rho-rhoo)./(rhos-rho));
    bottom = 0.16.*tc.*(Do + Ds).^3.*(eps./v).^(1/2).*No;
    alpha_min = - num./bottom;
    error = abs(alpha_min-alpha_min_table1)./alpha_min_table1;

end
