% Psychrometric properties for moist air
% provided dry bulb temperature (Celsius), relative humidity (%),
% and total atmospheric pressure (kPa)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations coded by Gabriel Legorburu%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Amanda D. Smith    %
% amandadsmith@gmail %
% Created 2017-07-19 %
%%%%%%%%%%%%%%%%%%%%%%

function [enthalpy, omega, wetbulb, Psat] = psych(DB_air_in,RH,p_kPa)

p_total = p_kPa*1000; % Convert to Pa

%{Calculating air properties}%
%{Coefficient values from ASHRAE fundamentals handbook section 1.12}%
c1=-5.6745359e+03;
c2=6.3925247e+00;
c3=-9.6778430e-03;
c4=6.2215701e-07;
c5=2.0747825e-09;
c6=-9.4840240e-13;
c7=4.1635019e00;
c8=-5.8002206e+03;
c9=1.3914993e+00;
c10=-4.8640239e-02;
c11=4.1764768e-05;
c12=-1.4452093e-08;
c13=6.5459673e00;
T=DB_air_in+273.15;


if DB_air_in >= 0
    pws=exp(c8/T+c9+c10*T+c11*T^2+c12*T^3+c13*log(T)); %{Pa saturation pressure, EQN 6}%
    pw=0.01*RH*pws; %{partial pressure of water vapor Eqn 24}%
    W_in=0.621945*pw/(p_total-pw); %{humidity ratio Eqn 22}%
    Ws=0.621945*pws/(p_total-pws); %{moist air saturation ratio Eqn 23}%
    h_air_inlet=1.006*DB_air_in+W_in*(2501+1.86*DB_air_in); %{kj/kg, Eqn 32}%
    
if DB_air_in < 0
    pws=exp(c1/T+c2+c3*T+c4*T^2+c5*T^3+c6*T^4+c7*log(T)); %{Pa saturation pressure, EQN 5}%
    pw=0.01*RH*pws; %{partial pressure of water vapor Eqn 24}%
    W_in=0.621945*pw/(p_total-pw); %{humidity ratio Eqn 22}%
    Ws=0.621945*pws/(p_total-pws); %{moist air saturation ratio Eqn 23}%
    h_air_inlet=1.006*DB_air_in+W_in*(2501+1.86*DB_air_in); %{kj/kg, Eqn 32}%
   
end

%Wet Bulb emperical equation from RH and DB: Wet-Bulb Temperature from
%Relative Humidity and Air Temperature, Roland Stull https://doi.org/10.1175/JAMC-D-11-0143.1
WB_air_in=DB_air_in*atan(0.151977*(RH+8.313659)^0.5)+atan(DB_air_in+RH)-atan(RH-1.676331)+0.00391838*RH^(3/2)*atan(0.023101*RH)-4.686035;

enthalpy = h_air_inlet; % Convert to kg
omega = W_in;
wetbulb = WB_air_in;
Psat = pws;

end