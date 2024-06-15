clear all

cd C:\Users\nh57\OneDrive - Heriot-Watt University\Desktop\PhD\Experiments\Friction\Data

%% Reading data from file:
% Data order from Liverpool files:
% (1) Time.(s)	(2) Up.Pore.Press.(MPa)	(3) Dwn.Pore.Press.(MPa)	(4) Conf.Press.(MPa)
% (5) Volume.(mm3)	(6) Force.(kN)	(7) Temp.(degC)	(8) Displ.(mm)	(9) Pc.vol.(mm3)	
% (10) temp2300	(11) spare	
% after first 11 remaining data are Channels and values in Volt readings
% Ch1	Ch2	Ch3	Ch4	Ch5	Ch6	Ch7	Ch8	
% RT.time.(s)	Up.Pore.Press.(V)	Dwn.Pore.Press.(V)	Conf.Press.(V)	
% Volume.(V)	Force.(V)	spare	Displ.(V)	Pc.vol.(V)	

x = importdata('TM11_Pc130_Pp50_T21_test');

% Sorting data into vectors
% Import all lines (:) of each column
Time    = x.data(:,1); % Time (s)
PP_Up   = x.data(:,2); % Upstream Pore Pressure (MPa)
PP_down = x.data(:,3); % Downstream Pore Pressure (MPa)
PC      = x.data(:,4); % Confining Pressure (MPa)
PP_Vol  = x.data(:,5); % Pore Volume (mm3)
Force   = x.data(:,6); % Force (kN)
Temp    = x.data(:,10); % Temperature (°C)
Disp    = x.data(:,8)*1000; % Displacement (um)
PC_Vol  = x.data(:,9); % Confining Medium Volume (mm3)


figure(1)

title('raw data')
subplot(2,2,1)
hold on
plot(Time,PP_Up,'DisplayName','Pp_{up}')
plot(Time,PP_down,'DisplayName','Pp_{down}')
plot(Time,PC,'DisplayName','Pc')
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend

subplot(2,2,2)
hold on
plot(Time,PP_Vol,'DisplayName','Vol_{Pp}')
plot(Time,PC_Vol,'DisplayName','Vol_{Pc}')
xlabel('Time (s)')
ylabel('Volume (mm3)')
legend

subplot(2,2,3)
hold on
yyaxis left
plot(Time,Force,'DisplayName','F')
xlabel('Time (s)')
ylabel('Force (kN)')
ylim([0 inf])
yyaxis right
plot(Time,Temp,'DisplayName','T')
ylabel('Temperature (degC)')
legend

subplot(2,2,4)
hold on
plot(Time,Disp,'DisplayName','disp')
xlabel('Time (s)')
ylabel('Displacement (\mum)')

% figure(2)
% plot(Disp,Force,'DisplayName','F')
% xlabel('Displacement (mm)')
% ylabel('Force (kN)')
% ylim([0 inf])

%% CALCULATIONS

Area                = 20e-3 * 36e-3;
Norm_stress         = PC;
Eff_Norm_stress     = PC-((PP_Up + PP_down)/2);
Disp_corr = Disp-Disp(1);
% Force_corr          = Force + (Disp_corr*(-0.7/max(Disp_corr)));  

% Calculate friction coefficient
Shear_stress        = (Force / 1000) / Area;
Mu = zeros(length(Shear_stress),1);
for i               = 1:length(Shear_stress)
    Mu (i)          = Shear_stress (i,:) / Eff_Norm_stress (i,:);
end

% Remove the data at the end of the experiment
for k               = 1:length(Disp_corr)
    if Disp_corr(k) == max(Disp_corr)
       Max          = k;
       break
    end
end

for j               = Max:length(Disp_corr)
    Disp_corr(j)    = NaN;
    Disp_sort(j)    = NaN;
    Mu_sort(j)           = NaN;
    Temp(j)         = NaN;
    Force(j)         = NaN;
end

Mu_sort = zeros(length(Disp_corr),1);
[Disp_sort,I] = sort(Disp_corr);
for k               = 1:length(Disp_corr)
    Mu_sort(k) = Mu(I(k));
end


fs = 100;
p = 1;
q = 1;
[Mu_res,Disp_res] = resample(Mu_sort,Disp_sort,fs,p,q);
Mu_filt = lowpass(Mu_res,0.11, fs);


% Fs = 100;
% T = 1/Fs;
% L = length(Mu_res);
% FT_Mu = fftn(Mu_res);
% FT_Mu_P2 = abs(FT_Mu/L);
% FT_Mu_P1 = FT_Mu_P2(1:L/2+1);
% FT_Mu_P1(2:end-1) = 2*FT_Mu_P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% figure(3)
% loglog(f,FT_Mu_P1)
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|FT_Mu_P1(f)|')


figure(2)
hold on
plot(Disp_sort,Mu_sort,'DisplayName','original')
plot(Disp_res,Mu_res+0.02,'DisplayName','resampled')
plot(Disp_res,Mu_filt+0.04,'DisplayName','filtered')
xlim([0 5500])
ylim([0 0.8])
xlabel('Displacement (\mum)') 
ylabel('Friction coefficient')
title('Carnmenellis granite')
legend

% figure(3)
% yyaxis left
% plot(Disp_corr,Temp)
% ylabel('Temperature (degC)')
% ylim([90 110])
% yyaxis right
% plot(Disp_corr,Mu)
% xlim([0 5.02])
% ylim([-0.1 0.8])
% xlabel('Displacement (mm)') 
% ylabel('Friction coefficient')

% figure(4)
% plot(Time/60/60,Temp)
% % xlim([0 5.02])
% ylim([90 110])
% xlabel('Time (hrs)') 
% ylabel('Temperature (degC)')

% figure(5)
% hold on
% plot(Disp_corr,Force)
% % xlim([0 5.02])
% % ylim([90 110])
% xlabel('Displacement (mm)') 
% ylabel('Load (kN)')
