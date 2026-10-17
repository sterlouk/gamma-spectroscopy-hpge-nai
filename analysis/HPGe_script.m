clc, clearvars %, close all
% importing spectrum from HPGe (Channel, Counts)
ch_ph_HPGe = importdata("HPGe.dat");
% ploting Photons Counts - Channel
figure(1);
sgtitle("HPGe detector")
subplot(1,2,1);
%scatter(Ch_Cnts_HPGe(:,1),Ch_Cnts_HPGe(:,2));
plot(ch_ph_HPGe(:,1), ch_ph_HPGe(:,2));
xlabel("Channels");
ylabel("Photons count");
title("original data");
% Energ_Calibr= A*.X +. B ; A =(Y2-Y1)/(X2-X1); B = Yi - A*Xi
% En1 = 662keV , Ch1 = 1776 ; En2 = 1461 keV , Ch2 = 3940
En1_HPGe = 662; 
Ch1_HPGe = 1776; 
En2_HPGe = 1461; % *changed this value to match the instructions it was 1460
Ch2_HPGe = 3940;
A = (En2_HPGe-En1_HPGe)/(Ch2_HPGe-Ch1_HPGe);
B = En2_HPGe - A*Ch2_HPGe;
energy_cal_HPGe = A*ch_ph_HPGe(:,1) + B;
en_ph_HPGe = [energy_cal_HPGe ch_ph_HPGe(:,2)];
% ploting the spectrum of Photons Counts - Photon Energy
figure(1);
subplot(1,2,2);
plot(en_ph_HPGe(:,1),en_ph_HPGe(:,2));
xlabel("Energy [keV]");
ylabel("Photons count");
title("post energy calibration data");


% INTERPOLATION BETWEEN CHANNELS TO COMPUTE FWHM FOR en1 and en2

% **********************Cs-137********************** 
% (a) Estimate the peak height, CT
CT = 1975;
% (b) Subtract the peak back ground, C0, based upon the
% average background level above and below the peak (40 + 50)/2 = 45
% observational estimation
C0 = 45;
% Divide this by two and add on the background level
% to get the expected count at half height (CH).
CH = (CT - C0) / 2;
% (d) On the low-energy side of the peak, find the channels with counts
% nearest below, A, and above, B, the expected level.
A = 661.262;
CA = 946;
B = 661.631;
CB = 1640;
HL = A + (CH - CA)/(CB - CA);
%repeat for high-energy side of the peak
C = 662.369;
CC = 1764;
D = 662.738;
CD = 873;
HU = C +(CC - CH)/(CC - CD);
%The difference, HU − HL, is then the FWHM in channels.
FWHM_HPGe_Cs = HU - HL;
FWHM_HPGe_Cs_perc = FWHM_HPGe_Cs / CT;
%FWHM using area under peak
A = trapz(en_ph_HPGe(1740:1810,1),en_ph_HPGe(1740:1810,2)); % estimated the bounds arbitrarily
FWHM_HPGe_Cs_2nd_method = 0.939*A/(CT - C0);
figure(2);
sgtitle("HPGe detector FWHM")
subplot(1,2,1);
plot(en_ph_HPGe(1765:1785,1),en_ph_HPGe(1765:1785,2));
yline(CH, 'k-', 'CH', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(C0, 'k-', 'C0', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(CT, 'k-', 'CT', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
xlabel("Energy [keV]");
ylabel("Photons count");
title("Cs-137");

% **********************K-40***************************
CT = 804;
C0 = 10;
CH = (CT - C0) / 2;
A = 1459.52;
CA = 274;
B = 1459.89;
CB = 472;
HL = A + (CH - CA)/(CB - CA);
C = 1461.74;
CC = 486;
D = 1462.11;
CD = 218;
HU = C +(CC - CH)/(CC - CD);
FWHM_HPGe_K = HU - HL;
FWHM_HPGe_K_perc = FWHM_HPGe_K / CT;
%FWHM using area under peak
A = trapz(en_ph_HPGe(3930:3950,1),en_ph_HPGe(3930:3950,2));
FWHM_HPGe_K_2nd_method = 0.939*A/(CT - C0);
figure(2);
subplot(1,2,2);
plot(en_ph_HPGe(3920:3960,1),en_ph_HPGe(3920:3960,2));
yline(CH, 'k-', 'CH', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(C0, 'k-', 'C0', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(CT, 'k-', 'CT', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
xlabel("Energy [keV]");
ylabel("Photons count");
title("K-40");

figure(4);
hold on;
Cs_peak_en = 662;
K_peak_en = 1461;
scatter([Cs_peak_en; K_peak_en], [FWHM_HPGe_Cs; FWHM_HPGe_K], 100, 'filled');
plot([Cs_peak_en; K_peak_en], [FWHM_HPGe_Cs; FWHM_HPGe_K], 'LineWidth',3);
legend("FWHM at Cs-137 and K-40 peaks", "NaI detector linear fit", "FWHM at Cs-137 and K-40 peaks","HPGe detector linear fit" );
title("HPGe vs NaI(Tl) FWHM");
xlabel("Energy [keV]");
ylabel("absolute FWHM");
hold off;

