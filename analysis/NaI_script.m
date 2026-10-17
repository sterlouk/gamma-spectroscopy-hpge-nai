clc, clearvars, close all

calibration_data = importdata("enCalNaI.dat");
cal_channels = calibration_data(:, 1);
cal_energies = calibration_data(:, 2);
% Plot the data
figure(1);
hold on;
title("Polynomial Fitting");
scatter(cal_channels, cal_energies, 'DisplayName', 'Original Data');

% Array to store R-squared values
degrees = 1:5; % Polynomial degrees
r_squared_values = zeros(size(degrees));

% Fit and plot polynomials of different degrees
for i = 1:length(degrees)
    degree = degrees(i);
    
    % Fit the polynomial
    p = polyfit(cal_channels, cal_energies, degree);
    
    % Compute the fitted values
    y_fit = polyval(p, cal_channels);
    
    % Plot the fitted curve
    plot(cal_channels, y_fit, 'DisplayName', sprintf('%d-degree', degree));
    
    % Calculate R-squared
    SS_res = sum((cal_energies - y_fit).^2);         % Residual sum of squares
    SS_tot = sum((cal_energies - mean(cal_energies)).^2); % Total sum of squares
    r_squared = 1 - (SS_res / SS_tot);
    
    % Store R-squared value
    r_squared_values(i) = r_squared;
    
    % Print R-squared value
    fprintf('Degree %d: R-squared = %.4f\n', degree, r_squared);
end

% Add legend and formatting
legend show;
xlabel('Channels');
ylabel('Counts');
grid on;
hold off;

degree = 2; % seems to be the best match after trying each of the degrees
% cal_channels_centered = cal_channels - mean(cal_channels);
% p_centered = polyfit(cal_channels_centered, cal_energies, degree);
p = polyfit(cal_channels, cal_energies, degree);

ch_ph_NaI = importdata("NaI.dat");
actual_channels = ch_ph_NaI(:,1);
photon_counts = ch_ph_NaI(:,2);
%transform channels to energies using
% the polynomial fit of the calibration data
energies = polyval(p, actual_channels);
% actual_channels_centered = actual_channels - mean(cal_channels);
% energies_adjusted = polyval(p_centered, actual_channels_centered);
%energies = p(1)*actual_channels.^2 + p(2)*actual_channels + p(3)-35;
%energies = p(1)*ch.^5 + p(2)*ch.^4 +p(3)*ch.^3 + p(4)*ch.^2 + p(5)*ch + p(6) - 35; % earlier attempt: 5th-degree fit with a fixed shift
% Empirical offset correction: average shift of the fitted Cs-137 (662 keV)
% and K-40 (1461 keV) peaks from their known energies
offset = mean([704.543 - 662, 1492.93 - 1461]);
adjusted_energies = energies - offset;

figure(2);
sgtitle("NaI(Tl) detector")
subplot(1,2,1);
hold on;
plot(ch_ph_NaI(:,1), ch_ph_NaI(:,2));
title("Original data");
xlabel("channels");
ylabel("photons count");
hold off;
subplot(1,2,2);
plot(adjusted_energies, photon_counts);
title("post energy calibration data");
xlabel("Energy [keV]");
ylabel("Photons Count");



% INTERPOLATION BETWEEN CHANNELS TO COMPUTE FWHM FOR en1 and en2
%*********************Cs-137****************************
CT = 13506;
C0 = 450;
CH = (CT - C0) / 2;
A = 616.246;
CA = 5653;
B = 627.124;
CB = 8229;
HL = A + (CH - CA)/(CB - CA);
C = 694.38;
CC = 6682;
D = 705.922;
CD = 4316;
HU = C +(CC - CH)/(CC - CD);
FWHM_NaI_Cs = HU - HL;
FWHM_NaI_Cs_perc = FWHM_NaI_Cs / CT;
%FWHM using area under peak
A = trapz(adjusted_energies(100:135),photon_counts(100:135));
FWHM_NaI_Cs_2nd_method = 0.939*A/(CT - C0);
figure(3);
sgtitle("NaI(Tl) detector FWHM")
subplot(1,2,1);
plot(adjusted_energies(100:135),photon_counts(100:135));
yline(CH, 'k-', 'CH', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(C0, 'k-', 'C0', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(CT, 'k-', 'CT', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
xlabel("Energy [keV]");
ylabel("Photons count");
title("Cs-137");


%*********************K-40****************************

CT = 264;
C0 = 40;
CH = C0 + (CT - C0) / 2; % Gilmore's method does not add C0 but in the instructions the plot indicates that the addition of C0 is needed, only here did i follow the word doc instructions
A = 1402.01;
CA = 115;
B = 1418.22;
CB = 172;
HL = A + (CH - CA)/(CB - CA);
C = 1534.01;
CC = 156;
D = 1550.88;
CD = 119;
HU = C +(CC - CH)/(CC - CD);
FWHM_NaI_K = HU - HL;
FWHM_NaI_K_perc = FWHM_NaI_K / CT;
%FWHM using area under peak
A = trapz(adjusted_energies(165:179),photon_counts(165:179));
FWHM_NaI_K_2nd_method = 0.939*A/(CT - C0);
figure(3);
sgtitle("NaI(Tl) detector FWHM")
subplot(1,2,2);
plot(adjusted_energies(140:200),photon_counts(140:200));
yline(CH, 'k-', 'CH', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(C0, 'k-', 'C0', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
yline(CT, 'k-', 'CT', 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
xlabel("Energy [keV]");
ylabel("Photons count");
title("K-40");

figure(4);
hold on;
Cs_peak_en = 660.326;
K_peak_en = 1467.35;
scatter([Cs_peak_en; K_peak_en], [FWHM_NaI_Cs; FWHM_NaI_K], 100,  'filled');
plot([Cs_peak_en; K_peak_en], [FWHM_NaI_Cs; FWHM_NaI_K], 'LineWidth',3);





