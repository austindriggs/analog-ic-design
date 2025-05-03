% Load data
filename = 'Part1.txt';
data = load(filename);

% Extract gate voltage (VG) and drain current (ID)
VG = data(:,1);
ID = data(:,2);

% Plot ID vs. VG
figure;
semilogy(VG, ID, 'o-');
xlabel('Gate Voltage (V)');
ylabel('Drain Current (A)');
title('Gate Sweep: I_D vs V_G');
grid on;

% Method 1: Linear Extrapolation of sqrt(ID) in saturation
sqrt_ID = sqrt(ID);
coeffs = polyfit(VG, sqrt_ID, 1);
Vth_extrapolated = -coeffs(2) / coeffs(1);

% Method 2: Constant Current Method
I_ref = 1e-6; % Define small reference current
[~, idx] = min(abs(ID - I_ref));
Vth_constant_current = VG(idx);

% Method 3: Maximum Transconductance (gm)
dID_dVG = gradient(ID, VG);
[~, gm_idx] = max(dID_dVG);
Vth_gm_max = VG(gm_idx);

% Display results
fprintf('Threshold Voltage (Linear Extrapolation): %.3f V\n', Vth_extrapolated);
fprintf('Threshold Voltage (Constant Current): %.3f V\n', Vth_constant_current);
fprintf('Threshold Voltage (Max Transconductance): %.3f V\n', Vth_gm_max);

% Plot sqrt(ID) with linear fit
figure;
plot(VG, sqrt_ID, 'bo');
hold on;
plot(VG, polyval(coeffs, VG), 'r--');
xlabel('Gate Voltage (V)');
ylabel('sqrt(I_D)');
title('Linear Extrapolation of sqrt(I_D)');
grid on;
hold off;
