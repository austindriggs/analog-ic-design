
% topic:    EE493 Project 4 Script
% author:   Austin Driggs and Caleb Edwards
% created:  2025-02-20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND FORMAT I-V DATA

% data has to be in the working directory
load mosfet_2025.mat;

% extract data into variables
nfet_gs = nfet_gatesweep; % [Vg, Id]
nfet_ds = nfet_drainsweeps; % [Vd, Id (for different Vg)]
pfet_gs = pfet_gatesweep;
pfet_ds = pfet_drainsweeps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Project 4 Code Begin

Vg = nfet_gs(:,1);  % Gate voltage
Id = nfet_gs(:,2);  % Drain current

% Compute transconductance (gm = dId/dVg)
gm = gradient(Id) ./ gradient(Vg);

% Ensure positive values for Id and gm for log log plots to not have the
% negative error
idx = Id > 0 & gm > 0;

% Apply the valid idx range to Id and gm
Vg = Vg(idx);
Id = Id(idx);
gm = gm(idx);



% Plot Transconductance vs. Bias Current
figure;
loglog(Id, gm, 'o-');
xlabel('Bias Current (A)');
ylabel('Transconductance (S)');
title('Transconductance vs. Bias Current log-log');
grid on;

% Plot Transconductance Efficiency vs. Bias Current
figure;
loglog(Id, (gm ./ Id), 'o-');
xlabel('ln(Id) (A)');
ylabel('ln(gm/Id) (S)');
title('Transconductance Efficiency vs. Bias Current');
grid on;

%% Finding Ith
% Define the x-limits for the two regions
x1_start = 2.49977e-10;
x1_end = 1.96282e-07;
x2_start = 5.1611e-06;
x2_end = 0.0013855;

% Extract the indices for the regions based on x-limits
ind1 = find(Id >= x1_start & Id <= x1_end);
ind2 = find(Id >= x2_start & Id <= x2_end);

% Fit a linear model to the first region (x1_start to x1_end)
p1 = polyfit(log(Id(ind1)), log(gm(ind1) ./ Id(ind1)), 1); % Linear fit
% Fit a linear model to the second region (x2_start to x2_end)
p2 = polyfit(log(Id(ind2)), log(gm(ind2) ./ Id(ind2)), 1); % Linear fit

% Calculate the intersection point of the two lines
intersection_log_id = (p2(2) - p1(2)) / (p1(1) - p2(1)); % x value of intersection
threshold_current = exp(intersection_log_id); % Threshold current (Id)

% Display the threshold current
disp(['Threshold Current (I_th) = ', num2str(threshold_current)]);

% Create extended ranges for the lines to meet
% Extend the range of Id to plot the lines beyond the data points
extended_Id1 = linspace(min(Id(ind1)), threshold_current, 100); % Range for first line (extend to threshold)
extended_Id2 = linspace(threshold_current, max(Id(ind2)), 100); % Range for second line (extend from threshold)

% Plot Transconductance Efficiency vs. Bias Current on a separate plot
figure;
loglog(Id, (gm ./ Id), 'o-', 'DisplayName', 'Data points');
hold on;

% Plot the two fitted lines on the graph (extended)
plot(extended_Id1, exp(polyval(p1, log(extended_Id1))), 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit (Region 1)');
plot(extended_Id2, exp(polyval(p2, log(extended_Id2))), 'g-', 'LineWidth', 2, 'DisplayName', 'Linear Fit (Region 2)');

% Plot the intersection point
plot(threshold_current, exp(polyval(p1, log(threshold_current))), 'bo', 'MarkerSize', 10, 'DisplayName', ['Threshold Current (I_{th}) = ', num2str(threshold_current, '%.2e')]);

% Labels and title
xlabel('ln(Id) (A)');
ylabel('ln(gm/Id) (S)');
title('Transconductance Efficiency vs. Bias Current with Extended Lines');
grid on;
legend('show');
hold off;

%% Finding Vth
% When Ith = 2.111e-06
% Looking at the gate voltage where the threshold current is will give us
% the threshold voltage
Vth = 0.745;
disp('Threshold Voltage = ');
disp(Vth);

%% Plotting the Graphs with requirements

%% Plot Transconductance vs. Bias Current
figure;
loglog(Id, gm, 'o-'); 
hold on;

% Indicate threshold current
xline(threshold_current, '--k', 'Ith');

% Annotate regions
text(1e-9, max(gm) * (0.02), 'Subthreshold', 'Color', 'blue', 'FontSize', 12);
text(1e-5, max(gm) * (0.02), 'Above Threshold', 'Color', 'blue', 'FontSize', 12);

xlabel('Bias Current (A)');
ylabel('Transconductance (S)');
title('Transconductance vs. Bias Current (log-log)');
grid on;
hold off;

%% Plot Transconductance Efficiency vs. Bias Current
figure;
loglog(Id, (gm ./ Id), 'o-');
xlabel('ln(Id) (A)');
ylabel('ln(gm/Id) (S)');
title('Transconductance Efficiency vs. Bias Current');
grid on;
hold on;

% Add vertical line at Ith to indicate threshold boundary
xline(threshold_current, '--k', 'Ith', 'LabelVerticalAlignment', 'middle');
xline(threshold_current./10, '--b'); % Denotes the two decades of moderate inversion around Ith
xline(threshold_current.*10, '--b');

% Annotate inversion regions at the top center of the plot
% Get the current x-axis limits to center the annotations
xLimits = xlim;

% Annotate inversion regions
text(min(Id) * 2, max(gm ./ Id) * 0.6, 'Weak Inversion', 'Color', 'b', 'FontSize', 12);  % Move right
text(threshold_current - 12e-07, max(gm ./ Id) * 0.3, 'Moderate Inversion:', 'Color', 'k', 'FontSize', 12);  % Move left
text(threshold_current - 12e-07, max(gm ./ Id) * 0.25, 'Two Decades around Ith', 'Color', 'k', 'FontSize', 12);  % Move left
text(threshold_current + 500e-07 , max(gm ./ Id) * 0.6, 'Strong Inversion', 'Color', 'r', 'FontSize', 12);


hold off;
%%

disp('For the transconductance plot the subthreshold region slope should scale with the bias and the theoretical value should be 1/2. The above threshold slope scales faster and should be around 1 theoretically.');
disp('For the transconductance efficiency plot the slope in the weak inversion region should theoretically be around 0 and the slope of the strong inversion region should be around -1/2.')





%% Slope plots

%Bias vs Gm

% Convert data to log scale
log_Id = log10(Id);
log_gm = log10(gm);

% Ensure valid index selection
idx1 = find(Id >= 1 & Id <= threshold_current);
idx2 = find(Id >= threshold_current & Id <= 3e-3);

% Check if sufficient points exist
if length(idx1) < 2
    warning('Not enough points for subthreshold fit. Expanding range.');
    idx1 = find(Id >= min(Id) & Id <= threshold_current); % Expand range
end

if length(idx2) < 2
    warning('Not enough points for above-threshold fit. Expanding range.');
    idx2 = find(Id >= threshold_current & Id <= max(Id)); % Expand range
end

% Perform curve fitting only if enough points exist
if length(idx1) >= 2
    p1 = polyfit(log_Id(idx1), log_gm(idx1), 1);
    slope_subthreshold = p1(1);
else
    slope_subthreshold = NaN;
    warning('Subthreshold region fit failed due to insufficient points.');
end

if length(idx2) >= 2
    p2 = polyfit(log_Id(idx2), log_gm(idx2), 1);
    slope_above_threshold = p2(1);
else
    slope_above_threshold = NaN;
    warning('Above-threshold region fit failed due to insufficient points.');
end

% Display slopes
fprintf('Slope in subthreshold region: %.3f\n', slope_subthreshold);
fprintf('Slope in above-threshold region: %.3f\n', slope_above_threshold);

% Plot trend lines on log-log plot
figure;
loglog(Id, gm, 'o-'); 
hold on;
xline(threshold_current, '--k', 'Ith');

% Add fitted lines if valid
if ~isnan(slope_subthreshold)
    fit_gm1 = 10.^(p1(1) * log_Id + p1(2));
    loglog(Id(idx1), fit_gm1(idx1), 'r-', 'LineWidth', 2); % Subthreshold fit
end

if ~isnan(slope_above_threshold)
    fit_gm2 = 10.^(p2(1) * log_Id + p2(2));
    loglog(Id(idx2), fit_gm2(idx2), 'g-', 'LineWidth', 2); % Above-threshold fit
end

% Annotate regions
if ~isnan(slope_subthreshold)
    text(1e-9, max(gm) * (0.02), sprintf('Subthreshold \n(Slope: %.2f)', slope_subthreshold), 'Color', 'blue', 'FontSize', 12);
end

if ~isnan(slope_above_threshold)
    text(1e-5, max(gm) * (0.02), sprintf('Above Threshold \n(Slope: %.2f)', slope_above_threshold), 'Color', 'blue', 'FontSize', 12);
end

xlabel('Bias Current (A)');
ylabel('Transconductance (S)');
title('Transconductance vs. Bias Current (log-log)');
grid on;
hold off;



%Efficiecny vs Id
%% Finding Ith
% Define the x-limits for the two regions
x1_start = 2.49977e-10;
x1_end = 1.96282e-07;
x2_start = 5.1611e-06;
x2_end = 0.0013855;

% Extract the indices for the regions based on x-limits
ind1 = find(Id >= x1_start & Id <= x1_end);
ind2 = find(Id >= x2_start & Id <= x2_end);

% Fit a linear model to the first region (x1_start to x1_end)
p1 = polyfit(log(Id(ind1)), log(gm(ind1) ./ Id(ind1)), 1); % Linear fit
% Fit a linear model to the second region (x2_start to x2_end)
p2 = polyfit(log(Id(ind2)), log(gm(ind2) ./ Id(ind2)), 1); % Linear fit

% Extract the slopes from the fits
slope1 = p1(1); % Slope of the first fit
slope2 = p2(1); % Slope of the second fit

% Display the slopes
disp(['Slope of the first region: ', num2str(slope1)]);
disp(['Slope of the second region: ', num2str(slope2)]);

% Calculate the intersection point of the two lines
intersection_log_id = (p2(2) - p1(2)) / (p1(1) - p2(1)); % x value of intersection
threshold_current = exp(intersection_log_id); % Threshold current (Id)

% Display the threshold current
disp(['Threshold Current (I_th) = ', num2str(threshold_current)]);

% Create extended ranges for the lines to meet
% Extend the range of Id to plot the lines beyond the data points
extended_Id1 = linspace(min(Id(ind1)), threshold_current, 100); % Range for first line (extend to threshold)
extended_Id2 = linspace(threshold_current, max(Id(ind2)), 100); % Range for second line (extend from threshold)

% Plot Transconductance Efficiency vs. Bias Current on a separate plot
figure;
loglog(Id, (gm ./ Id), 'o-', 'DisplayName', 'Data points');
hold on;

% Plot the two fitted lines on the graph (extended) with slope values in the legend
plot(extended_Id1, exp(polyval(p1, log(extended_Id1))), 'r-', 'LineWidth', 2, 'DisplayName', ['Linear Fit (Region 1), Slope = ', num2str(slope1, '%.2f')]);
plot(extended_Id2, exp(polyval(p2, log(extended_Id2))), 'g-', 'LineWidth', 2, 'DisplayName', ['Linear Fit (Region 2), Slope = ', num2str(slope2, '%.2f')]);

% Exclude the intersection point from the plot
% Commenting out the intersection plot line
% plot(threshold_current, exp(polyval(p1, log(threshold_current))), 'bo', 'MarkerSize', 10, 'DisplayName', ['Threshold Current (I_{th}) = ', num2str(threshold_current, '%.2e')]);

% Labels and title
xlabel('ln(Id) (A)');
ylabel('ln(gm/Id) (S)');
title('Transconductance Efficiency vs. Bias Current with Extended Lines');
grid on;
legend('show');
hold off;


