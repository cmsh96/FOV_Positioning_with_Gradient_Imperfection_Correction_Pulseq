%% Estimation of actual gradient using GIRF

% load GIRF
load('...\GIRF\spher\x\Results29993\SphericalHarmonics_X_1.mat'); % load the girf signal form the path
girfX_FT = GIRF_FT;
girfX_f = freqFull;
clear GIRF_FT freqFull roTime

load('...\GIRF\spher\y\Results29993\SphericalHarmonics_Y_1.mat');
girfY_FT = GIRF_FT;
girfY_f = freqFull;
clear GIRF_FT freqFull roTime

load('...\GIRF\spher\z\Results29993\SphericalHarmonics_Z_1.mat');
girfZ_FT = GIRF_FT;
girfZ_f = freqFull;
clear GIRF_FT freqFull roTime

% load spiral sequence
spiral_seq=mr.Sequence();
spiral_seq.read('...\sp64_iso.seq'); % load the multi-shot (64) spiral sequence
g_ideal = spiral_seq.waveforms_and_times;
gx_ideal = g_ideal{1}; gy_ideal = g_ideal{2}; gz_ideal = g_ideal{3};
clear g_ideal spiral_seq

% prepare g_ideal: add zero at the beginning (or end if necessary)
t_end = max( [max(gx_ideal(1,end)), max(gy_ideal(1,end)), max(gz_ideal(1,end))] );
if gx_ideal(1,1)~=0
    gx_ideal = [[0;0], gx_ideal]; end
if gx_ideal(1,end)~=t_end
    gx_ideal = [gx_ideal, [0;0]]; end

if gy_ideal(1,1)~=0
    gy_ideal = [[0;0], gy_ideal]; end
if gy_ideal(1,end)~=t_end
    gy_ideal = [gy_ideal, [0;0]]; end

if gz_ideal(1,1)~=0
    gz_ideal = [[0;0], gz_ideal]; end
if gz_ideal(1,end)~=t_end
    gz_ideal = [gz_ideal, [0;0]]; end

dt_grad = 5e-6;  % gradient raster time 5 microseconds

% prepare gxyz_cmd (uniform time)
time_vec = 0:dt_grad:t_end;
gx_ideal_ready = zeros(2,length(time_vec)); gx_ideal_ready(1,:) = time_vec;
gy_ideal_ready = zeros(2,length(time_vec)); gy_ideal_ready(1,:) = time_vec;
gz_ideal_ready = zeros(2,length(time_vec)); gz_ideal_ready(1,:) = time_vec;
gx_ideal_ready(2,:) = interp1(gx_ideal(1,:), gx_ideal(2,:), time_vec, 'linear');
gy_ideal_ready(2,:) = interp1(gy_ideal(1,:), gy_ideal(2,:), time_vec, 'linear');
gz_ideal_ready(2,:) = interp1(gz_ideal(1,:), gz_ideal(2,:), time_vec, 'linear');

G_actual = apply_girf_diag(gx_ideal_ready, gy_ideal_ready, gz_ideal_ready, dt_grad, girfX_FT, girfX_f, girfY_FT, girfY_f, girfZ_FT, girfZ_f);

% Plot ideal vs predicted actual gradients
figure;
subplot(311);
plot(G_actual.t, gx_ideal_ready(2,1:end-1), 'b', G_actual.t, real(G_actual.gx_act), 'r');
xlabel('Time [us]'); ylabel('G_x [Hz/m]');
legend('Ideal','Predicted actual'); grid on;

subplot(312);
plot(G_actual.t, gy_ideal_ready(2,1:end-1), 'b', G_actual.t, real(G_actual.gy_act), 'r');
xlabel('Time [us]'); ylabel('G_y [Hz/m]');
legend('Ideal','Predicted actual'); grid on;

subplot(313);
plot(G_actual.t, gz_ideal_ready(2,1:end-1), 'b', G_actual.t, real(G_actual.gz_act), 'r');
xlabel('Time [us]'); ylabel('G_z [Hz/m]');
legend('Ideal','Predicted actual'); grid on;

