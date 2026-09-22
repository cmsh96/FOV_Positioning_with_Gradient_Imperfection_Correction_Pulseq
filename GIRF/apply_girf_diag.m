% calculate the actual gradient waveform based on the nominal gradient
% waveforms and the measured GIRF signal
function G_actual = apply_girf_diag(gx_ideal, gy_ideal, gz_ideal, dt_grad, girfX_FT, girfX_f, girfY_FT, girfY_f, girfZ_FT, girfZ_f)
    gx_ideal = gx_ideal(2,1:end-1);
    gy_ideal = gy_ideal(2,1:end-1);
    gz_ideal = gz_ideal(2,1:end-1);
    N = length(gx_ideal);
    GX_ideal = fftshift(fft(gx_ideal));
    GY_ideal = fftshift(fft(gy_ideal));
    GZ_ideal = fftshift(fft(gz_ideal));
    f_grad = ((0:N-1)-N/2) / (N*dt_grad);
    % f_grad = ((0:N-1)-N/2+1/2) / (N*dt_grad);
    
    girfX_FT = girfX_FT(2,:); 
    girfX_FT = [0, girfX_FT, 0]; girfX_f = [min(f_grad), girfX_f, max(f_grad)];
    % girfX_FT = [girfX_FT(1), girfX_FT, girfX_FT(end)]; girfX_f = [min(f_grad), girfX_f, max(f_grad)];
    
    girfY_FT = girfY_FT(3,:); 
    girfY_FT = [0, girfY_FT, 0]; girfY_f = [min(f_grad), girfY_f, max(f_grad)];
    % girfY_FT = [girfY_FT(1), girfY_FT, girfY_FT(end)]; girfY_f = [min(f_grad), girfY_f, max(f_grad)];

    girfZ_FT = girfZ_FT(4,:); 
    girfZ_FT = [0, girfZ_FT, 0]; girfZ_f = [min(f_grad), girfZ_f, max(f_grad)];
    % girfZ_FT = [girfZ_FT(1), girfZ_FT, girfZ_FT(end)]; girfZ_f = [min(f_grad), girfZ_f, max(f_grad)];

    Hxx = interp_complex(girfX_f(:), girfX_FT, f_grad);
    Hyy = interp_complex(girfY_f(:), girfY_FT, f_grad);
    Hzz = interp_complex(girfZ_f(:), girfZ_FT, f_grad);
    GX_act = Hxx .* GX_ideal;
    GY_act = Hyy .* GY_ideal;
    GZ_act = Hzz .* GZ_ideal;
    gx_act = ifft(fftshift(GX_act));
    gy_act = ifft(fftshift(GY_act));
    gz_act = ifft(fftshift(GZ_act));
    G_actual.gx_act = gx_act; G_actual.gy_act = gy_act; G_actual.gz_act = gz_act;
    % % Plot nominal vs predicted actual k-space trajectory
    % kx_act = cumsum(gx_act) * dt_grad;
    % ky_act = cumsum(gy_act) * dt_grad;
    % kz_act = cumsum(gz_act) * dt_grad;
    % G_actual.t = (0:N-1) * dt_grad;
    % G_actual.kx_act = kx_act; G_actual.ky_act = ky_act; G_actual.kz_act = kz_act;
end
function yq = interp_complex(x, y, xq)
    yr = interp1(x, real(y), xq, 'linear');
    yi = interp1(x, imag(y), xq, 'linear');
    yq = yr + 1i*yi;
end


% 
% % kopt=interp1(t_smooth, ka', t_grad)';
% 
% [t_smooth_u, ia] = unique(t_smooth, 'stable');
% ka_u = ka(:, ia);
% kopt = interp1(t_smooth_u, ka_u', t_grad)';