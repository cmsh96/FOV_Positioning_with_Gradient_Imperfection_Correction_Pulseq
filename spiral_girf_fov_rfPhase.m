%% girf-corrected gradients
load G_sp64.mat; % load the predicted gradient waveforms by the girf signal and the nominal gradient waveforms
% correct it by shifting 5/2=2.5us (half adc dwell time)
G_sp64_25uscor = zeros(4, size(G_sp64,2)+1);
G_sp64_25uscor(1,:) = [0, 2.5e-6+G_sp64(1,1:end-1), G_sp64(1,end)];
G_sp64_25uscor(2,:) = [0, G_sp64(2,1:end-1), (G_sp64(2,end-1)+G_sp64(2,end))/2];
G_sp64_25uscor(3,:) = [0, G_sp64(3,1:end-1), (G_sp64(3,end-1)+G_sp64(3,end))/2];
G_sp64_25uscor(4,:) = [0, G_sp64(4,1:end-1), (G_sp64(4,end-1)+G_sp64(4,end))/2];

%% original sequence
seq=mr.Sequence();
seq.read('...\sp64_iso.seq');

%% apply fov positioning: fov-shifted sequence
delta_r = [16,24,0]*1e-3;
T = mr.TransformFOV_sp('translation', delta_r);
transformed_seq = T.applyToSeq(seq);
% prepare gradients: add [0,0] at the begining of the grad waveforms and
% make pp
grad_waveforms = transformed_seq.waveforms_and_times;
gx_ideal = grad_waveforms{1,1}; 
if gx_ideal(1,1)~=0
    gx_ideal = [[0;0], gx_ideal]; end
[gx_ideal_pp, gx_ideal_int_pp] = make_pp(gx_ideal);

gy_ideal = grad_waveforms{1,2};  
if gy_ideal(1,1)~=0
    gy_ideal = [[0;0], gy_ideal]; end
[gy_ideal_pp, gy_ideal_int_pp] = make_pp(gy_ideal);

gz_ideal = grad_waveforms{1,3};  
if gz_ideal(1,1)~=0
    gz_ideal = [[0;0], gz_ideal]; end
[gz_ideal_pp, gz_ideal_int_pp] = make_pp(gz_ideal);

%% calculate pp from girf-corrected grad waveform
gx_actual = [G_sp64_25uscor(1,:);G_sp64_25uscor(2,:)];
[gx_actual_pp, gx_actual_int_pp] = make_pp(gx_actual);

gy_actual = [G_sp64_25uscor(1,:);G_sp64_25uscor(3,:)];
[gy_actual_pp, gy_actual_int_pp] = make_pp(gy_actual);

gz_actual = [G_sp64_25uscor(1,:);G_sp64_25uscor(4,:)];
[gz_actual_pp, gz_actual_int_pp] = make_pp(gz_actual);

%% calculate the residual adc phase between the nominal and actual grad waveforms and use it later in image recon
gdDelayCorrected_seq = mr.Sequence(seq.sys); 
rf_center_phase = [];
rf_counter = 0;
adc_counter = 0;
for iB=1:length(seq.blockEvents)
    iB
    % B = transformed_seq.getBlock(iB);
    B = seq.getBlock(iB);
    if ~isempty(B.rf)
        rf_counter = rf_counter + 1;
        if mod(rf_counter,2)==0
            rf = B.rf;
            if iB==1
                t_rf = rf.center + rf.delay; % t_rf = rf.t + rf.delay;
            else
                t_rf = rf.center + rf.delay + sum(transformed_seq.blockDurations(1:iB-1)); % t_rf = rf.t + rf.delay + sum(transformed_seq.blockDurations(1:iB-1));
            end
            rf_phase_diff = 2*pi*( local_frac( (ppval(gx_actual_int_pp, t_rf) - ppval(gx_ideal_int_pp, t_rf)) * delta_r(1) + ...
                                    (ppval(gy_actual_int_pp, t_rf) - ppval(gy_ideal_int_pp, t_rf)) * delta_r(2) + ...
                                    (ppval(gz_actual_int_pp, t_rf) - ppval(gz_ideal_int_pp, t_rf)) * delta_r(3) ) );
            rf_center_phase = [rf_center_phase, rf_phase_diff];
            % rf.signal = rf.signal .* exp(1i * rf_phase_diff);
            % B.rf = rf;
        end
    end

    if ~isempty(B.adc)
        adc_counter = adc_counter + 1;
        adc = B.adc;
        if iB==1
            t_adc = adc.dwell*(0.5:adc.numSamples-0.5) + adc.delay;
        else
            t_adc = adc.dwell*(0.5:adc.numSamples-0.5) + adc.delay + sum(transformed_seq.blockDurations(1:iB-1));
        end
        adc_phase_diff = 2*pi*( local_frac( (ppval(gx_actual_int_pp, t_adc) - ppval(gx_ideal_int_pp, t_adc)) * delta_r(1) + ...
                                 (ppval(gy_actual_int_pp, t_adc) - ppval(gy_ideal_int_pp, t_adc)) * delta_r(2) + ...
                                 (ppval(gz_actual_int_pp, t_adc) - ppval(gz_ideal_int_pp, t_adc)) * delta_r(3) ) ) - rf_center_phase(adc_counter);
        if isempty(adc.phaseModulation)
            adc.phaseModulation = adc_phase_diff(:);     % radians
        else
            adc.phaseModulation = adc.phaseModulation + adc_phase_diff(:); % TODO: add adc phase offset
        end
        B.adc = adc;
    end 

        gdDelayCorrected_seq.addBlock(B);
end

%% 1- add the residual adc phase from gdDelayCorrected_seq to the raw data of transformed_seq 
%  2- use the corrected trajectory of original seq after using girf and
%  then do image recon: recon_spiral.m



function out=local_frac(in)
    out = in-floor(in);
end


function Mod = accurate_mod_pp(breaks, coefs, t, shift) 
Mod = zeros(size(t));
A = [];
i_breaks = 0;
for c=1:length(t)
    % c
    index = find(breaks <= t(c), 1, 'last');
    t0 = t(c);
    area = 0;

    while index > (i_breaks + 1)
        i_breaks = i_breaks +1;
        AB = coefs(i_breaks,1) * shift * ( (breaks(i_breaks+1)-breaks(i_breaks))^2-(breaks(i_breaks)-breaks(i_breaks))^2 );
        CD = coefs(i_breaks,2) * shift * ( (breaks(i_breaks+1)-breaks(i_breaks))-(breaks(i_breaks)-breaks(i_breaks)) );
        A = [A,local_frac( local_frac(AB) + local_frac(CD) )];
    end
    
    if t0==breaks(i_breaks+1)
        area =  local_frac(sum(A));
    else
        AB_n = coefs(i_breaks+1,1) * shift * ( (t0-breaks(i_breaks+1))^2-(breaks(i_breaks+1)-breaks(i_breaks+1))^2 );
        CD_n = coefs(i_breaks+1,2) * shift * ( (t0-breaks(i_breaks+1))-(breaks(i_breaks+1)-breaks(i_breaks+1)) );
        area =  local_frac( local_frac(sum(A)) + local_frac(AB_n) + local_frac(CD_n) );
    end
    Mod(c) = area;
end
end


