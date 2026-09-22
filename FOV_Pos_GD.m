%% Run the original sequence (2d radial)
fov = 256e-3;
sliceThickness = 3e-3; 

% run the radial gre or spiral sequence and save the seq file:
% rad_gre.m     % radial gre sequence
% spiral64.m    % 64-shot spiral sequence

fpo = make_fovpos_phaseCorr_rfPhaseCorr_rgre(o, delta_r, arbitrary_delay, fov, sliceThickness);

%% slow down your sequence by arbitrary delay (1e-5)
arbitrary_delay = 1e-5;
[o, ok1] = make_slow_seq(seq, arbitrary_delay, fov, sliceThickness);%, adcSamplesPerSegment);
o.setDefinition('Name', 'o');
o.write('o.seq');

% FOV Positioning Parameters
Dx = arbitrary_delay;  Dy = Dx;  Dz = Dx;
delta_r = [16, 24, 0] * 1e-3; 
T = mr.TransformFOV('translation', delta_r); % rad gre
% T = mr.TransformFOV_sp('translation', delta_r); % spiral

%% First Path: o > do > fdo
[do, ok2] = make_delayed_seq(o, arbitrary_delay, fov, sliceThickness);%, adcSamplesPerSegment);
fdo = T.applyToSeq(do);
do.setDefinition('Name', 'do');
do.write('do.seq');
fdo.setDefinition('Name', 'fdo');
fdo.write('fdo.seq'); 

%% Second Path: o > fo > dfo
fo = T.applyToSeq(o);
[dfo, ok3] = make_delayed_seq(fo, arbitrary_delay, fov, sliceThickness);%, adcSamplesPerSegment);
fo.setDefinition('Name', 'fo');
fo.write('fo.seq');
dfo.setDefinition('Name', 'dfo');
dfo.write('dfo.seq');

%% Third Path: o > (fp)o > d(fp)o
fpo = make_fovpos_phaseCorr(o, delta_r, arbitrary_delay, fov, sliceThickness);%, adcSamplesPerSegment);
% fpo = make_fovpos_phaseCorr_rfPhaseCorr(o, delta_r, arbitrary_delay, fov, sliceThickness);
% fpo = make_fovpos_phaseCorr_rfPhaseCorr_rgre(o, delta_r, arbitrary_delay, fov, sliceThickness);
[dfpo, ok4] = make_delayed_seq(fpo, arbitrary_delay, fov, sliceThickness);%, adcSamplesPerSegment);
fpo.setDefinition('Name', 'fpo');
fpo.write('fpo.seq');
dfpo.setDefinition('Name', 'dfpo');
dfpo.write('dfpo.seq');


%% helper functions

% we shift rf and adc back by arbitrary_delay
% function [delayed_seq, ok] = make_delayed_seq(seq, arbitrary_delay, fov, sliceThickness, adcSamplesPerSegment) % for spiral only
function [delayed_seq, ok] = make_delayed_seq(seq, arbitrary_delay, fov, sliceThickness) % gre rad 
    delayed_seq = mr.Sequence(seq.sys);
    for iB=1:length(seq.blockEvents)
        B = seq.getBlock(iB);
        if ~isempty(B.rf) || ~isempty(B.adc)
            if ~isempty(B.rf)
                B.rf.delay = B.rf.delay - arbitrary_delay;
            end
            if ~isempty(B.adc)
                B.adc.delay = B.adc.delay - arbitrary_delay;
            end
        end
        % B=rmfield(B,'gexternal');
        delayed_seq.addBlock(B);
    end
    % figure; seq.plot('showBlocks', 1, 'timeDisp','us', 'stacked', 1);
    % figure; delayed_seq.plot('showBlocks', 1, 'timeDisp','us', 'stacked', 1);
    %
    [ok, error_report]=delayed_seq.checkTiming;
    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end
    delayed_seq.setDefinition('FOV', [fov fov sliceThickness]);
    % delayed_seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % for spiral only
end







% this code reads an example sequence and add a delay at the begining of the blocks that have rf or adc
% function [slow_seq, ok] = make_slow_seq(seq, arbitrary_delay, fov, sliceThickness, adcSamplesPerSegment) % for spiral    
function [slow_seq, ok] = make_slow_seq(seq, arbitrary_delay, fov, sliceThickness) % for gre radial
    slow_seq = mr.Sequence(seq.sys);
    % original_seq.read('rgre_iso.seq');
    % arbitrary_delay = 1e-5;
    for iB=1:length(seq.blockEvents)
        B = seq.getBlock(iB);
        if ~isempty(B.rf) || ~isempty(B.adc)
            if ~isempty(B.rf)
                B.rf.delay = B.rf.delay + arbitrary_delay;
            end
            if ~isempty(B.adc)
                B.adc.delay = B.adc.delay + arbitrary_delay;
            end
            if ~isempty(B.gx)
                B.gx.delay = B.gx.delay + arbitrary_delay;
            end
            if ~isempty(B.gy)
                B.gy.delay = B.gy.delay + arbitrary_delay;
            end
            if ~isempty(B.gz)
                B.gz.delay = B.gz.delay + arbitrary_delay;
            end
            B.blockDuration = B.blockDuration + arbitrary_delay;
        end
        % B=rmfield(B,'gexternal');
        slow_seq.addBlock(B);
    end
    % figure; seq.plot('showBlocks', 1, 'timeDisp','us', 'stacked', 1); 
    % figure; slow_seq.plot('showBlocks', 1, 'timeDisp','us', 'stacked', 1);
    %
    [ok, error_report]=slow_seq.checkTiming;
    if (ok)
        fprintf('Timing check passed successfully\n');
    else
        fprintf('Timing check failed! Error listing follows:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end
    slow_seq.setDefinition('FOV', [fov fov sliceThickness]);
    % slow_seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % for spiral only
end







% function gdDelayCorrected_Transformed_seq = make_fovpos_phaseCorr(seq, delta_r, arbitrary_delay, fov, sliceThickness, adcSamplesPerSegment)
function gdDelayCorrected_Transformed_seq = make_fovpos_phaseCorr(seq, delta_r, arbitrary_delay, fov, sliceThickness)
gdDelayCorrected_Transformed_seq = mr.Sequence(seq.sys); 
% Dx = arbitrary_delay;  Dy = Dx;  Dz = Dx;
Dx = arbitrary_delay(1);  Dy = arbitrary_delay(2);  Dz = arbitrary_delay(3);

% T = mr.TransformFOV_f2('translation', delta_r);
T = mr.TransformFOV_sp('translation', delta_r);

transformed_seq = T.applyToSeq(seq);
grad_waveforms = transformed_seq.waveforms_and_times;
gx_ideal = grad_waveforms{1,1};  [gx_ideal_pp, gx_ideal_int_pp] = make_pp(gx_ideal);
gy_ideal = grad_waveforms{1,2};  [gy_ideal_pp, gy_ideal_int_pp] = make_pp(gy_ideal);
gz_ideal = grad_waveforms{1,3};  [gz_ideal_pp, gz_ideal_int_pp] = make_pp(gz_ideal);

% shift the grads based on their delay
gx_actual_pp = gx_ideal_pp;  gx_actual_pp.breaks = gx_actual_pp.breaks + Dx;  gx_actual_int_pp = fnint(gx_actual_pp);
gy_actual_pp = gy_ideal_pp;  gy_actual_pp.breaks = gy_actual_pp.breaks + Dy;  gy_actual_int_pp = fnint(gy_actual_pp);
gz_actual_pp = gz_ideal_pp;  gz_actual_pp.breaks = gz_actual_pp.breaks + Dz;  gz_actual_int_pp = fnint(gz_actual_pp);

% extract rf and adc samples of the transformed sequence and calculate the
% phase difference based on the actual and ideal grads and add it to the
% sequence
%
for iB=1:length(transformed_seq.blockEvents)
    % B = transformed_seq.getBlock(iB);
    B = seq.getBlock(iB);
    if ~isempty(B.rf)
        rf = B.rf;
        if iB==1
            t_rf = rf.t + rf.delay;
        else
            t_rf = rf.t + rf.delay + sum(transformed_seq.blockDurations(1:iB-1));
        end
        rf_phase_diff = 2*pi*(  (ppval(gx_actual_int_pp, t_rf) - ppval(gx_ideal_int_pp, t_rf)) * delta_r(1) + ...
                                (ppval(gy_actual_int_pp, t_rf) - ppval(gy_ideal_int_pp, t_rf)) * delta_r(2) + ...
                                (ppval(gz_actual_int_pp, t_rf) - ppval(gz_ideal_int_pp, t_rf)) * delta_r(3)  );
        rf.signal = rf.signal .* exp(1i * rf_phase_diff);
        B.rf = rf;
    end
    if ~isempty(B.adc)
        adc = B.adc;
        if iB==1
            t_adc = adc.dwell*(0.5:adc.numSamples-0.5) + adc.delay;
        else
            t_adc = adc.dwell*(0.5:adc.numSamples-0.5) + adc.delay + sum(transformed_seq.blockDurations(1:iB-1));
        end
        adc_phase_diff = 2*pi*(  (ppval(gx_actual_int_pp, t_adc) - ppval(gx_ideal_int_pp, t_adc)) * delta_r(1) + ...
                                 (ppval(gy_actual_int_pp, t_adc) - ppval(gy_ideal_int_pp, t_adc)) * delta_r(2) + ...
                                 (ppval(gz_actual_int_pp, t_adc) - ppval(gz_ideal_int_pp, t_adc)) * delta_r(3)  );
        if isempty(adc.phaseModulation)
            adc.phaseModulation = adc_phase_diff(:);     % radians
        else
            adc.phaseModulation = adc.phaseModulation + adc_phase_diff(:);
        end
        B.adc = adc;
    end 
    % B = rmfield(B,'gexternal');
    gdDelayCorrected_Transformed_seq.addBlock(B);
end
gdDelayCorrected_Transformed_seq.setDefinition('FOV', [fov fov sliceThickness]);
% gdDelayCorrected_Transformed_seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % for spiral only
end