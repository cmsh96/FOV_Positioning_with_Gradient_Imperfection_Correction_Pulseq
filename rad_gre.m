% radial gre sequence
seq=mr.Sequence();              % Create a new sequence object
fov=256e-3;    
Nx=256;             % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=3e-3;            % slice
TE=8e-3;                        % TE; give a vector here to have multiple TEs (e.g. for field mapping)
TR=20e-3;                       % only a single value for now
Nr=180;                         % number of radial spokes

Ndummy=20;                      % number of dummy scans
delta= pi / Nr;                 % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it
% delta= 2*pi / Nr; 

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment

% set system limits
sys = mr.opts('MaxGrad', 60, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% Create alpha-degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',4e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
%% option 1
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',6.4e-3/5,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2-deltak/2,'Duration',2e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',2e-3,'system',sys);

%% option 2
% deltak=1/fov;
% gx2 = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',6.4e-3/10,'system',sys);
% adc = mr.makeAdc(2*Nx,'Duration',gx2.flatTime,'Delay',gx2.riseTime,'system',sys);
% gx = gx2;
% gxPre = mr.makeTrapezoid('x','Area',-gx.area/2-deltak/2,'Duration',2e-3,'system',sys);
% gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',2e-3,'system',sys);

%%

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',0.5*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

rf_phase=0;
rf_inc=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOV Positioning Parameters
delta_r = [0, 0, 0]; 
T = mr.TransformFOV('translation',delta_r);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=(-Ndummy):Nr
    for c=1:length(TE)
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        %
        seq.addBlock(T.applyToBlock(rf,gz));

        phi=delta*(i-1);
        events = mr.rotate('z',phi,gxPre,gzReph);
        seq.addBlock(T.applyToBlock(events));
        seq.addBlock(mr.makeDelay(delayTE(c)));

        if (i>0)
            events = mr.rotate('z',phi,gx,adc);
            seq.addBlock(T.applyToBlock(events));
        else
            events = mr.rotate('z',phi,gx);
            seq.addBlock(T.applyToBlock(events));
        end
        events = mr.rotate('z',phi,gxSpoil,gzSpoil);
        seq.addBlock(T.applyToBlock(events));
        seq.addBlock(mr.makeDelay(delayTR));
    end
end


%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%%
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'seq');
% seq.write('seq.seq')       % Write to pulseq file

%seq.install('siemens');

% return;
%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

figure;plot3(ktraj(1,:),ktraj(2,:),ktraj(3,:),'b');
hold on;plot3(ktraj_adc(1,:),ktraj_adc(2,:),ktraj_adc(3,:),'r.');

figure;plot3(ktraj_adc(1,:),ktraj_adc(2,:),ktraj_adc(3,:),'r.');

