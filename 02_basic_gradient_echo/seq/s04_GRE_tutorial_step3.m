% Define FOV and resolution
fov = 256e-3;
sliceThickness = 5e-3;
Nx = 256;
Ny = Nx;

% Define sequence parameters
TE = 8e-3;
TR = 22e-3;
alpha=30;

% set system limits
sys = mr.opts('MaxGrad',12,'GradUnit','mT/m',...
    'MaxSlew',100,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6);

% Create a new sequence object
seq=mr.Sequence(sys);

% Create slice selective alpha-pulse and corresponding gradients
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180, 'Duration', 4e-3,...
    'SliceThickness', sliceThickness, 'apodization', 0.5,'timeBwProduct', 4, ...
    'system' ,sys, 'use', 'excitation');

% Define other gradients and ADC events
deltak = 1/fov; % Pulseq toolbox defaults to k-space units of m^-1
gx = mr.makeTrapezoid('x', 'FlatArea', Nx*deltak, 'FlatTime', 6.4e-3);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
gxPre = mr.makeTrapezoid('x', 'Area', -gx.area/2, 'Duration', 2e-3);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;


% Calculate timing
delayTE = round((TE - mr.calcDuration(gxPre) - mr.calcDuration(gz)/2 ...
                    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = round((TR - mr.calcDuration(gxPre) - mr.calcDuration(gz) ...
                    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;

spoilArea=4*gx.area(); % 4 "looks" good
% Add spoilers in read, refocus in phase and spoiler in slice
gxPost = mr.makeTrapezoid('x', 'Area', spoilArea, 'system', sys); % we pass 'system' here to calculate shortest time gradient
gyPost = mr.makeTrapezoid('y', 'Area', -max(phaseAreas(:)), 'Duration', 2e-3);
gzPost = mr.makeTrapezoid('z', 'Area', spoilArea, 'system', sys);

delayTR = delayTR - mr.calcDuration(gxPost, gyPost, gzPost);

% Loop over phase encodes and define sequence blocks
for i=1:Ny
    % Vary RF phase quasi-randomly
    rand_phase = mod(117*(i^2 + i + 2), 360)*pi/180;
    [rf, gz] = mr.makeSincPulse(alpha*pi/180, 'Duration', 4e-3,...
                                'SliceThickness', 5e-3, ...
                                'apodization', 0.5, ...
                                'timeBwProduct', 4, ...
                                'system', sys, ...
                                'phaseOffset', rand_phase, ...
                                'use', 'excitation');
    seq.addBlock(rf, gz);
    gyPre = mr.makeTrapezoid('y', 'Area', phaseAreas(i), 'Duration', 2e-3);
    seq.addBlock(gxPre, gyPre, gzReph);
    seq.addBlock(mr.makeDelay(delayTE));
    seq.addBlock(gx, adc);
    gyPost = mr.makeTrapezoid('y', 'Area', -gyPre.area, 'Duration', 2e-3);
    % Add spoilers in read and slice and may be in phase
    seq.addBlock(gxPost, gyPost, gzPost);
    seq.addBlock(mr.makeDelay(delayTR));
end

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% export definitions
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'DEMO_gre3');

seq.write(['DEMO_gre3.seq'])       % Write to pulseq file

seq.plot('timeRange', [0 2*TR])

% do not run the rest of the script automatically
return

%% plot gradients to check for gaps and optimality of the timing
wave_data=seq.waveforms_and_times(); 
% plot the entire gradient shape
figure; plot(wave_data{1}(1,:),wave_data{1}(2,:)); xlabel('time /s'); ylabel('gradient /(Hz/m)');
hold on; plot(wave_data{2}(1,:),wave_data{2}(2,:));
plot(wave_data{3}(1,:),wave_data{3}(2,:));
legend('G_x', 'G_y', 'G_z');

%% calculate k-space trajectory 

[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
title('k-space vector components as functions of time');
legend('k_x', 'k_y', 'k_z');

figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
axis('equal'); % enforce aspect ratio for the correct trajectory display
title('k-space trajectory (k_x/k_y)');

%%

% [ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

% plot k-spaces
% time_axis=(1:(size(ktraj,2)))*sys.gradRasterTime;
% figure; plot(time_axis, ktraj'); % plot the entire k-space trajectory
% hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%% very optional step, slow but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);

%% listen to the sequence
seq.sound();
