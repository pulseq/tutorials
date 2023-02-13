system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
adcDur=2.56e-3; 
rfDur1=3e-3;
rfDur2=1e-3;
TR=150e-3;
TE=54e-3;
spA=1000; % spoiler area in 1/m (=Hz/m*s) % MZ: need 5000 for my oil phantom

sliceThickness=3e-3;            % slice
fov=256e-3; Nx=256;             % Define FOV and resolution
Nr=256;                         % number of radial spokes
Ndummy=10;                      % number of dummy scans
delta=pi/Nr;                    % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it
starting_angle=0;               % the orientation of the first projection

% Create 90 degree slice selection pulse and gradient and even the
% refocusing gradient; we will not use it however but will subtract its
% area from the first spoiler
[rf_ex, gs, gsr] = mr.makeSincPulse(pi/2,system,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.4,'timeBwProduct',4);

% Create non-selective refocusing pulse
rf_ref = mr.makeBlockPulse(pi,'Duration',rfDur2, 'system', system, 'PhaseOffset' ,pi/2 , 'use', 'refocusing'); % needed for the proper k-space calculation
    
% calculate spoiler gradient
g_sp1=mr.makeTrapezoid('z','Area',spA+gsr.area,'system',system);
rf_ref.delay=max(mr.calcDuration(g_sp1),rf_ref.delay);

g_sp2=mr.makeTrapezoid('z','Area',spA,'system',system);

% Define delays and ADC events
deltak=1/fov;
gr = mr.makeTrapezoid('x',system,'FlatArea',Nx*deltak,'FlatTime',adcDur);
adc = mr.makeAdc(Nx,system,'Duration',adcDur,'delay',gr.riseTime);

grPredur = 5e-3; % use a fixed time to make this gradient visible on the plot
grPre = mr.makeTrapezoid('x',system,'Area',gr.area/2+deltak/2,'Duration',grPredur); 

delayTE1=TE/2-(mr.calcDuration(gs)-mr.calcRfCenter(rf_ex)-rf_ex.delay)-rf_ref.delay-mr.calcRfCenter(rf_ref)-grPredur;
delayTE2=TE/2-mr.calcDuration(rf_ref)+rf_ref.delay+mr.calcRfCenter(rf_ref)-adc.delay-adcDur/2; % this is not perfect, but -adcDur/2/Nx  will break the raster alignment
delayTR=TR-mr.calcDuration(gs)-grPredur-delayTE1-mr.calcDuration(rf_ref)-delayTE2-mr.calcDuration(gr);

assert(delayTE1>=0);
assert(delayTE2>mr.calcDuration(g_sp2));
assert(delayTR>=0);

% Loop over repetitions and define sequence blocks
for i=(1-Ndummy):Nr
    seq.addBlock(rf_ex,gs);
    seq.addBlock(mr.rotate('z',delta*(i-1),grPre)); 
    seq.addBlock(mr.makeDelay(delayTE1));  
    seq.addBlock(rf_ref,g_sp1);
    seq.addBlock(g_sp2, mr.makeDelay(delayTE2));
    if (i>0)
        seq.addBlock(mr.rotate('z',delta*(i-1),adc,gr));  
    else
        seq.addBlock(mr.rotate('z',delta*(i-1),gr));  
    end
    seq.addBlock(mr.makeDelay(delayTR));  
end

seq.plot();

% check whether the timing of the sequence is compatible with the scanner
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'se_rad');

seq.write('se_radial.seq')       % Write to pulseq file
%seq.install('siemens');    % copy to scanner

% calculate k-space but only use it to check timing
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay',[0 0 0]*1e-6); % play with anisotropic trajectory delays -- zoom in to see the trouble ;-)

if Ndummy==0
    assert(abs(t_refocusing(1)-t_excitation(1)-TE/2)<1e-6); % check that the refocusing happens at the 1/2 of TE
    assert(abs(t_adc(Nx/2)-t_excitation(1)-TE)<adc.dwell); % check that the echo happens as close as possible to the middle of the ADC elent
end

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kz-axis
title('k-vector components as functions of time'); xlabel('time /s'); ylabel('k-component /m^-^1');
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('2D k-space trajectory'); xlabel('k_x /m^-^1'); ylabel('k_y /m^-^1');
