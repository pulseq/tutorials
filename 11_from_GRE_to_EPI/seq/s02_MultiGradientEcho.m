% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% basic parameters
seq=mr.Sequence(sys);           % Create a new sequence object
fov=256e-3; Nx=256; Ny=Nx;      % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=3e-3;            % slice
TR=21e-3;                       % TR, a single value
TE=[4 9 15]*1e-3;               % give a vector here to have multiple TEs

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
rfDuration=3e-3;
roDuration=3.2e-3;              % not all values are possible, watch out for the checkTiming output

% Create alpha-degree slice selection pulse and corresponding gradients 
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov; % Pulseq default units for k-space are inverse meters
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys); % Pulseq default units for gradient amplitudes are 1/Hz
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys); % if no 'Duration' is provided shortest possible duration will be used
gxFlyBack = mr.makeTrapezoid('x','Area',-gx.area,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',sys);      % 2 cycles over the voxel size in X
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys); % 4 cycles over the slice thickness

% Calculate timing (need to decide on the block structure already)
helperT=ceil((gz.fallTime + gz.flatTime/2 + mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
for c=1:length(TE)
    delayTE(c)=TE(c) - helperT;
    helperT=helperT+delayTE(c)+mr.calcDuration(gx);
end
assert(all(delayTE(1)>=mr.calcDuration(gxPre,gzReph)));
assert(all(delayTE(2:end)>=mr.calcDuration(gxFlyBack)));
delayTR=round((TR - mr.calcDuration(gz) - sum(delayTE) ...
    - mr.calcDuration(gx)*length(TE))/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

% initialize the RF spoling counters 
rf_phase=0;
rf_inc=0;

% define sequence blocks
for i=1:Ny % loop over phase encodes
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);
    %
    seq.addBlock(rf,gz);
    gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
    for c=1:length(TE) % loop over TEs
        if (c==1)
            seq.addBlock(mr.align('left', mr.makeDelay(delayTE(c)),gyPre,gzReph,'right',gxPre)); 
        else
            seq.addBlock(mr.align('left', mr.makeDelay(delayTE(c)),'right',gxFlyBack));
        end
        seq.addBlock(gx,adc);
        % to check/debug TE calculation with seq.testReport() comment out
        % the above line and uncommend the line below this comment block; 
        % change 'c==3' statement to select the echo to test
        % if c==3, seq.addBlock(gx,adc); else, seq.addBlock(gx); end
    end
    gyPre.amplitude=-gyPre.amplitude; % better to use mr.scaleGrad(gyPre,-1);
    seq.addBlock(mr.makeDelay(delayTR),gxSpoil,gyPre,gzSpoil)
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

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'mgre');
seq.write('mgre.seq')       % Write to pulseq file
%seq.install('siemens');

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5]*TR);

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);
