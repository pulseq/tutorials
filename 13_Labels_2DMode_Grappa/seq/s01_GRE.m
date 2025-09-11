%% 2D GRE sequence with some time optimisation 
% This sequence has a reasonably optimal timing and serves as a starting
% point for this tutorial. Please study it to understand single steps.
% Amongtst other things, it demonstrates how combined readout and spoiler
% gradinets van be realised. It also demonstrates the use of the mr.align()
% functionality.

%% set system limits and parameters
sys = mr.opts('MaxGrad', 22, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6) ;

fov = 256e-3 ;                     % Define FOV 
Nx = 256 ; Ny = 256 ;              % resolution / matrix size
alpha = 10 ;                       % flip angle
Nslices = 3;                       % number of slices
sliceThickness = 3e-3 ;            % slice thickness & gap
sliceGap = 1e-3 ;
TR = 20e-3 ;                       % repetition time TR
TE = 8e-3 ;                        % echo time TE
roDuration = 5.12e-3 ;             % ADC duration
rfSpoilingInc = 117;               % RF spoiling increment

seq = mr.Sequence(sys) ;           % Create a new sequence object

%% Create alpha-degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'use','excitation','system',sys);

%% Define other gradients and ADC events
deltak = 1/fov ;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys) ;
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys) ;
gxPre = mr.makeTrapezoid('x','Area',-gx.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gx.riseTime),'system',sys) ;
phaseAreas = ((0:Ny-1) - Ny/2) * deltak ;
gy = mr.makeTrapezoid('y','Area',max(abs(phaseAreas)),'Duration',mr.calcDuration(gxPre,gzReph),'system',sys) ;

%% gradient spoiling in X and Z directions, refocusing in Y
% we cut the RO gradient into two parts for the optimal spoiler timing
[gx1, ~] = mr.splitGradientAt(gx, gx.riseTime + gx.flatTime) ;
% gradient spoiling
gxSpoil = mr.makeExtendedTrapezoidArea(gx.channel, gx.amplitude, 0, 2*Nx*deltak, sys) ;
gxSpoil.delay = mr.calcDuration(gx1) ;
gx_add = mr.addGradients({gx1, gxSpoil},'system',sys) ;

gzSpoil = mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys) ;
gzSpoil.delay = max(mr.calcDuration(gx1), gzSpoil.delay) ;

gyr=mr.scaleGrad(gy,-1);
gyr.delay = mr.calcDuration(gx1);

%% slice positions
slicePositions=(sliceThickness+sliceGap)*((0:(Nslices-1)) - (Nslices-1)/2);
%slicePositions=slicePositions([1:2:Nslices 2:2:Nslices]); % reorder slices for an interleaved acquisition (optional)

%% Calculate timing
delayTE = ceil((TE - gz.fallTime - gz.flatTime + rf.center ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime) * seq.gradRasterTime ;
delayTR = round((TR - mr.calcDuration(gz) - delayTE ...
    - mr.calcDuration(gx_add, gzSpoil, gyr))/seq.gradRasterTime)*seq.gradRasterTime ;
assert(delayTE >= mr.calcDuration(gxPre,gy,gzReph) ) ;
assert(delayTR >= 0 ) ;

%% precalculate phase encoding gradients
clear gyPre gyReph; % needed in this tutorial because these arrays are incompatible between sequences
for iY = Ny:-1:1
    gtmp = mr.scaleGrad(gy, phaseAreas(iY)/gy.area);
    gyPre(iY) = gtmp;
    gtmp = mr.scaleGrad(gyr, phaseAreas(iY)/gy.area);
    gyReph(iY) = gtmp;
end

%% Loop over phase encodes and define sequence blocks
% initialize RF spoiling
rf_phase=0;
rf_inc=0;

tic ;
    
for s=1:Nslices
    for cPE=1:Ny
        % configure RF for the slice position
        rf.freqOffset   = gz.amplitude*slicePositions(s);
        rf.phaseOffset  = -2*pi*rf.freqOffset*rf.center ... % compensate for the slice-offset induced phase
                        + rf_phase ; % and RF spoiling (vary RF phase pseudo-randomly)
        adc.phaseOffset = rf_phase ; % make ADC phase follow the RF phase
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        %
        seq.addBlock(rf, gz) ;
        % we put all objects into a single align function to remain compatible with seq.addBlock(); 
        % actually ...(delayTE,'right',gxPre) would be sufficient because all other objects are 
        % already on the left by default; we also align the phase-encoding gradient to the center 
        % just to demonstrate capabilities of the mr.align() function, even if it makes no sense from the MR physics perspective
        seq.addBlock(mr.align('left', delayTE, gzReph, 'center', gyPre(cPE), 'right', gxPre)); 
        seq.addBlock(gx_add, adc, gzSpoil, gyReph(cPE) ) ;
        seq.addBlock(delayTR) ; 
    end
end
toc ;

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
seq.setDefinition('FOV', [fov fov (max(slicePositions)-min(slicePositions)+sliceThickness)]) ;
seq.setDefinition('Name', 'gre') ;

seq.write('gre.seq')       % Write to pulseq file

%% plot sequence and k-space diagrams
seq.plot('timeRange', [0 5]*TR) ;

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);

