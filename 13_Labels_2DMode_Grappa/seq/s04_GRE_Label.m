%% 2D GRE sequence with data labels
% It demonstrates the use of the LABEL extension to produce raw data
% reconstructable by the integrated ICE reconstruction on the
% scanner and also make the raw data convertable to ISMRM-RD e.g. to be 
% used with Gadgetron

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

%% prepare label events
% Set PAT scan flags (correspond to Mdh.setPATRefScan, Mdh.setPATRefAndImaScan)
lblSetRefScan = mr.makeLabel('SET','REF', true) ;
lblSetRefAndImaScan = mr.makeLabel('SET','IMA', true) ;
lblResetRefScan = mr.makeLabel('SET','REF', false) ;
lblResetRefAndImaScan = mr.makeLabel('SET','IMA', false) ;
lblIncSlc = mr.makeLabel('INC','SLC', 1) ;

%% accelerate computations
% preregister constant objects to accelerate computations
% this is not necessary, but accelerates the sequence creation
gz.id = seq.registerGradEvent(gz) ;
gzReph.id = seq.registerGradEvent(gzReph) ;
gxPre=mr.align('left',delayTE,'right',gxPre); % update gxPre by aligning it to the right of the TE filling delay; this needs to be done prior to the pre-registering
gxPre.id = seq.registerGradEvent(gxPre) ;
gx_add.id = seq.registerGradEvent(gx_add) ;
gzSpoil.id = seq.registerGradEvent(gzSpoil) ;
% the phase of the RF object will change, therefore we only per-register the shapes
[~, rf.shapeIDs] = seq.registerRfEvent(rf) ; 
% rf.id = seq.registerRfEvent(rf) ; % NO GO EXAMPLE!!!

% precalculate and register phase encoding gradients
clear gyPre gyReph; % needed in this tutorial because these arrays are incompatible between sequences
for iY = Ny:-1:1
    gtmp = mr.scaleGrad(gy, phaseAreas(iY)/gy.area);
    gtmp.id = seq.registerGradEvent(gtmp) ;
    gyPre(iY) = gtmp;
    gtmp = mr.scaleGrad(gyr, phaseAreas(iY)/gy.area);
    gtmp.id = seq.registerGradEvent(gtmp) ;
    gyReph(iY) = gtmp;
end

% also pre-register non-changing labeling events
lblSetRefScan.id = seq.registerLabelEvent(lblSetRefScan);
lblSetRefAndImaScan.id = seq.registerLabelEvent(lblSetRefAndImaScan);
lblResetRefScan.id = seq.registerLabelEvent(lblResetRefScan);
lblResetRefAndImaScan.id = seq.registerLabelEvent(lblResetRefAndImaScan);
lblIncSlc.id = seq.registerLabelEvent(lblIncSlc);

%% Loop over phase encodes and define sequence blocks
% initialize RF spoiling
rf_phase=0;
rf_inc=0;

tic ;

for s=1:Nslices
    % reset the line counter prior to the (new) PE loop
    seq.addBlock(mr.makeLabel('SET', 'LIN', 0)); 
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
        seq.addBlock(delayTE, gzReph, gyPre(cPE), gxPre);
        seq.addBlock(gx_add, adc, gzSpoil, gyReph(cPE) ) ;
        seq.addBlock(delayTR, mr.makeLabel('INC','LIN', 1)) ; % we do not pre-register the line increment event now because with the iPAT implementation it will change inside the loop
    end
    seq.addBlock(lblIncSlc);
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
seq.setDefinition('Name', 'gre_s04') ;

% switch the 2D mode
seq.setDefinition('SlicePositions', slicePositions) ; 
seq.setDefinition('SliceThickness', sliceThickness) ;
seq.setDefinition('SliceGap', sliceGap) ;

seq.write('gre_s04.seq')       % Write to pulseq file

%% evaluate label settings more specifically
adc_lbl = seq.evalLabels('evolution','adc') ;
figure ; 
hold on ;
plot(adc_lbl.LIN) ;
plot(adc_lbl.SLC) ;

legend('LIN','SLC') ;
title('evolution of labels/counters') ;

%% plot sequence and k-space diagrams
seq.plot('timeRange', [0 5]*TR) ;

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b') ; % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

return; % skip the seq.testReport() call, it is the same as for s01 & s02...

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  

rep = seq.testReport;
fprintf([rep{:}]);

