% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% basic parameters
seq=mr.Sequence(sys);           % Create a new sequence object
fov=256e-3; Nx=128; Ny=128;     % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=3e-3;            % slice
TR=30e-3;                       % TR, a single value
TE=20e-3;                       % only a single TE is accepted now
nSeg=5;                         % number of segments
%TODO: increase nSeq, (remember to increase TE), see how we hit the "blip limit"

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
rfDuration=3e-3;
roDuration=640e-6;              % not all values are possible, watch out for the checkTiming output

% Create alpha-degree slice selection pulse and corresponding gradients 
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov; % Pulseq default units for k-space are inverse meters
gxp = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys); % Pulseq default units for gradient amplitudes are 1/Hz
gxm=mr.scaleGrad(gxp,-1);
adc = mr.makeAdc(Nx,'Duration',gxp.flatTime,'Delay',gxp.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gxp.area/2,'system',sys); % if no 'Duration' is provided shortest possible duration will be used

% with segmentation it gets trickier
if mod(nSeg,2)==0, warning('for even number of segments additional steps are required to avoid the segment edge hitting the k-space center, expect artifacts...'); end % the code will work but the images are likely to be affected by the discontinuity at the center of k-space
phaseAreas = ((0:floor(Ny/nSeg)-1)-Ny/2)*deltak;

% calculate the blip gradient
gyBlip = mr.makeTrapezoid('y','Area',floor(Ny/nSeg)*deltak,'Delay',gxp.riseTime+gxp.flatTime,'system',sys);
if mr.calcDuration(gyBlip)-mr.calcDuration(gxp)<gyBlip.fallTime
% QC: If the peak of gyBlip is within the event duration of gxp, move the peak to the end of gxp by adding  (gyBlip.fallTime - mr.calcDuration(gyBlip)+mr.calcDuration(gxp)) to gyBlip's delay.
    gyBlip.delay=gyBlip.delay+mr.calcDuration(gxp)-mr.calcDuration(gyBlip)+gyBlip.fallTime; % adjust the delay to have the cut line before the beginning of the ramp down of the blip
end
gyBlip_parts = mr.splitGradientAt(gyBlip, mr.calcDuration(gxp), sys);
gyBlip_parts(2).delay=0; % need to reset the delay because we will use this part in a different block
% for all echoes other than the first we may need to copy the objects
gxp0=gxp;
adc0=adc;
gyBlip_part_tmp=gyBlip_parts(1);
if mr.calcDuration(gyBlip_parts(2))>gxp.riseTime
    % need to update the copy to allow for the split gradient to ramp down
    gxp.delay=mr.calcDuration(gyBlip_parts(2))-gxp.riseTime;
    gxm.delay=mr.calcDuration(gyBlip_parts(2))-gxm.riseTime;
    adc.delay=adc.delay+gxp.delay;
    gyBlip_part_tmp.delay=gyBlip_part_tmp.delay+gxp.delay;
end
% now for inner echos create a special gy gradient, that will ramp down to 0, stay at 0 for a while and ramp up again
gyBlip_down_up=mr.addGradients({gyBlip_parts(2), gyBlip_part_tmp}, sys); % QC: gyBlip_parts(2) is inserted before gxp. gyBlip_part_tmp with a long delay is inserted after gxp.
% copy for readability
gyBlip_up=gyBlip_parts(1);
gyBlip_down=gyBlip_parts(2);

% gradient spoiling
if mod(length(TE),2)==0, spSign=-1; else, spSign=1; end
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak*spSign,'system',sys);      % 2 cycles over the voxel size in X
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys); % 4 cycles over the slice thickness

% Calculate timing (need to decide on the block structure already)
delayTE=TE - ceil((gz.fallTime + gz.flatTime/2 + nSeg*mr.calcDuration(gxp0)/2 + floor((nSeg-1)/2)*gxp.delay)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTE>=mr.calcDuration(gxPre,gzReph)));
delayTR=round((TR - mr.calcDuration(gz) - delayTE ...
    - nSeg*mr.calcDuration(gxp0) - floor((nSeg-1)/2)*gxp.delay)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

% initialize the RF spoling counters 
rf_phase=0;
rf_inc=0;

% define sequence blocks
for i=1:length(phaseAreas) % loop over phase encodes
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    adc0.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);
    %
    seq.addBlock(rf,gz);
    gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
    seq.addBlock(mr.align('left', mr.makeDelay(delayTE),gyPre,gzReph,'right',gxPre));
    for s=1:nSeg % loop over segments
        if s==1
            seq.addBlock(gxp0,adc0,gyBlip_up);
        else
            if mod(s,2)==0, gx=gxm; else, gx=gxp; end
            if s~=nSeg
                seq.addBlock(gx,adc,gyBlip_down_up);
            else
                seq.addBlock(gx,adc,gyBlip_down);
            end
        end
    end
    gyPost = mr.makeTrapezoid('y','Area',-phaseAreas(i)-gyBlip.area*(nSeg-1),'Duration',mr.calcDuration(gxPre),'system',sys);
    seq.addBlock(mr.makeDelay(delayTR),gxSpoil,gyPost,gzSpoil)
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
seq.setDefinition('Name', 'seg-gre');
seq.write('seg-gre.seq')       % Write to pulseq file
%seq.install('siemens');

%% plot sequence and k-space diagrams

%seq.plot('timeRange', [0 2]*TR);
seq.plot('timeDisp','us','showBlocks',1,'timeRange',[0 2]*TR); %detailed view

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
