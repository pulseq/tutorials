% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% basic parameters
seq=mr.Sequence(sys);           % Create a new sequence object
fov=256e-3; Nx=128; Ny=Nx;      % Define FOV and resolution
alpha=90;                       % flip angle
sliceThickness=3e-3;            % slice
%TR=21e-3;                      % ignore TR, go as fast as possible
%TE=60e-3;                      % ignore TE, go as fast as possible
% TODO: run this EPI sequence with and without PE, to calibrate the delay
% TODO: change MaxGrad/MaxSlew/roDuration to see what happends to the
% stimulation (e.g. 80/200/500)

% more in-depth parameters
pe_enable=1;                    % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
rfDuration=3e-3;
roDuration=640e-6;              % not all values are possible, watch out for the checkTiming output

% Create alpha-degree slice selection pulse and corresponding gradients 
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

% define the output trigger to play out with every slice excitation
trig=mr.makeDigitalOutputPulse('ext1','duration', 100e-6,'delay', rf.delay+mr.calcRfCenter(rf)); % possible channels: 'osc0','osc1','ext1'

% Define other gradients and ADC events
deltak=1/fov; % Pulseq default units for k-space are inverse meters
kWidth = Nx*deltak;

% start with the blip
blip_dur = ceil(2*sqrt(deltak/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2; % we round-up the duration to 2x the gradient raster time
gyBlip = mr.makeTrapezoid('y',sys,'Area',-deltak,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end, each equal to a half of blip_dur
% the area between the blips should be equal to kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slew rate and then scale down the amlitude to fix the area 
extra_area=blip_dur/2*blip_dur/2*sys.maxSlew;
gx = mr.makeTrapezoid('x',sys,'Area',kWidth+extra_area,'duration',roDuration+blip_dur);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx=mr.scaleGrad(gx,kWidth/actual_area);
%adc = mr.makeAdc(Nx,'Duration',roDurtion,'Delay',blip_dur/2,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys); % if no 'Duration' is provided shortest possible duration will be used
gyPre = mr.makeTrapezoid('y','Area',(Ny/2-1)*deltak,'system',sys);

% calculate ADC - it is quite trickly
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be quite different from Nx and
% readoutTime/Nx, respectively. 
adcDwellNyquist=deltak/gx.amplitude; % dwell time on the top of the plato
% round-down dwell time to sys.adcRasterTime (100 ns)
adcDwell=floor(adcDwellNyquist/sys.adcRasterTime)*sys.adcRasterTime;
adcSamples=floor(roDuration/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell);
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % Pulseq (andsiemens) define the samples to happen in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)/sys.rfRasterTime)*sys.rfRasterTime; 
          % above we adjust the delay to align the trajectory with the gradient.
          % We have to aligh the delay to seq.rfRasterTime (1us) 
          % this rounding actually makes the sampling points on odd and even readouts
          % to appear misalligned. However, on the real hardware this misalignment is
          % much stronger anyways due to the grdient delays

% finish the blip gradient calculation
% split the blip into two halves and produce a combined synthetic gradient
gyBlip_parts = mr.splitGradientAt(gyBlip, blip_dur/2, sys);
[gyBlip_up,gyBlip_down,~]=mr.align('right',gyBlip_parts(1),'left',gyBlip_parts(2),gx);
% now for inner echos create a special gy gradient, that will ramp down to 0, stay at 0 for a while and ramp up again
gyBlip_down_up=mr.addGradients({gyBlip_down, gyBlip_up}, sys);

% pe_enable support
gyBlip_up=mr.scaleGrad(gyBlip_up,pe_enable);
gyBlip_down=mr.scaleGrad(gyBlip_down,pe_enable);
gyBlip_down_up=mr.scaleGrad(gyBlip_down_up,pe_enable);
gyPre=mr.scaleGrad(gyPre,pe_enable);

% gradient spoiling
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys); % 4 cycles over the slice thickness

% skip timing (TE/TR calculation), we'll accept the shortest TE/TR

% define sequence blocks
seq.addBlock(gz,trig);
seq.addBlock(mr.align('left',gyPre,gzReph,'right',gxPre));
for i=1:Ny % loop over phase encodes
    if i==1
        seq.addBlock(gx,gyBlip_up,adc); % Read the first line of k-space with a single half-blip at the end
    elseif i==Ny
        seq.addBlock(gx,gyBlip_down,adc); % Read the last line of k-space with a single half-blip at the beginning
    else
        seq.addBlock(gx,gyBlip_down_up,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
    end 
    gx = mr.scaleGrad(gx,-1);   % Reverse polarity of read gradient),'Duration',mr.calcDuration(gxPre),'system',sys);
end
seq.addBlock(gzSpoil);

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
seq.setDefinition('Name', 'epi');
seq.write('epi.seq')       % Write to pulseq file
%seq.install('siemens');

%% plot sequence and k-space diagrams

%seq.plot();
seq.plot('timeDisp','us','showBlocks',1); %detailed view

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');

%% PNS calc

[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
%[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_GC98SQ.asc'); % aera-xq
%[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2298_2250V_793A_SC72CD_EGA.asc'); % TERRA-XR

if (pns_ok)
    fprintf('PNS check passed successfully\n');
else
    fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);
