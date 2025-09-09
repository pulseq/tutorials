% this is an experimental spiral sequence

fov=256e-3; Nx=96; Ny=Nx;  % Define FOV and resolution
sliceThickness=3e-3;       % slice thinckness
Nslices=1;
adcOversamplingFactor=2;            % oversampling along the readout (trajectory) dimension; I would say it needs to be at least 2
phi=pi/2;                  % orientation of the readout e.g. for interleaving

% Set system limits
sys = mr.opts('MaxGrad',22,'GradUnit','mT/m',...
    'MaxSlew',160,'SlewUnit','T/m/s',...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6, 'adcDeadTime', 10e-6,...
    'adcSamplesLimit', 8192);  % adcSamplesLimit is important on many Siemens platforms

seq=mr.Sequence(sys);      % Create a new sequence object
warning('OFF', 'mr:restoreShape'); % restore shape is not compatible with spirals and will throw a warning from each plot() or calcKspace() call

% Create fat-sat pulse 
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
B0=2.89; % 1.5 2.89 3.0
sat_ppm=-3.35;
sat_freq=sat_ppm*1e-6*B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqPPM',sat_ppm,'use','saturation');
rf_fs.phasePPM=-2*pi*rf_fs.freqPPM*rf_fs.center; % compensate for the frequency-offset induced phase    

gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create 90 degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(pi/2,'system',sys,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'use','excitation');

% define k-space parameters
deltak=1/fov;
kRadius = round(Nx/2) ; % the radius of the spiral in k-space, determined by the user-defiled resolution; it is also the number of windings of the spiral.
kSamples = round(2*pi*kRadius)*adcOversamplingFactor ; % the number of samples of the outest cycle of the spiral

% calculate a raw Archimedian spiral trajectory
% the single-shot Archimedian spiral is a clasical spiral trajectory, with the k-space radius
% |K| propotional to its Azimuth angle (phi): K(t) = phi(t)/(2*pi*FOV).
clear ka;
ka(kRadius*kSamples+1) = 1i ; % initialize as complex
for c = 0:kRadius*kSamples % total number of k-space points
    r = deltak*c/kSamples ;
    a = mod(c,kSamples)*2*pi/kSamples;
    ka(c+1) = r*exp(1i*a) ; % in the polar coordinate, k(t) = r(t) * exp(1i*phi(t))
end
ka=[real(ka); imag(ka)]; % convert complex numbers to k-vectors
% calculate gradients and slew rates
[ga, sa] = mr.traj2grad(ka,'RasterTime',sys.gradRasterTime,'firstGradStepHalfRaster',true,'conservativeSlewEstimate',true);

% limit analysis
safety_margin = 0.94 ; % we need that, otherwise we just about violate the slew rate due to the rounding errors
dt_gabs = abs(ga(1,:) + 1i*ga(2,:))/(sys.maxGrad*safety_margin)*sys.gradRasterTime ; % dt in case of decreased g
dt_sabs = sqrt(abs(sa(1,:)+1i*sa(2,:))/(sys.maxSlew*safety_margin))*sys.gradRasterTime ; % dt in case of decreased slew

dt_opt=max([dt_gabs;dt_sabs]) ; % select the larger dt in case of decreased g and s

% apply the lower limit not to lose the trajectory detail
dt_min = 4*sys.gradRasterTime/kSamples; % we want at least 4 points per revolution (cycle)
dt_opt0 = dt_opt ;
dt_opt(dt_opt<dt_min) = dt_min ;

figure;plot([dt_opt0; dt_opt]');title('combined time stepping');

t_smooth = [0 cumsum(dt_opt,2)] ;

t_end = t_smooth(end)-0.5*sys.gradRasterTime ;
t_grad = [0 (0.5+(0:floor(t_end/sys.gradRasterTime)))*sys.gradRasterTime] ;

kopt=interp1(t_smooth, ka', t_grad)';

% analyze what we've got
fprintf('duration orig %d us\n', round(1e6*sys.gradRasterTime*length(ka)));
fprintf('duration opt %d us\n', round(1e6*sys.gradRasterTime*length(kopt)));

[gos, sos]=mr.traj2grad(kopt,'RasterTime',sys.gradRasterTime,'firstGradStepHalfRaster',false);

figure;plot([gos;abs(gos(1,:)+1i*gos(2,:))]');
hold on; yline(sys.maxGrad,'--'); title('gradient with the abs constraint');

figure;plot([sos;abs(sos(1,:)+1i*sos(2,:))]');
hold on; yline(sys.maxSlew,'--'); title('slew rate with the abs constraint')

% Define gradients and ADC events
spiral_grad_shape=gos;

% calculate ADC
% round-down dwell time to 10 ns
adcTime = sys.gradRasterTime*size(spiral_grad_shape,2);
% actually it is trickier than that: the (Siemens) interpreter sequence 
% per default will try to split the trajectory into segments with the number of samples <8192
% and every of these segments will have to have duration aligned to the
% gradient raster time

adcSamplesDesired=kRadius*kSamples; 
adcDwell=round(adcTime/adcSamplesDesired/sys.adcRasterTime)*sys.adcRasterTime; 
adcSamplesDesired=ceil(adcTime/adcDwell);
[adcSegments,adcSamplesPerSegment]=mr.calcAdcSeg(adcSamplesDesired,adcDwell,sys); 

adcSamples=adcSegments*adcSamplesPerSegment;

% we would like to sample the point k=0 with the first ADC sample (i.e. at 
% t=adcDwell/2), so we advance the ADC and round the delay to the RF raster time
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',round((mr.calcDuration(gzReph)-adcDwell/2)/sys.rfRasterTime)*sys.rfRasterTime);%lims.adcDeadTime);

% extend spiral_grad_shape by repeating the last sample
% this is needed to accomodate for the ADC tuning delay
spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)];

% readout grad 
gx = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:),'Delay',mr.calcDuration(gzReph),'first',0,'last', spiral_grad_shape(1,end),'system',sys);
gy = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:),'Delay',mr.calcDuration(gzReph),'first',0,'last', spiral_grad_shape(2,end),'system',sys);

% spoilers
gz_spoil=mr.makeTrapezoid('z',sys,'Area',deltak*Nx*4);
gx_spoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(1,end),0],'system',sys); %todo: make a really good spoiler
gy_spoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(2,end),0],'system',sys); %todo: make a really good spoiler

% because of the ADC alignment requirements the sampling window possibly
% extends past the end of the trajectory (these points will have to be
% discarded in the reconstruction, which is no problem). However, the
% ramp-down parts and the Z-spoiler now have to be added to the readout
% block otherwise there will be a gap inbetween
% gz_spoil.delay=mr.calcDuration(gx);
% gx_spoil.delay=gz_spoil.delay;
% gy_spoil.delay=gz_spoil.delay;
% gx_combined=mr.addGradients([gx,gx_spoil], lims);
% gy_combined=mr.addGradients([gy,gy_spoil], lims);
% gz_combined=mr.addGradients([gzReph,gz_spoil], lims);
 
% Define sequence blocks
for s=1:Nslices
    seq.addBlock(rf_fs,gz_fs); % fat-sat    
    rf.freqOffset=gz.amplitude*sliceThickness*(s-1-(Nslices-1)/2);
    seq.addBlock(rf,gz);
    seq.addBlock(mr.rotate('z',phi,gzReph,gx,gy,adc,'system',sys));
    seq.addBlock(mr.rotate('z',phi,gx_spoil,gy_spoil,gz_spoil,'system',sys));
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

%
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'spiral');
seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % this is important for making the sequence run automatically on siemens scanners without further parameter tweaking

seq.write('spiral.seq');   % Output sequence for scanner

% the sequence is ready, so let's see what we got 
seq.plot();             % Plot sequence waveforms

%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); title('2D k-space');

% seq.install('siemens');
return

%% calculate PNS and check whether we are still within limits

[pns,tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
