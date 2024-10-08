% this is an experimental spiral sequence

fov=256e-3; Nx=96; Ny=Nx;  % Define FOV and resolution
sliceThickness=3e-3;       % slice thinckness
Nslices=1;
Oversampling=2;            % oversampling along the readout (trajectory) dimension; I would say it needs to be at least 2
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
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',sys,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);

% define k-space parameters
deltak=1/fov;
kRadius = round(Nx/2);
kSamples=round(2*pi*kRadius)*Oversampling;

% calculate a raw Archimedian spiral trajectory
clear ka;
ka(kRadius*kSamples+1)=1i; % init as complex
for c=0:kRadius*kSamples
    r=deltak*c/kSamples;
    a=mod(c,kSamples)*2*pi/kSamples;
    ka(c+1)=r*exp(1i*a);
end
ka=[real(ka); imag(ka)];
% calculate gradients and slew rates
[ga, sa]=mr.traj2grad(ka);

% limit analysis
safety_magrin=0.94; % we need that  otherwise we just about violate the slew rate due to the rounding errors
dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(sys.maxGrad*safety_magrin)*sys.gradRasterTime;
dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(sys.maxSlew*safety_magrin))*sys.gradRasterTime;
dt_smooth=max([dt_gabs;dt_sabs]);

% apply the lower limit not to lose the trajectory detail
dt_min=4*sys.gradRasterTime/kSamples; % we want at least 4 points per revolution
dt_smooth0=dt_smooth;
dt_smooth(dt_smooth<dt_min)=dt_min;
t_smooth=[0 cumsum(dt_smooth,2)];
kopt_smooth=interp1(t_smooth, ka', (0:floor(t_smooth(end)/sys.gradRasterTime))*sys.gradRasterTime)';

% analyze what we've got
fprintf('duration orig %d us\n', round(1e6*sys.gradRasterTime*length(ka)));
fprintf('duration opt %d us\n', round(1e6*sys.gradRasterTime*length(kopt_smooth)));

[gos, sos]=mr.traj2grad(kopt_smooth);

figure;plot([gos;abs(gos(1,:)+1i*gos(2,:))]');title('gradient with abs(grad|slew) constraint')
figure;plot([sos;abs(sos(1,:)+1i*sos(2,:))]');title('slew rate with abs(grad|slew) constraint')

% Define gradients and ADC events
spiral_grad_shape=gos;
% Create 90 degree slice selection pulse and gradient
[rf, gz] = mr.makeSincPulse(pi/2,'system',sys,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
gzReph = mr.makeTrapezoid('z',sys,'Area',-gz.area/2);

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
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',mr.calcDuration(gzReph));%lims.adcDeadTime);

% extend spiral_grad_shape by repeating the last sample
% this is needed to accomodate for the ADC tuning delay
spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)];

% readout grad 
gx = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:),'Delay',mr.calcDuration(gzReph));
gy = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:),'Delay',mr.calcDuration(gzReph));

% spoilers
gz_spoil=mr.makeTrapezoid('z',sys,'Area',deltak*Nx*4);
gx_spoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(1,end),0]); %todo: make a really good spoiler
gy_spoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(2,end),0]); %todo: make a really good spoiler

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
    seq.addBlock(mr.rotate('z',phi,gzReph,gx,gy,adc));
    seq.addBlock(mr.rotate('z',phi,gx_spoil,gy_spoil,gz_spoil));
    %seq.addBlock(gx_combined,gy_combined,gz_combined,adc);
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
