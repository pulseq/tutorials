% set system limits (slew rate 130 and max_grad 30 work on Prisma)
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 10e-6,'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);           % Create a new sequence object
fov=240e-3; Nx=240;             % Define FOV and resolution
alpha=5;                        % flip angle
sliceThickness=3e-3;            % slice
% TE & TR are as short as possible derived from the above parameters and
% the system specs above
Nr=30;                         % number of radial spokes
Ndummy=1;                      % number of dummy scans
delta= pi / Nr;                 % angular increment; try golden angle pi*(3-5^0.5) or 0.5 of it
% more in-depth parameters
ro_dur=1200e-6;                 % RO duration
ro_os=2;                        % readout oversampling
ro_spoil=0.5;                   % additional k-max excursion for RO spoiling
rf_dur = 600e-6;
sl_spoil=2;                     % spoil area compared to the slice thickness



% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment

% Create alpha-degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rf_dur,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',2,'system',sys);
gzReph.delay=mr.calcDuration(gz);
gzComb=mr.addGradients({gz, gzReph}, 'system', sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','Amplitude',Nx*deltak/ro_dur,'FlatTime',ceil(ro_dur/sys.gradRasterTime)*sys.gradRasterTime,'system',sys);
adc = mr.makeAdc(Nx*ro_os,'Duration',ro_dur,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.amplitude*(ro_dur/Nx/ro_os*(Nx*ro_os/2-0.5)+0.5*gx.riseTime),'system',sys); % 0.5 is necessary to account for the Siemens sampling in the center of the dwell periods
% start gxPre at least right after the RF pulse and when possible end it at the same time as the end of the slice refocusing gradient 
[gxPre,~,~]=mr.align('right', gxPre, 'right', gzComb, 'left', mr.makeDelay(mr.calcDuration(rf)+mr.calcDuration(gxPre)));

% gradient spoiling
if sl_spoil>0
    sp_area_needed=sl_spoil/sliceThickness-gz.area/2;
    gzSpoil=mr.makeTrapezoid('z','Area',sp_area_needed,'system',sys,'Delay',gx.riseTime+gx.flatTime);
else
    gzSpoil=[];
end

if ro_spoil>0
    ro_add_time=ceil(((gx.area/Nx*(Nx/2+1)*ro_spoil)/gx.amplitude)/sys.gradRasterTime)*sys.gradRasterTime;
    gx.flatTime=gx.flatTime+ro_add_time; % careful, areas stored in the object are now wrong
end


% we don't calculate timing but just accept what is achievable
TR=0;
TE=0;

% start the sequence
rf_phase=0;
rf_inc=0;
for i=(1-Ndummy):3
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);
    %
    phi=delta*(i-1);
    seq.addBlock(mr.rotate('z',phi,rf,gzComb,gxPre));
    if TE<=0, TE=seq.duration+adc.delay+adc.dwell*(adc.numSamples/2+0.5); end
    if i>0 
        seq.addBlock(mr.rotate('z',phi,gx,adc,gzSpoil));
    else
        seq.addBlock(mr.rotate('z',phi,gx,gzSpoil));
    end
    if TR<=0, TR=seq.duration; end
end

% check whether the timing of the sequence is compatible with the scanner
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
seq.setDefinition('Name', 'gre_rad');

seq.write('gre_rad.seq')       % Write to pulseq file
%seq.install('siemens');

%% plot sequence and k-space diagrams

%seq.plot('timeRange',[0 5]*TR);
seq.plot();

% trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing1] = seq.calculateKspacePP();
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay',[0 0 0]*1e-6); % play with anisotropic trajectory delays -- zoom in to see the trouble ;-)

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
title('k-vector components as functions of time'); xlabel('time /s'); ylabel('k-component /m^-^1');
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('2D k-space trajectory'); xlabel('k_x /m^-^1'); ylabel('k_y /m^-^1');

return;

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slew rate limits  

rep = seq.testReport; 
fprintf([rep{:}]); 

