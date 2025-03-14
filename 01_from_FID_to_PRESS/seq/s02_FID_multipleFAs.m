system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
Nx=8192;
Nrep=1;
adcDur=256e-3; 
rfDur=1000e-6; % increased to 1ms
TR=5000e-3;    % increased to 5s avoid T1 saturation
TE=30e-3; % the used coil has a very long switching time!
flip_angles=9;
% todo: change flip_angles to 18 and pi/2 to pi below

% Create non-selective pulse 
rf = mr.makeBlockPulse(pi/2,'Duration',rfDur, 'system', system, 'use', 'excitation');
    
% Define delays and ADC events
adc = mr.makeAdc(Nx,'Duration',adcDur, 'system', system, 'delay', TE-rf.shape_dur/2-system.rfRingdownTime);

delayTR=TR-mr.calcDuration(rf);
assert(delayTR>=0);

for f=1:flip_angles
    rf = mr.makeBlockPulse(pi/flip_angles*f,'Duration',rfDur, 'system', system, 'use', 'excitation');
    % Loop over repetitions and define sequence blocks
    for i=1:Nrep
        seq.addBlock(rf);
        seq.addBlock(adc,mr.makeDelay(delayTR));
    end
end

seq.setDefinition('Name', 'fid-mfa');
seq.write('fid-mfa.seq')       % Write to pulseq file
%seq.install('siemens');    % copy to scanner

seq.plot();
