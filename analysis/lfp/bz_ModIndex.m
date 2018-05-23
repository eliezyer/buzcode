function [comod] = bz_ModIndex(lfp,phaserange,amprange,flagPlot)
%This function calculates the comodulogram of phase-amplitude between
%phaserange (lower frequencies) to amplitude range (higher frequencies)
%
%%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%    phaserange     [min:steps:max] array of frequencies to filter for the 
%                   phase signal
%    amprange       [min:stepsmax] array of frequencies range for wavelets 
%                   for the power signal
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       not implemented yet.... 'nfreqs' for amplitude signal
%                               'ncyc' for wavelet parameters
%                               'interval' interval on which to calculate
%                               other filter parameters
%                               'phasesignal' option for other signal to be
%                                               used as phase information
%    =========================================================================
%
%OUTPUT
%   comod               Modulation Index matrix between phase frequency and
%                       amplitude
%
%Dependencies
%   bz_Filter
%   bz_WaveSpec
%
% Implemented by Eliezyer de Oliveira, 2018
% Last Update: 23/05/2018
%



nfreqs = length(amprange);

%% Filter LFP for the phase
for bnd = 1:length(phaserange)-1
    filtered_phase(bnd,:) = bz_Filter(lfp,'passband',phaserange(bnd:bnd+1),'filter','fir1');
end

%% Wavelet Transform LFP in intervals
wavespec_amp = bz_WaveSpec(lfp,'frange',[amprange(1) amprange(end)],'nfreqs',nfreqs);
%[ampfreqs,~,spec_int]
%spec_int = cellfun(@(X) abs(X),spec_int,'UniformOutput',false);
wavespec_amp.data = log10(abs(wavespec_amp.data)); %log scaled for normality
wavespec_amp.mean = mean(wavespec_amp.data,1);


%% Bin phase and power
numbins = 50;
phasebins = linspace(-pi,pi,numbins+1);
phasecenters = phasebins(1:end-1)+(phasebins(2)-phasebins(1));

for idx = 1:length(filtered_phase)
    [phasedist,~,phaseall] = histcounts(filtered_phase(idx).phase,phasebins);
    
    phaseAmp = zeros(numbins,nfreqs);
    for bb = 1:numbins
        phaseAmp(bb,:) = mean(wavespec_amp.data(phaseall==bb,:),1)./wavespec_amp.mean;
    end
    
    phaseAmp = phaseAmp./sum(phaseAmp,1);
    comod(:,idx) = sum(phaseAmp.*log(phaseAmp./(ones(numbins,size(phaseAmp,2))/numbins)))/log(numbins);
end


 
ampfreqs = wavespec_amp.freqs;
%% Plot
if flagPlot
figure
    imagesc(phaserange,log2(ampfreqs),comod);
    colormap jet
    hold on
    xlabel('Frequency phase');
    ylabel('Frequency amplitude')
    LogScale('y',2)
    colorbar
    axis xy

end

