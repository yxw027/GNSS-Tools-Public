function obs = measMask(obs,PARAMS)


if PARAMS.measUse.dfOnly
   obs.range.obs(obs.range.sig < 100) = 0;
end

% remove doppler measurements that don't have associated pseudoranges.
% actually also only keep L1 doppler measurements
obs.doppler.obs(2:end,:) = 0;
obs.doppler.obs(:,~any(obs.range.obs)) = 0;

% obs.range.obs(obs.range.ind == 2) = 0;

% obs.range.obs(obs.range.sig > 100) = 0;

obs.range.obs([3 5],:) = 0;
% obs.range.obs([4 6],:) = 0;


if PARAMS.measUse.L1_THRESH > 0 || PARAMS.measUse.L2_THRESH > 0
    % SNR mask
    
    
    
end


end