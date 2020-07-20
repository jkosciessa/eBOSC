function [detected] = eBOSC_episode_sparsefreq(cfg, detected, TFR)

    freqWidth = (2/cfg.eBOSC.wavenumber)*cfg.eBOSC.F;
    lowFreq = cfg.eBOSC.F-(freqWidth/2);
    highFreq = cfg.eBOSC.F+(freqWidth/2);

    for indF = 1:numel(cfg.eBOSC.F)
        if ~isempty(find(cfg.eBOSC.F<=lowFreq(indF), 1, 'last'))
            fmat(indF,1) = find(cfg.eBOSC.F<=lowFreq(indF), 1, 'last')+1; % first freq falling into range
        else fmat(indF,1) = 1;
        end
        if ~isempty(find(cfg.eBOSC.F>=highFreq(indF), 1, 'first'))
            fmat(indF,3) = find(cfg.eBOSC.F>=highFreq(indF),1, 'first')-1; % last freq falling into range
        else fmat(indF,3) = numel(cfg.eBOSC.F);
        end
    end; fmat(:,2) = (1:numel(cfg.eBOSC.F))'; clear indF;
    range = diff(fmat,[],2);
    range = max(range,[],1);

    % initialize variables
    % append search space (i.e. range at both ends. first index refers to lower range
    tmp_B    = [zeros(range(1,1),size(TFR,2)); TFR.*detected; zeros(range(1,2),size(TFR,2))];
    %tmp_B    = TFR.*detected;
    detected = zeros(size(detected));

    for f = range(1,1)+1:size(tmp_B,1)-range(1,2) % loop across frequencies. note that indexing respects the appended segments

        r = 1;
        % encode detected positions at current frequency
        tmp_det(r,:) = tmp_B(f,:) ~= 0;
        % encode detected positions within LOWER range where current freq
        % point is larger
        for r = 1:range(1,1)
            tmp_det(r,:) = tmp_B(f,:) ~= 0 & tmp_B(f,:) >= tmp_B(f-r,:);
        end
        % encode detected positions within HIGHER range where current freq
        % point is larger
        for r = 1:range(1,2)
            r2 = range(1,1)+r;
            tmp_det(r2,:) = tmp_B(f,:) ~= 0 & tmp_B(f,:) >= tmp_B(f+r,:);
        end; clear r r2;

        detected(f,:) = mean(tmp_det,1) == 1;

    end; clear f;

    % remove padded zeros
    detected = detected(range(1,1)+1:size(tmp_B,1)-range(1,2),:);

end
