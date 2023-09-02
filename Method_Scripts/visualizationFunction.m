function visualizationFunction(centralFreqs, freqs, positiveSpectrum, BW_factor, method)
    % Visualization function for DB-VMD and VMD in progress or in the end
    K = length(centralFreqs);
    clf
    h = nan(K+1,1);
    for k=1:K
        hold on
        h(k) = plot(freqs, abs(1 ./ (1 + BW_factor * (freqs - centralFreqs(k)).^2)), "black");
    end
    hold on
    h(K+1) = plot(freqs, abs(positiveSpectrum) / max(abs(positiveSpectrum)), 'b', 'LineWidth', 0.1);
    set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    xlim([0,1/2])
    ylim([0,1.01])
    legend(h(K:K+1),"Wiener filters", "Spectrum")
    xlabel("Normalized positive frequencies")
    ylabel("Normalized Magnitude")
    title("Visualization of " + method)
    pause(0.1)
end
