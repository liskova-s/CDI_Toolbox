function error = error_object_domain(g_k,mask)
    % RMS
    % G...calculated FT, F...measured pattern
    num = abs(g_k.*~mask).^2;
    denom = abs(g_k).^2;
    error = sqrt(sum(num, 'all') / sum(denom, 'all'));
end