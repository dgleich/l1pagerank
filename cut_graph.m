function [A,alpha,beta] = cut_graph(G,S,alpha,beta)
n = size(G,1);

% cases for S
% S is a list of indices
% S is a logical indicator 
Sind = [];
if islogical(S) && numel(S) == n
    % Do nothing, 
    Sind = S;
elseif numel(S) < n
    Sind = false(n,1);
    Sind(S) = 1;
else
    % numel(S) == n and S is not a logical
    if max(S) > 1
        % S is an indicator
        Sind = false(n,1);
        Sind(S) = 1;
    else
        Sind = logical(S);
    end
end

Sbar = ~Sind;

d = full(sum(G,1))';
volS = sum(d(Sind));
volSbar = sum(d(Sbar));

assert(volS < volSbar);

if nargin < 3
    % use Reid's alpha, beta
    cut = volS - Sind'*G*Sind;    
    alpha = cut / volS; % see assert above
    beta = alpha*volS/volSbar;
end

ws = alpha*d.*Sind;
wsbar = beta*d.*Sbar;

A = general_cut_graph(G,ws,wsbar);