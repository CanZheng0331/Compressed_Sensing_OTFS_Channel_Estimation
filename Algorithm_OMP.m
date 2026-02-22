function hEst = Algorithm_OMP(Psi, yPilot, params)
    % Input:
    %   Psi    - Sensing matrix 
    %   yPilot - Received pilot signal
    %   params - Parameter structure
    % Output:
    %   hEst   - Estimated channel response
    
    % --- OMP Core Computational Logic ---
    % FUNCTION solves yPilot = Psi hEst, takes the input parameters yPilot, Psi, params.taps where
    % yPilot is the output field, Psi is the dataset field and params.taps is the sparsity. It
    % return the solution of hEst.
    
    xbeg = zeros(size(Psi, 2), 1); % Initialize a vector of all zeros with the same number of columns as the dictionary matrix
    support = [];                % Initialize the support set
    temp = yPilot;               % Initialize the residual vector to the original observed signal y
    count = 1;                   % Counter
    
    % Iterate (params.taps) times, finding the atom that best matches the residual in each iteration
    while count < params.taps + 1
        ST = abs(Psi' * temp);     % Calculate the inner product (i.e., correlation) between all columns of the Psi and the current residual
        [~, b] = max(ST);        % Find the index b of the column with the highest relevance (ignoring the maximum value itself)
        support = [support b];   % Add the index to the support set
        
        % Least square
        xfinal = Psi(:, support) \ yPilot; 
        
        % Update residual
        temp = yPilot - Psi(:, support) * xfinal;
        count = count + 1;
    end
    
    % Reconstruct the final sparse signal
    hEst = xbeg;
    t = support';
    hEst(t) = xfinal;            % Fill the non-zero estimate into the corresponding index position
end