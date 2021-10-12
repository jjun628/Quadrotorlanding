function [phi_e,the_e,psi_e,d] = rotation_estimation_func(mm,mmC2)
%ROTATION_ESTIMATION_FUNC
%   This function takes two vectors (original vector, vector after the rotation) and
%   calculates rotation and translation
%   q-method is used described in 
%   https://rotations.berkeley.edu/estimating-rotations-and-translations-from-optical-targets/

    r = mm;    % mark vectors before rotation
    M = mmC2;    % mark vectors after the rotation
    
    alpha = [1 1 1 1];    % weight
    weightr = alpha*r'; weightr = weightr'./sum(alpha);
    weightM = alpha*M'; weightM = weightM'./sum(alpha);

    %% Define B,z,K
    %% 1. Define B matrix
    B = zeros(3,3);
    for j = 1:size(B,1)   %loop over rows
        for k = 1:size(B,1)   %loop over columns
            B1 = zeros(size(alpha));
            %Calculate the term in bracket, call that B1:
            for m = 1:size(alpha,2)   %to over alpha to calculate summation
                B1(m) = alpha(m)*M(j,m)*r(k,m);
            end
            B1 = sum(B1,2)./sum(alpha);
            B(j,k) = B1 - weightM(j)*weightr(k);
        end
    end

    %% 2. Define z
    z = zeros(3,size(alpha,2));
    for j = 1:size(alpha,2)   %loop over alpha to calculate summation
        z(:,j) = cross(alpha(j)*M(:,j), r(:,j));
    end
    z = sum(z,2)./sum(alpha);
    z = cross(weightM,weightr) - z;

    %% 3. Define K
    K = [B + B' - trace(B)*eye(size(B)) z];
    K = [K; z' trace(B)];

    %% Find maximum eigenvalue of K and corresponding eigenvector:
    [dum idx] = max(eig(K));
    [eigV lam] = eig(K);
    qe = eigV(:,idx);    % where q = [e1 e2 e3 e0]

    %% Reconstruct rotation matrix using Euler-Rodrigues parameters:
    R_method = (qe(4)^2 - qe(1)^2 - qe(2)^2 - qe(3)^2)*eye(3,3) + ...
                2*qe(1:3)*qe(1:3)' + ...
                2*[0 -qe(4)*qe(3) qe(4)*qe(2); qe(4)*qe(3) 0 -qe(4)*qe(1); -qe(4)*qe(2) qe(4)*qe(1) 0];
                   
    %% calculate translation vector
    d = weightM - R_method*weightr;
    
    %% calculate rotation angles
    phi_e = -atan2(R_method(3,2), R_method(3,3));
    the_e = -atan2(-R_method(3,1),sqrt(R_method(3,2)^2+R_method(3,3)^2));
    psi_e = -atan2(R_method(2,1), R_method(1,1));
end

