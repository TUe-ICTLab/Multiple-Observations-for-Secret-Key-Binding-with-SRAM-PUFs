%% calculate mutual information
% calculate the basic statistics of the SRAM PUF given parameter settings
% lambda1 and lambda2
% this script is used for producing Fig. 2 and Fig. 4 in the Multiple
% Observations paper
clear all
close all

% settings
rng(210);              % Set RNG state for repeatability
lambda1 = 0.51;
lambda2 = 0; % should be zero for no bias, i.e., Pr(X=1)=.5
t_max = 20; % max number of enrollments

[~,pdftheta,~] = generate_cdf_theta(lambda1,lambda2);
theta = pdftheta(2,:);
p_theta = pdftheta(1,:);
% H(X)
p1 = [theta;1-theta];
pp = sum(p1.*p_theta,2); % average over theta distr
HX = entropyLK(pp);

% I(X;\Theta) (max achievable rate)
HX_T = sum(entropyLK(p1).*p_theta,2); % H(X|\Theta)
IX_T = HX - HX_T; % I(X;\Theta)

% now start a loop to calculate the joint entropies
pold = p1;
HXX = zeros(1,t_max+1);HXX(1) = HX;
for i=1:t_max % kan sowieso tot 20
    pnew = calcnewdistr(theta,pold);
    pp = sum(pnew.*p_theta,2); % average over theta distr
    HXX(i+1) = entropyLK(pp);
    pold = pnew;
end
IX_Y = HX(1)- diff(HXX);
 figure;plot(IX_T*ones(size(IX_Y)),'--');
hold on;plot(IX_Y);title('Mutual information after \tau enrollments');

 xlabel('\tau enrollments');
 ylabel('I(X_1 .. X_\tau;Z)');grid on;
 legend('I(X;\Theta)','Location','southeast');
 text(t_max-1,IX_T,sprintf('bound: %0.4f',IX_T),'Color','r') 
 figure;
 plot(HXX-HX(1));hold on;plot(log2(1:t_max))
 legend('H(X_1 ... X_t|Y)','log2(M)');

function y = calcnewdistr(theta,theta_old)
% given the old theta.^k * (1-theta).^n-k generate new ones
theta = repmat(theta,size(theta_old,1),1);
y = repmat(theta_old,2,1).*[theta;(1-theta)];
end





function [cdftheta,pdftheta,STATS] = generate_cdf_theta(lambda1,lambda2)
    % generate cdf and pdf of one-probabilities theta
    % 1st row is cdf and pdf
    % 2nd row is values
    % statistics = [P1, Pe_avg, Pe_dom]
    % Pe_avg = average error probability
    % Pe_dom = avg. error prob. w.r.t. dominant value
    % P1 = average one-probability
    nstepstheta = 201; % number of steps in the cdf
    theta = linspace(0,1,nstepstheta);
    cdftheta = normcdf(lambda1*norminv(theta,0,1)-lambda2,0,1);
    pdftheta = [diff(cdftheta);theta(1:end-1)+diff(theta(1:2))/2];
    pdftheta(1,:) = pdftheta(1,:)/sum(pdftheta(1,:)); % normalize
    cdftheta = [cdftheta;theta];
    
    % verify by plotting the pdf
    if 1
        figure;plot(pdftheta(2,:),pdftheta(1,:));
        xlabel('theta');ylabel('p(theta)');
    end
    % average error probability between two observations of a cell
    Pe_avg = sum(2.*pdftheta(2,:).*(1-pdftheta(2,:)).*pdftheta(1,:));
    % average error probability w.r.t. dominant value
    Pe_dom = sum(min([pdftheta(2,:);1-pdftheta(2,:)],[],1).*pdftheta(1,:));
    P1 = sum(pdftheta(2,:).*pdftheta(1,:)); % Pr(X=1), bias
    STATS = [P1,Pe_avg,Pe_dom];
end

function entr = entropyLK(p)
% ENTROPY: calculate the binary entropy of distribution given by p
% p should be column vector
entr = -p.*log2(p);
entr(p==0|p==1) = 0;
entr = sum(entr,1);
end