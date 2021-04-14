function [n_errors,sim_settings] = sim_MOhds(n_sim,R,t_max)
% Simulate MO helper data scheme with LDPC code of rate R using n_sim Monte
% carlo simulations. Evaluate FER for a range of t_min till t_max
% enrollment observations
sim_settings.n_sim = n_sim;
sim_settings.t_max = t_max; % max number of enrollment observations
sim_settings.start_date = date;


%% prepare for the simulations
lambda1 = 0.51;
lambda2 = 0;
k = 128; % 128 bit key
n = ceil(k/R); % number of SRAM cells
R = k/n;
 
% generate the pdf of theta
[cdftheta,pdftheta,stats] = generate_cdf_theta(lambda1,lambda2);

% generate LLRTABLE
LLRTable = generate_LLR(pdftheta,t_max);

% settings LDPC (general)
% DL-SCH coding parameters
cbsInfo = nrDLSCHInfo(k,R);
 
%% start simulations
tic
n_errors = zeros(1,t_max);
for i = 1:n_sim
    % Random transport block data generation
    in = randi([0 1],k,1,'int8');
    % encode
    chIn = LDPC_encoder(cbsInfo,n,in);
    %% CHANNEL
    % generate an SRAMPUF according to the cdf
    SRAMPUF = generate_SRAMPUF(cdftheta,n);
    % first enrollment
    y = int8(rand(length(SRAMPUF),1)<SRAMPUF);
    M = xor(y,chIn);
    for n_enroll = 1:t_max
        % reconstruction
        y = int8(rand(length(SRAMPUF),1)<SRAMPUF);
        linearInd = sub2ind(size(LLRTable), ...
                repmat(n_enroll,n,1), M+1, y+1);
        LLR = LLRTable(linearInd);
        % decode
        out = LDPC_decoder(cbsInfo,k,R,LLR);
        tbErr = sum(xor(out,in),1); % number of decoder errors
        if tbErr>0
            n_errors(n_enroll) = n_errors(n_enroll)+1;
        end
        % helper data update (new enrollment)
        M = M+xor(y,chIn);
    end
    
end

if 0 % plot
    figure;
    t_enroll = (1:t_max)';
    FER = n_errors'/n_sim;
    plot(t_enroll,FER,'-*');
    set(gca, 'YScale', 'log')
    xlabel('t enrollment observations');ylabel('FER')
    title(['MO scheme, R = ',num2str(R)])
    grid on
    pause(1)
end

ldpc.k = k;
ldpc.n = n;
ldpc.R = R;
sim_settings.ldpc = ldpc;

sim_settings.hds = 'MO'; % method
sim_settings.R = R;
sim_settings.hours_simulation_time = toc/(60*60);
sim_settings.model_statistics = stats;
sim_settings.model_l1 = lambda1;
sim_settings.model_l2 = lambda2;
sim_settings.cbsInfo = cbsInfo;
end

%% LLR table construction
function LLRTable = generate_LLR(pdftheta,t_max)
    nstepstheta = size(pdftheta,2);
    % first gerenate k-ones pmf
    THETAs = repmat(reshape(pdftheta(2,:),1,1,[]),t_max+1,t_max+1+1,1);
    PTHETAs = repmat(reshape(pdftheta(1,:),1,1,[]),t_max+1,t_max+1+1,1);
    k = 0:t_max+1;t = 1:t_max+1;
    Ks = repmat(k,t_max+1,1,nstepstheta);
    Ts = repmat(t',1,t_max+1+1,nstepstheta);
    KONEs = tril(sum(THETAs.^Ks.*(1-THETAs).^(Ts-Ks).*PTHETAs,3),1);
    % then generate LLR table
    LLRTable = zeros(t_max,t_max+1,2);
    for t_enroll = 1:t_max
        n_obs = t_enroll+1; % total number of observations (enroll + 1)
        m = (0:t_enroll); % ones in helper data
        for y = 0:1 % ones in decoder obs
            LLRTable(t_enroll,1:length(m),y+1) = log(KONEs(n_obs,(m+y+1)))...
                -log(KONEs(n_obs,(t_enroll-m+y+1))); % +1 is because indexing starts at 1
        end
    end
end

%% LDPC encoder + decoder functions
function chIn = LDPC_encoder(cbsInfo,outlen,in)
    % settings
    RV = 0; % Redundancy version, 0-3
    modulation = 'BPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
    nlayers = 1;           % Number of layers, 1-4 for a transport block
    
    % Transport block CRC attachment
    tbIn = nrCRCEncode(in,cbsInfo.CRC);
    % Code block segmentation and CRC attachment
    cbsIn = nrCodeBlockSegmentLDPC(tbIn,cbsInfo.BGN);
    % LDPC encoding
    enc = nrLDPCEncode(cbsIn,cbsInfo.BGN);
    % Rate matching and code block concatenation
    chIn = nrRateMatchLDPC(enc,outlen,RV,modulation,nlayers);
end
function out = LDPC_decoder(cbsInfo,k,R,chOut)
    % settings
    RV = 0; % Redundancy version, 0-3
    modulation = 'BPSK';   % Modulation scheme, QPSK, 16QAM, 64QAM, 256QAM
    nlayers = 1;           % Number of layers, 1-4 for a transport block
    
    % Rate recovery
    raterec = nrRateRecoverLDPC(chOut,k,R,RV,modulation,nlayers);
    % LDPC decoding
    decBits = nrLDPCDecode(raterec,cbsInfo.BGN,25);
    % Code block desegmentation and CRC decoding
    [blk,~] = nrCodeBlockDesegmentLDPC(decBits,cbsInfo.BGN,k+cbsInfo.L);
    %disp(['CRC error per code-block: [' num2str(blkErr) ']'])
    % Transport block CRC decoding
    [out,~] = nrCRCDecode(blk,cbsInfo.CRC);
end

%% SRAMPUF generation
function SRAMPUF = generate_SRAMPUF(cdftheta,n)
    % generate SRAMPUF with n cells, and with one-probabilities according 
    % to given distribution
    SRAMPUF = zeros(n,1);
    RAND = rand(n,1);
    theta = cdftheta(2,:);cdftheta = cdftheta(1,:);
    dtheta = diff(theta(1:2))/2;
    for cdf_step = 1:length(cdftheta)
        select = (RAND<=cdftheta(cdf_step));
        SRAMPUF(select) = theta(cdf_step)-dtheta;
        RAND(select)= 2; %any larger number than 1
    end
%    % plot the distribution
%     if 0
%         pdftheta = [diff(cdftheta);theta(1:end-1)+diff(theta(1:2))/2];
%         pdftheta(1,:) = pdftheta(1,:)/sum(pdftheta(1,:)); % normalize
%         figure;histogram(SRAMPUF,length(theta));xlabel('theta');
%         title('verify that generated SRAM cells follow target distribution')
%         hold on;
%         plot(pdftheta(2,:),n*pdftheta(1,:));
%     end
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
%     % verify by plotting the pdf
%     if 0
%         figure;plot(pdftheta(2,:),pdftheta(1,:));
%         xlabel('theta');ylabel('p(theta)');
%         
%     end
    % average error probability between two observations of a cell
    Pe_avg = sum(2.*pdftheta(2,:).*(1-pdftheta(2,:)).*pdftheta(1,:));
    % average error probability w.r.t. dominant value
    Pe_dom = sum(min([pdftheta(2,:);1-pdftheta(2,:)],[],1).*pdftheta(1,:));
    P1 = sum(pdftheta(2,:).*pdftheta(1,:));
    STATS = [P1,Pe_avg,Pe_dom];
end

