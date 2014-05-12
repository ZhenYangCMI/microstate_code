function icatb_post_process_dfnc(param_file)
%% Post process dFNC
%


icatb_defaults;

%% Select dFNC file
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('title', 'dFNC Parameter File', 'typeEntity', 'file', 'typeSelection', 'single', 'filter', '*dfnc.mat');
    drawnow;
    if (isempty(param_file))
        error('dFNC parameter file is not selected');
    end
end

load(param_file);

outputDir = fileparts(param_file);

if (isempty(outputDir))
    outputDir = pwd;
end

cd (outputDir);

outputDir = pwd;

%% Initialise vars
TR = dfncInfo.TR;

num_clusters = 6;
kmeans_max_iter = 150;
dmethod = 'city';
try
    num_clusters = dfncInfo.postprocess.num_clusters;
    kmeans_max_iter = dfncInfo.postprocess.kmeans_max_iter;
    dmethod = dfncInfo.postprocess.dmethod;
catch
end


numParameters = 1;

inputText(numParameters).promptString = 'Select number of clusters to compute k-means on dynamic FNC correlations';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(num_clusters);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'num_clusters';
inputText(numParameters).enable = 'on';

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Select max iterations for computing k-means';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(kmeans_max_iter);
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'kmeans_max_iter';
inputText(numParameters).enable = 'on';

numParameters = numParameters + 1;

opts = char('City', 'sqEuclidean', 'Hamming', 'Correlation', 'Cosine');
vals = find(strcmpi(opts, dmethod) == 1);
if (isempty(vals))
    vals = 1;
end
vals = vals(1);
inputText(numParameters).promptString = 'Select distance method for computing k-means';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = opts;
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'dmethod';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = vals;

% Input dialog box
answer = icatb_inputDialog('inputtext', inputText, 'Title', 'Select Post-processing options', 'handle_visibility',  'on');

if (isempty(answer))
    error('Post-processing options are not selected');
end

drawnow;

for n = 1:length(answer)
    tmp = answer{n};
    if (isnumeric(tmp))
        tmp = num2str(tmp);
    else
        tmp = ['''', tmp, ''''];
    end
    str = [inputText(n).tag, '=', tmp, ';'];
    eval(str);
end

dfncInfo.postprocess.num_clusters = num_clusters;
dfncInfo.postprocess.kmeans_max_iter = kmeans_max_iter ;
dfncInfo.postprocess.dmethod = dmethod;

M = length(dfncInfo.outputFiles);

SP = cell(M,1);
k1_peaks = zeros(1,M);

disp('Computing spectra on FNC correlations ...');

%% Frequency analysis
for nR = 1:length(dfncInfo.outputFiles)
    
    current_file = fullfile(outputDir,  dfncInfo.outputFiles{nR});
    load(current_file, 'FNCdyn');
    
    if (nR == 1)
        Nwin = size(FNCdyn, 1);
        Fs = 1/TR;
        nfft = 2^(nextpow2(Nwin)+1); % pad to the next power of 2
        df=Fs/nfft;
        f=0:df:Fs;    % all possible frequencies
        f=f(1:nfft);
        fpass = [0.0, 1/(2*TR)];
        findx=find(f>=fpass(1) & f<=fpass(end)); % just one-half of the spectrum
        f = f(findx);
        dfncInfo.freq = f;
        FNCdynflat = zeros(M, Nwin, size(FNCdyn, 2));
    end
    
    
    FNCdynflat(nR, :, :) = FNCdyn;
    
    DEV = std(FNCdyn, [], 2);
    [xmax, imax, xmin, imin] = icatb_extrema(DEV); % find the extrema
    pIND = sort(imax);
    k1_peaks(nR) = length(pIND);
    SP{nR} = FNCdyn(pIND, :);
    
    %% Compute spectra
    FNCdyn = FNCdyn - repmat(mean(FNCdyn, 1), Nwin, 1);
    FNCdyn = FNCdyn.*repmat(icatb_hamming(Nwin), 1, size(FNCdyn, 2));
    S = fft(FNCdyn, nfft, 1)/Fs;
    S = S(findx,:);
    spectra_fnc = sqrt(S.*conj(S));
    
    icatb_save(current_file, 'spectra_fnc', '-append');
    
    %tmp_amp = squeeze(trapz(f, spectra_fnc));
    tmp_amp = squeeze(std(spectra_fnc));
    tmp_cm = squeeze(mean(spectra_fnc.*repmat(f(:), [1, size(spectra_fnc, 2)])));
    tmp_cm = tmp_cm ./ mean(spectra_fnc);
    
    if (nR == 1)
        FNCamp = zeros(size(tmp_amp));
        FNCcm = zeros(size(tmp_cm));
    end
    
    FNCamp = FNCamp + tmp_amp;
    FNCcm = FNCcm + tmp_cm;
    
    clear FNCdyn;
end

FNCamp = FNCamp / M;
FNCcm = FNCcm / M;


post_process_file = fullfile(outputDir,  [dfncInfo.prefix, '_post_process.mat']);
icatb_save(post_process_file, 'FNCamp', 'FNCcm');

clear FNCamp FNCcm;

fprintf('\n');

disp('Computing k-means on FNC correlations ...');

%% Cluster
SPflat = cell2mat(SP);

clear SP;

try
    [IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, num_clusters, 'distance', dmethod, 'Replicates', 1, 'MaxIter', kmeans_max_iter, 'Display', 'iter', 'empty', 'drop');
catch
    [IDXp, Cp, SUMDp, Dp] = icatb_kmeans(SPflat, num_clusters, 'distance', dmethod, 'Replicates', 1, 'MaxIter', kmeans_max_iter, 'Display', 'iter', 'empty', 'drop');
end


clear SPflat IDXp SUMDp;

%% Use as a starting point to cluster all the data
FNCdynflat = reshape(FNCdynflat, M*Nwin , size(FNCdynflat, 3));
try
    [IDXall, Call, SUMDall, Dall] = kmeans(FNCdynflat, num_clusters, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
        'empty', 'drop', 'Start', Cp);
catch
    [IDXall, Call, SUMDall, Dall] = icatb_kmeans(FNCdynflat, num_clusters, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
        'empty', 'drop', 'Start', Cp);
end

clusterInfo.Call = Call;
clusterInfo.Cp = Cp;
clusterInfo.IDXall = IDXall;
clusterInfo.SUMDall = SUMDall;
clusterInfo.Dall = Dall;

icatb_save(post_process_file, 'clusterInfo', '-append');

icatb_save(param_file, 'dfncInfo');

disp('Done');

fprintf('\n');
