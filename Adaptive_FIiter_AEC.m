clc;clear all;close all;
%% Input Parameters
order = 320;
miu = 0.01;
block_size = 2;
lambda = 0.1;
files = dir(['E:\AEC-Challenge\AEC-Challenge\datasets\real\','*_farend_singletalk_lpb.wav']);
len = length(files);
SER = 0;
%% Filter Step
tic
for i = 1:100
    % Path
    lpb_path = ['E:\AEC-Challenge\AEC-Challenge\datasets\real\', files(i).name];
    mic_path = strrep(lpb_path,'lpb','mic');
    err_path = strrep(lpb_path,'lpb','err');
    out_path = strrep(lpb_path,'lpb','out');
    
    if exist(mic_path) && exist(lpb_path)
        % Read wav
        [lpb_wav, ~] = audioread(lpb_path);
        [mic_wav, ~] = audioread(mic_path);

        % Align
        [lpb_wav, mic_wav] = align(lpb_wav, mic_wav);

        % Filt
%         [error, output, weights] = LMS(mic_wav, lpb_wav, miu, order);
%         [error, output, weights] = Block_LMS(mic_wav, lpb_wav, miu, order, block_size);
%         [error, output, weights] = NLMS(mic_wav, lpb_wav, order, 0.3, 0.1);
%         [error, output, weights] = VSS_LMS(mic_wav, lpb_wav, order);
%         [error, output, weights] = RLS(mic_wav, lpb_wav, order, lambda);
%         [error, output, weights] = AP(mic_wav, lpb_wav, miu, order,2,1e-4);
%         [error, output, weights] = PNLMS(mic_wav, lpb_wav, miu, order);
%         S = SAFinit(160,0.001,8,8,640);
%         [error,output] = SAFadapt(lpb_wav,mic_wav,S);
        [error, output, weights] =  FDAF(mic_wav,lpb_wav,0.001,0.001, order, 0);

        % plot
%         figure
%         subplot 211
%         plot(output,'r');hold on;
%         plot(mic_wav,'blue');hold on;
%         title(lpb_path)
%         legend('estmic','mic')
%         subplot 212
%         plot(error,'green')

        % SER
        10*log10(sum(mic_wav.^2)/(sum(error.^2)+1e-9))
        SER = SER + 10*log10(sum(mic_wav.^2)/sum(error.^2));
        % Write wav
        audiowrite(err_path, error, 16000);
%         audiowrite(out_path, output, 16000);
    end
end
toc
SER_Mean = SER ./ 100
