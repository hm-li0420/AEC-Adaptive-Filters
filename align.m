function [wav1, wav2] = align(wav1, wav2)
    if length(wav1) > length(wav2)
        wav1 = wav1(1:length(wav2));
    elseif length(wav1) < length(wav2)
        wav2 = wav2(1:length(wav1));
    end
end