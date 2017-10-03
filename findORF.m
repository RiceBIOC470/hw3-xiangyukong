function [ORFlength, startpos, stoppos] = findORF(dnaseq)
newseq=upper(dnaseq)
startcodon=strfind(newseq,'ATG');
stopcodon=[strfind(newseq,'TAA'),strfind(newseq,'TAG'),strfind(newseq,'TGA')];
%Function to find the length of the longest open reading frame of a
%sequences called dnaseq

%fill in here.
%I still can't quite understand hw1, so I used ouyang27's code.
firstStopCodon = zeros(1, length(startcodon));
N = length(newseq);
for ii = 1:length(startcodon)
    Olengths = stopcodon - startcodon(ii)+3;
    goodlength = 1e8;
    goodind = 0;
    for jj = 1:length(Olengths)
        if Olengths(jj) > 0 && ...
                mod(Olengths(jj),3) == 0 && ...
                Olengths(jj) < goodlength
            goodlength = Olengths(jj);
            goodind = jj;
        end
    end
    if goodind > 0
        firstStopCodon(ii) = stopcodon(goodind);
    else
        firstStopCodon(ii) = startcodon(ii);
    end
end
ORF = firstStopCodon - startcodon + 3;
[M, I] = max(ORF);
if M > 3
    ORFlength = M;
    startpos = startcodon(I);
    stoppos = firstStopCodon(I);
else
    ORFlength = NaN;
    startpos = NaN;
    stoppos = NaN;
end