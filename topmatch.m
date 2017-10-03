function result=topmatch(accession,N);
gb_data=getgenbank(accession);
seq_match=gb_data.Sequence;
[requestID,requestTime]=blastncbi(seq_match,'blastn');
blast_data=getblast(requestID,'WaitTime',requestTime);
Hit_accession={};
datawant=blast_data.Hits(1:N);
for ii=1:N;
    extract=strsplit(datawant(ii).Name,'|');
    Hit_accession{1,ii}=extract(4);
end
result=Hit_accession;
end
