function result=species_match(accession);
gb_data=getgenbank(accession);
seq_match=gb_data.Sequence;
[requestID,requestTime]=blastncbi(seq_match,'blastn');
blast_data=getblast(requestID,'WaitTime',requestTime);
hm='';
nhm='';
n=1
while n<51
if contains(blast_data.Hits(n).Name,'Homo sapiens')||contains(blast_data.Hits(n).Name,'human');
    extract=strsplit(blast_data.Hits(n).Name,'|');
    hm=extract(4);
    return;
    if n>50;
        hm='No good human match';
    end
else n=n+1
end
end
i=1;
while i<51;
if  contains(blast_data.Hits(i).Name,'Homo sapiens')||contains(blast_data.Hits(i).Name,'human');
    i=i+1;
else extract1=strsplit(blast_data.Hits(i).Name,'|');
    nhm=extract1(4);return;
end
end
result=['human_match',hm,'non_human_match',nhm];
end


 

