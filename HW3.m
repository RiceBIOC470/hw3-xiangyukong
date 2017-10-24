% GB comments
1.	80 Need to write out the final alignment 
2a. 100
2b. 100
2c. 100
3a 100 
3b. 100
3c. 100  	
Overall: 97

%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

MAPK3_stat=fastaread('sequence0.fasta');
ERK2_stat=fastaread('sequence1.fasta');
seq1=MAPK3_stat.Sequence;
seq2=ERK2_stat.Sequence;
[MAPK3_ORF,start1,stop1]=findORF(seq1);
[ERK2_ORF,start2,stop2]=findORF(seq2);
seqr_MAP=seq1(start1:stop1);
seqr_ERK=seq2(start2:stop2);

[score,align,start]=swalign(seqr_MAP,seqr_ERK,'Alphabet','nt','Showscore','true');

ERK1(825:948)

% Part2. Perform an of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

MAPK3_protein=aminolookup(dna2protein(seqr_MAP,1));
ERK2_protein=aminolookup(dna2protein(seqr_ERK,1));
[score,align,start]=swalign(MAPK3_protein,ERK2_protein,'Alphabet','aa','Showscore','true')
%the last 1/3 of the ORF somewhat align

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 
 mRNA_mouse_ERK1=fastaread('mouseERK1mRNA.fasta');
 mRNA_mouse_ERK2=fastaread('mouseERK2mRNA.fasta');
 aa_mouse_ERK1=fastaread('mouseERK1protein.fasta');
 aa_mouse_ERK2=fastaread('mouseERK2protein.fasta');
 protein_mouse1=aa_mouse_ERK1.Sequence;
 protein_mouse2=aa_mouse_ERK2.Sequence;
[mouse1_ORF,start3,stop3]=findORF(mRNA_mouse_ERK1.Sequence);
[mouse2_ORF,start4,stop4]=findORF(mRNA_mouse_ERK2.Sequence);
seqr_mouse1=mRNA_mouse_ERK1.Sequence(start3:stop3);
seqr_mouse2=mRNA_mouse_ERK2.Sequence(start4:stop4);
[score,align,start]=swalign(seqr_MAP,seqr_mouse1,'Alphabet','nt','Showscore','true');
[score,align,start]=swalign(seqr_ERK,seqr_mouse2,'Alphabet','nt','Showscore','true');
[score,align,start]=swalign(MAPK3_protein,protein_mouse1,'Alphabet','aa','Showscore','true');
[score,align,start]=swalign(ERK2_protein,protein_mouse2,'Alphabet','aa','Showscore','true');
%ERK1 has a match on mRNA for about 120bp and 25aa; ERK2 has a match on
%mRNA for about 350bp and 75aa.
%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

 accession='NM_002746';
 topmatch(accession,3)

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 



% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 
NM_002746
human_match:{'NM_002746.2'}
non_human_match:{'X60188.1'}

AC_000171
%can't find sequence....



