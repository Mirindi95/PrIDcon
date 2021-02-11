class Contig:
    """Class containing all sequences related with a contig and some annotations."""
    
    codon_code = {"UUU":"F","UUC":"F","UUA":"L","UUG":"L","UCU":"S","UCC":"S","UCA":"S","UCG":"S","UAU":"Y","UAC":"Y","UAA":"_","UAG":"_","UGU":"C","UGC":"C","UGA":"_","UGG":"W","CUU":"L","CUC":"L","CUA":"L","CUG":"L","CCU":"P","CCC":"P","CCA":"P","CCG":"P","CAU":"H","CAC":"H","CAA":"Q","CAG":"Q","CGU":"R","CGC":"R","CGA":"R","CGG":"R","AUU":"I","AUC":"I","AUA":"I","AUG":"M","ACU":"T","ACC":"T","ACA":"T","ACG":"T","AAU":"N","AAC":"N","AAA":"K","AAG":"K","AGU":"S","AGC":"S","AGA":"R","AGG":"R","GUU":"V","GUC":"V","GUA":"V","GUG":"V","GCU":"A","GCC":"A","GCA":"A","GCG":"A","GAU":"D","GAC":"D","GAA":"E","GAG":"E","GGU":"G","GGC":"G","GGA":"G","GGG":"G"} # Stop codons are represented with the underscore (_) character.
    
    def __init__(self, contig_ID, sequence, splicing = True, minlength_orf = 30):
        self.__contig_ID = contig_ID
        self.__DNA = sequence.upper()
        self.__cDNA = self.get_cDNA(self.__DNA)
        self.__mRNA = self.get_mRNA(splicing)
        self.__reading_frames = self.get_readingframes()
        self.__ORFs = self.get_ORFs(minlength_orf)
        
    # ---------- returning attributes ----------
    
    def contig_sequence(self):
        return {self.__contig_ID: self.__DNA}
    
    def mrna_sequence(self):
        return self.__mRNA
    
    def found_orfs(self):
        return self.__ORFs
    
    # ---------- instance methods ----------
    
    def get_mRNA(self, splicing):
        """Returns the full transcription from both forward and reverse contig strands."""
        f_RNA = self.__DNA.translate(str.maketrans({"T": "U"}))
        r_RNA = self.__cDNA.translate(str.maketrans({"T": "U"}))
        if splicing:
            f_RNA = Contig.splice_sequence(f_RNA)
            r_RNA = Contig.splice_sequence(r_RNA)
        return {self.__contig_ID + '_forward': f_RNA, self.__contig_ID + '_reverse': r_RNA}
    
    def get_readingframes(self):
        """Returns all six reading frames in a list."""
        return self.readingframes(self.__mRNA[self.__contig_ID + '_forward']) + self.readingframes(self.__mRNA[self.__contig_ID + '_reverse'])
    
    def get_ORFs(self, orf_length):
        """This function returns a dictionary of possible ORFs, it uses ORF_ID as key and sequence as value."""
        orf_full = []
        for reading_frame in self.__reading_frames:
            orf_full += self.openreadingframes(reading_frame, orf_length)
        orf_full.sort(key=len, reverse=True)
        return {'{}-ORF_{}'.format(self.__contig_ID, i): orf for i, orf in enumerate(orf_full)}
    
    # ---------- static methods ----------
    
    @staticmethod
    def get_cDNA(sequence):
        """Returns the cDNA sequence in 5' - 3' orientation.
        
        Key Arguments:
            sequence -- DNA sequence.
        """
        return sequence[::-1].translate(str.maketrans("ATGC", "TACG"))
    
    @staticmethod
    def splice_sequence(rna_sequence):
        """This function deletes introns of rna_sequence with the recognition sitas 5'-GU-intron-AG-3' and returns the resulting sequence (mature rna)
        
        Key Argument:
            rna_sequence -- sequence of RNA.
        """
        start_introns, stop_introns = Contig.position_subsequences(rna_sequence, 'GU', 'AG')
        # min intron length = 30 doi: 10.1093/dnares/dsv028
        intron_sequences = Contig.extract_subsequence(rna_sequence, start_introns, stop_introns, 30)
        for i in intron_sequences:
            rna_sequence = rna_sequence.replace(i, '')
        return rna_sequence
    
    @staticmethod
    def readingframes(RNA_sequence):
        """returns a list of all possible reading frames of RNA_sequence"""
        output = [] # list of all three reading frames
        for f in range(3): # f is the "frame offset" (0, 1 and 2) 
            f_seq = '' # frame sequence
            for i in range(2 + f, len(RNA_sequence), 3): # 2 + f is the 3' codon end
                codon = RNA_sequence[i - 2: i + 1]
                f_seq += Contig.codon_code[codon]
            output.append(f_seq)
        return output
         
    @staticmethod
    def openreadingframes(polypep_seq, min_length):
        """Return a list of all possible reading frames (including nested) of the polypeptide sequence polypep_seq with minimum length."""
        min_length += 1 # consider the stop codon from the Orf length
        start_codons, stop_codons = Contig.position_subsequences(polypep_seq, 'M', '_')
        orf_list = Contig.extract_subsequence(polypep_seq, start_codons, stop_codons, min_length)
        return [orf[:-1] for orf in orf_list] # -1 to deleted the stop symbol '_'
    
    @staticmethod
    def extract_subsequence(sequence, head_list, tail_list, min_length = 0):
        """This function extracts the substrings delimited by the indices present in head_list and tail_list. It does not extract those sequences contained in a bigger subsequence. The output is returned as a list of substrings.
        
        Example:
            sequence = ABCDEFGHIJ
            head_list = [0, 3]
            tail_list = [5]
            output = ABCDEF (DEF not returned cause is contained in the output)
            
        Key Arguments:
            sequence -- string sequence.
            head_list -- indices list of subsequence head.
            tail_list -- indices list of subsequence tail.
            min_length -- minimum subsequence length.
        """
        output = []
        last_tail = None
        for head in head_list:
            tail_iter = iter(tail_list)
            current_tail = next(tail_iter, 'end')
            # Find the closest stop codon after start_c
            while current_tail != 'end' and current_tail <= head:
                current_tail = next(tail_iter, 'end')
            if current_tail == 'end':
                break
            # check if the current tails haven't been used for a longer subsequence and has the minimum length
            if current_tail != last_tail:
                substring = sequence[head: current_tail + 1]
                if len(substring) >= min_length:
                    output.append(substring)
            last_tail = current_tail
        return output
    
    @staticmethod
    def position_subsequences(sequence, start, stop):
        """Returns a tuple of list with positions of start and stop substrings respectively.
        
        Key Arguments:
            sequence -- string sequence.
            start -- head of searched substring.
            stop -- tail of searched substring.
        """
        pos_start = []
        pos_stop = []
        n_start = len(start)
        n_stop = len(stop)
        for i in range(len(sequence)):
            if sequence[i: i + n_start] == start:
                pos_start.append(i)
            if sequence[i: i + n_stop] == stop:
                pos_stop.append(i + n_stop - 1) # index of right limit
        return pos_start, pos_stop
