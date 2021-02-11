def read_fasta(filename):
    """Returns a list of tuples of each header and sequence in a fasta (or multifasta) file.
    first element in tuple is header and second the sequence.
    
    Key Arguments:
        filename -- fasta file.
    """
    tmp_seq = None
    seqs_list = []
    with open(filename, 'r') as fasta_file:
        for line in fasta_file:
            line = line.replace('\n','')
            if '>' in line:
                if tmp_seq != None:
                    seqs_list.append((hd, tmp_seq))
                tmp_seq = ''
                hd = line.replace('>','')
            else:
                tmp_seq += line
        seqs_list.append((hd, tmp_seq))
    try:
        assert len(seqs_list) > 0
    except AssertionError:
        print('The selected file is not a Fasta file.')
    else:
        return seqs_list

def write_fasta(outfile, seq_dict):
    """Writes fasta with dictionary where keys are headers and values sequences.
    
    Key Arguments:
        outfile.
    """
    step = 70
    with open(outfile, 'w') as file: 
        for header, sequence in seq_dict.items():
            sequence_list = [sequence[i - step: i] for i in range(step, len(sequence) + 1, step)]
            last = sequence[step * (len(sequence) // step):]
            if last != '':
                sequence_list.append(last)
            sequence = '\n'.join(sequence_list)
            file.write('>' + header + '\n' + sequence + '\n')     

def reads_generator(fasta_file, read_length, k):
    """This function simulates the reads generation from a fasta file with a coverage not less than 50.
    It will return a list of tuples. First element in tuple is read ID and second the sequence.
    
    Key Arguments:
        fasta_file -- fasta file.
        read_length -- size of reads.
    """
    reads_list = []
    overlap = k - 1
    input_header, input_seq = read_fasta(fasta_file)[0]
    n = len(input_seq)
    for i in range(0, n - overlap, read_length - overlap):
        read_seq = input_seq[i: i + read_length]
        reads_list.append(read_seq)
    return [('{}_{}'.format(input_header, i), read) for i, read in enumerate(reads_list)]

def write_fastq(reads_list, filename):
    """This function created a FASTQ file from a list of read generated by the reads_generator function.
    Key Arguments:
        reads_list -- list of reads generated with reads_generator.
        filename -- name of output file WITH EXTENSION.
    """
    with open(filename, 'w') as fastq_file:
        for read_id, read in reads_list:
            fastq_file.write('@{}\n'.format(read_id))
            fastq_file.write(read + '\n')
            fastq_file.write('+\n')
            fastq_file.write('I' * len(read) + '\n')  # max possible score

def read_fastq(filename):
    """This function reads a FASTQ file storing the read and its ID in a dictionary where keys are IDs and read value.
    This function does not consider + and score lines.
    
    Key Arguments:
        filename -- name of FASTQ input file.
    """
    reads_dict = dict()
    with open(filename, 'r') as fastq_file:
        for line in fastq_file:
            if '@' in line:
                reads_dict[line[1:].replace('\n', '')] = next(
                    fastq_file).replace('\n', '')
                next(fastq_file)
                next(fastq_file)
    return reads_dict