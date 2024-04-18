#!usr/bin/env python

class Sequence():
    '''
    Sequence generic class
    instantiate the elements of sequence object
    write method to calculate the GC content of the sequence
    write method for getting sequence id, sequence, length & quality line of sequence
    '''
    def __init__(self,seq_id,seq,info=None,qual=None):
        #initialize everything that the sequence object will contain
        #specify info and qual default values as none in case of fasta files
        #initialize seq_id, seq, info and qual
        self.seq_id = seq_id
        self.seq = seq
        self.info = info
        self.qual = qual

    def get_gc_content(self):
        '''
        Define a method to calculate the GC content of the sequence
        '''
        seq = self.seq
        total = len(seq)
        #count the no. of Gs and Cs and total them
        c = seq.count("C")
        g = seq.count("G")
        gc_total = g + c
        #Divide the total GC content with the length of the sequence
        #to get GC content of the sequence
        gc_content = gc_total / total
        return gc_content

    def get_length(self):
        '''
        Define a method to get the length of a sequence
        '''
        return len(self.seq)

    def get_sequence(self):
        '''
        Define a method to get the sequence from an object
        '''
        return self.seq

    def get_seqid(self):
        '''
        Define a method to get sequence id from an object
        '''
        return self.seq_id

    def get_qual(self):
        '''
        Define a method to get the quality line from an object
        '''
        return self.qual

    def get_qual_score(self):
        '''
        Define a method to calculate quality score of each sequence from the list of sequence objects
        '''
        #Create a dictionary containing key-value pairs from the two lists
        #the keys are Phred quality metrics in L1
        #the values are the scores assigned to each metric ranging from 0 to 42 in L2
        L1 = ['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5',
              '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
        L2 = list(range(42))
        qual_dic = dict(zip(L1, L2))
        #initialize the sum value to 0
        sum = 0
        #Initialize the average quality score to None
        avg_qual = None
        #Get the qual line of the read using the getter method defined in the module
        qual = self.get_qual()
        #loop over every character in the qual line
        for i in qual:
            #get the quality score of each character from the dictionary
            #add it to sum variable
            sum += qual_dic.get(i)
            #calculate the average qual by dividing the sum by length of the qual
        avg_qual = sum / len(qual)
        #return the value of average quality score
        return avg_qual


