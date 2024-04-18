#!usr/bin/env python
#Import all the modules in the main program
import sys
from FASTA_module import*
from FASTQ_module import *

if __name__ == '__main__':
    # Ask the user to enter a filename or filepath for processing
    x = input("Enter a filename/filepath for processing: ")
    #insert an error handling method to check if the file is valid
    try:
        #check if the input file ends with .fasta
        if x.endswith('.fasta'):
            #if the above condition is true, create an fasta object by calling the Fasta function on input file
            fasta_obj = Fasta(x)
            #Call the read function from FASTA module to read the fasta object
            fasta_obj.read_fasta(x)
            #create a list containing all nucleotides found in DNA/RNA fasta file
            nuc = ['A', 'T', 'G', 'C', 'N']
            #ask the user to input a minimum length threshold for the output sequences
            len_thresh = int(input("Specify a minimum length threshold: "))
            #print the metrics of the file
            fasta_obj.print_metrics()
            #to check if the fasta file is a nucleotide or a protein file
            #get the first read from the fasta object
            first_read = fasta_obj.fasta_reads[0]
            #get the sequence from the first read
            frst_seq = first_read.get_sequence()
            #initially assign the type of file to N to specify nucleotide
            fasta_file_type = 'N'
            #loop over every character in the frst_seq
            for i in frst_seq:
                #check if it is not in the nuc list
                if i not in nuc:
                    #if the condition is true, assign the type of file to P to specify peptide
                    fasta_file_type = 'P'
                    #break the loop after initial check
                    break
            #if the file type is of N print the average GC content of the file
            if fasta_file_type == 'N':
                print("The average GC content is: ", fasta_obj.average_gc())

            #Ask the user if he/she is specifying an output file to write the output sequences
            y = input("Do you want to write to an output file?: ")
            # if the answer is yes
            if y == "Y":
                #ask for the file path
                z = input("Provide file name or path: ")
                #write the sequences to file filtering the minimum length threshold provided by the user
                fasta_obj.write_output(z, len_thresh)
            #if there is no file
            elif y == "N":
                #loop over each read in the fasta object
                for read in fasta_obj.fasta_reads:
                    #check if the read length is greater than or equal to the threshold value
                    if read.get_length() >= len_thresh:
                        #print the seq_id and seq which pass the threshold criterion
                        print(read.get_seqid(), '\n', read.get_sequence())

        #if the file ends with .fastq, execute the below code
        elif x.endswith('.fastq'):
            #create an fastq object by calling the Fastq function on input file
            fastq_obj = Fastq(x)
            # Call the read function from FASTQ module to read the fastq object
            fastq_obj.read_fastq(x)
            #print the basic metrics of the fastq file
            fastq_obj.print_metrics()
            # ask the user to input a minimum length threshold for the output sequences
            len_thresh = int(input("Specify a minimum length threshold: "))
            # ask the user to input a minimum quality threshold value for the output sequences
            qual_thresh = float(input("Enter minimum quality threshold: "))
            #ask the user if he wants to write the data to an output file

            y = input("Do you want to write to an output file?: ")
            #if the answer is Y
            if y == "Y":
                #ask for the file path
                z = input("Provide file path: ")
                #write the sequences that have passed the length and quality threshold criteria to the output file
                fastq_obj.write_output(z, len_thresh, qual_thresh)
                #if the user is not specifying an output file
            elif y == 'N':
                #loop over the reads in the fastq object
                for read in fastq_obj.fastq_reads:
                    #check if the read passes the quality and length threshold criteria
                    if (read.get_qual_score() >= qual_thresh) and  (read.get_length() >= len_thresh):
                        #print the sequences that pass the above criteria
                        print(read.get_seqid(), '\n', read.get_sequence())

    #provide exception if the correct file name or file path is not entered
    except FileNotFoundError:
        print("The file that you entered does not exist. Please enter a valid file name with its absolute path")
    #handle the exception when str or special characters are entered instead of a numeric value
    except ValueError:
        print('Please enter a valid positive number')
    #handle exception for all other types of errors
    except:
        print('An unexpected error occured')
















