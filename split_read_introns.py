#Command line
#Importing necessary modules:
import sys
import re

#Function to ensure the code will take SAM and annotations file as input:
def argument_parser():
    
    #Making sure there are only three arguments in input:
    if len(sys.argv) != 3:
        
        #Printing the error message and exiting:
        print("Format: python3 Script.py SamFile.sam Gene_Location_File.txt")
        sys.exit(1)
    
    #returns the SAM file and the Gene locations file:
    return sys.argv[1], sys.argv[2]

#Gene File
#Extracting the necessary information from the Gene_locations file:

#Function to Separate the required data in Gene_locations file:
def gene_file_parser(gene_file):
    
    #Initializing an empty list to store the required data:
    gene_info = []
    
    #Opening the Gene location file:
    with open(gene_file,"r") as file:
        for line in file:
            
            #Extracting gene_id and locations:
            columns = line.strip().split('\t')
            gene_id = columns[0]
            locations = columns[2]
            
            #Extracting the coordinates and chromosomes:
            if ':' in locations and '..' in locations:
                
                #Splitting the chromosomes and coordinates:
                chromosomes, coordinate_and_strand = locations.split(':')
                #Removing commas from the coordinates and splitting coordinates and strand:
                coordinates, strand = coordinate_and_strand.split('(')
                coordinates = coordinates.replace(',','')
                #Splitting the Coordinates Into Start and End
                start, end = map(int, coordinates.split('..'))
                strand = strand.strip(')')
                
                #Appending the necessary information from the gene file:
                gene_info.append({'gene_id':gene_id, 'chromosomes':chromosomes, 'start':start, 'end':end, 'strand':strand})
                
            #Wrong format informer:
            else:
                print(f"Skipping wrongly formatted location: {gene_id}")
    
    return gene_info

#SAM File
#Extracting necessary information from SAM file:

#Function to separate the required data from the SAM file:
def sam_file_parser(sam_file):
    
    #Initializing empty dictionary to store necessary data:
    sam_info = {}
    
    #Opening the SAM_file:
    with open(sam_file, 'r') as file:
        
        #Iterating over the file line by line:
        for line in file:
            
            #Checking if line starts with @ and skipping header:
            if line.startswith('@'):
                continue
            
            #Splitting the file into columns:
            columns = line.strip().split('\t')
            
            #Storing all chromosome names
            chromosomes = columns[2]
            #Storing all start positions as integers:
            position = int(columns[3])
            #Storing CIGAR strings:
            cigar = columns[5]
            #Storing NH column to check if reads are split:
            NH_nos = int(columns[-1].split(':')[2])
            
            #For reads that have aligned only once:
            if NH_nos == 1:
                
                #To see if read is split:
                if 'N' in cigar:
                    
                    #If chromosome does not exist in the dict then we add it:
                    if chromosomes not in sam_info:
                        sam_info[chromosomes] = {}
                        
                    #Initializing the startposition of read and Extracting the M/N:
                    current_position = position
                    matches = re.finditer(r'(\d+)([MIDNS])', cigar)
                    
                    #Iterating over the matches:
                    for match in matches:
                        #Extracting the length and alignment status of each match:
                        length, char = int(match.group(1)), match.group(2)
                        
                        #If already matched we just add the length
                        if char == 'M':
                            current_position += length
                            
                        #If skipped previously we update the positions and increase the count of supporting reads:
                        elif char == 'N':
                            start = current_position
                            end = current_position + length
                            
                            #Checks to see if the junctions already exist in dict:
                            if (start, end) not in sam_info[chromosomes]:
                                #if not we initialize the count:
                                sam_info[chromosomes][(start, end)] = 0
                                
                            #Increases the count every time a supporting read is found:
                            sam_info[chromosomes][(start, end)] += 1
                            current_position += length
                            
    return sam_info 

#Function to Map Junctions obtained from the SAM file to Gene locations:
def junction_to_gene_mapper(gene_info, sam_info):
    #Initializing a dicitionary to store the mapped junctions
    junctions = {}
    
    #Iterating every gene:
    for gene in gene_info:
        
        #Storing chromosome data:
        chromosomes = gene['chromosomes']
        #Storing start and end positions
        start, end = gene['start'], gene['end']
        #Storing gene id
        gene_id = gene['gene_id']
        
        #Creating empy list for genes:
        if gene_id not in junctions:
            junctions[gene_id] = []
    
    #Checking if chromosome has junctions:
        if chromosomes in sam_info:
            for(junction_start, junction_end), count in sam_info[chromosomes].items():
            
                #Checking if junctions exist within the length of gene:
                if start <= junction_start and junction_end <= end:
                
                    #Appending the junction information to the gene:
                    junctions[gene_id].append({'junction_start': junction_start, 'junction_end': junction_end, 'count': count })
                
    return junctions

#Output file
#Writing a function to print output in required format:
def output_writer(output_file, junctions):
    
    #Making sure the file writes:
    try:
    #Opening the file to write too:
        with open(output_file, 'w') as file:
        
        #Iterating over required data:
            for gene_id, junction_list in junctions.items():
                for junction in junction_list:
                
                    #Writing the data in required format:
                    file.write(f"{gene_id}\t{junction['junction_start']}\t{junction['junction_end']}\t{junction['count']}\n")
                #Adding empty line in between genes:
                file.write("\n")
                
    #Error display:
    except Exception as e:
        print(f"Error writing to file: {e}")
        
#Main Function:
def main():
    
    #Calling functions for command line prompt:
    sam_file, gene_file = argument_parser()
    print(f"Parsing gene file: {gene_file}")
    print(f"Parsing SAM file: {sam_file}")
    
    #Calling functions to preprocess the necessary files:
    gene_info = gene_file_parser(gene_file)
    print(f"Processed {len(gene_info)} genes.")
    sam_info = sam_file_parser(sam_file)
    print(f"Processed {len(sam_info)} junctions.")
    
    #Calling a function to perform the mapping between genes and their Junctions
    junctions = junction_to_gene_mapper(gene_info, sam_info)
    
    #Writing the output
    #Hardcoding the name of the file:
    outputfile = "Split_read_introns"
    output_file = f"{outputfile}.txt"
    
    #calling the output writer function
    output_writer(output_file, junctions)
    print(f"Output written to {output_file}")
    
if __name__ == "__main__":
    main()                

        
