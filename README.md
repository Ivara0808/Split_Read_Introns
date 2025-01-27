# Split_Read_Introns
Split Read Introns is a Python script designed to analyze SAM alignment files and identify split reads that span introns. By parsing CIGAR strings and alignment positions, the script detects intronic junctions, maps them to genes, and counts the number of supporting reads. The results are saved in a structured output file.

Run the script using:
python3 split_read_introns.py mySamFile.sam myGeneTable.txt

