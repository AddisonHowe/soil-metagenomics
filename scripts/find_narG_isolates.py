import os
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_narG_sequences(annotations_folder, output_fasta):
    """
    Extract narG sequences from annotation files and create a FASTA file.
    
    Args:
        annotations_folder: Path to folder containing annotation .txt files
        output_fasta: Path for output FASTA file
    """
    
    narG_records = []
    
    # Get all text files in the annotations folder
    annotation_files = glob.glob(os.path.join(annotations_folder, "*.txt"))
    
    print(f"Found {len(annotation_files)} annotation files")
    
    for file_path in annotation_files:
        # Get the base filename (e.g., "ACM01" from "ACM01.txt")
        filename = os.path.basename(file_path)
        sample_id = filename.replace('.txt', '')
        
        try:
            with open(file_path, 'r') as f:
                headers = f.readline().strip().split('\t')
                
                # Find the column indices we need
                function_idx = headers.index('function')
                aa_sequence_idx = headers.index('aa_sequence')
                
                for line_num, line in enumerate(f, start=2):
                    fields = line.strip().split('\t')
                    
                    if len(fields) <= max(function_idx, aa_sequence_idx):
                        continue  # Skip incomplete lines
                    
                    function = fields[function_idx].lower()
                    aa_sequence = fields[aa_sequence_idx]
                    
                    # Check if this is a narG sequence
                    if ('nitrate reductase' in function and 
                        ('alpha' in function or 'narG' in function.lower())):
                        
                        # Skip if amino acid sequence is empty or too short
                        if aa_sequence and len(aa_sequence) > 10:
                            # Create sequence record
                            record = SeqRecord(
                                Seq(aa_sequence),
                                id=sample_id,
                                description=f"narG_{sample_id}_line_{line_num}"
                            )
                            narG_records.append(record)
                            print(f"Found narG in {filename}: {function}")
                            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
    
    # Write all narG sequences to FASTA file
    if narG_records:
        from Bio import SeqIO
        SeqIO.write(narG_records, output_fasta, "fasta")
        print(f"\nSuccessfully created {output_fasta} with {len(narG_records)} narG sequences")
    else:
        print("No narG sequences found!")
    
    return narG_records


if __name__ == "__main__":
    annotations_folder = "../data/isolate_seqs/annotations"  
    output_fasta = "../data/isolate_seqs/K00370_sequences.fasta"
    
    narG_sequences = extract_narG_sequences(annotations_folder, output_fasta)