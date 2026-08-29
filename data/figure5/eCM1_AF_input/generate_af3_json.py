import json
from pathlib import Path

# Define species and habitat mapping for job naming
mag_info = {
    "human706": ("S_aureus", "Skin"),
    "human199": ("S_epidermidis", "Skin"),
    "human545": ("C_accolens", "Skin"),
    "human491": ("D_pigrum", "Skin"),
    "human76": ("C_acnes", "Skin"),
    "human3694": ("S_aureus", "Nose"),
    "human3700": ("S_epidermidis", "Nose"),
    "human3518": ("C_accolens", "Nose"),
    "human3566": ("D_pigrum", "Nose"),
    "human3553": ("C_acnes", "Nose")
}

script_dir = Path(__file__).resolve().parent
fasta_file = script_dir / "metNIQ_sequences.fasta"
output_json = script_dir / "metNIQ_AF3_batch.json"

# Parse sequences
sequences = {} # MAG_ID -> {metQ: seq, metI: seq, metN: seq}
if fasta_file.exists():
    with open(fasta_file, 'r') as f:
        current_mag = None
        current_role = None
        for line in f:
            if line.startswith(">"):
                # Header format: >human706_metQ_p1899
                parts = line[1:].strip().split("_")
                current_mag = parts[0]
                current_role = parts[1]
                if current_mag not in sequences:
                    sequences[current_mag] = {}
            else:
                sequences[current_mag][current_role] = line.strip()

# Build JSON structure
jobs = []
for mag_id, (species, habitat) in mag_info.items():
    if mag_id not in sequences:
        print(f"Warning: No sequences found for {mag_id}")
        continue
    
    seq_data = sequences[mag_id]
    
    # Check if all components are present
    if not all(role in seq_data for role in ["metQ", "metI", "metN"]):
        print(f"Warning: Incomplete operon for {mag_id}")
        continue

    job = {
        "name": f"{species}_{habitat}_metNIQ",
        "modelSeeds": [],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": seq_data["metN"],
                    "count": 2
                }
            },
            {
                "proteinChain": {
                    "sequence": seq_data["metI"],
                    "count": 2
                }
            },
            {
                "proteinChain": {
                    "sequence": seq_data["metQ"],
                    "count": 1
                }
            }
        ],
        "dialect": "alphafoldserver",
        "version": 1
    }
    jobs.append(job)

with open(output_json, 'w') as f:
    json.dump(jobs, f, indent=2)

print(f"Successfully generated AF3 batch JSON with {len(jobs)} jobs at {output_json}")


