"""
============================================================
🧬 UNIFIED BIOINFORMATICS ANALYSIS PIPELINE
============================================================
Complete single-file version for Google Colab
Supports: FASTA, GenBank, PDB, Phylogenetic Trees
Author: Arin DAS
============================================================
"""

# ============================================================
# STEP 1: SETUP & INSTALLATION
# ============================================================
print("📦 Installing required packages...")
import subprocess
import sys

packages = ['biopython', 'py3Dmol', 'logomaker', 'matplotlib', 'pandas']
for package in packages:
    try:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', package])
    except subprocess.CalledProcessError:
        print(f"⚠️ Could not install {package}. If you're running locally, install it manually.")
print("✅ Installation complete!\n")

# ============================================================
# IMPORTS
# ============================================================
import os
import io
import re
import tempfile
import warnings
warnings.filterwarnings("ignore")

from Bio import SeqIO, Entrez, Phylo, motifs, AlignIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.Align import PairwiseAligner, MultipleSeqAlignment
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

import pandas as pd
import matplotlib.pyplot as plt
import py3Dmol
from IPython.display import display

try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except Exception:
    LOGOMAKER_AVAILABLE = False

# Try to import google.colab.files but allow local fallback
try:
    from google.colab import files
    COLAB_AVAILABLE = True
except Exception:
    COLAB_AVAILABLE = False
    # Provide a simple fallback for local environments
    def _local_upload_prompt():
        print("⚠️ Not in Colab. Please place files in the working directory and enter filenames when prompted.")
        return {}

# Set NCBI email
Entrez.email = "dasarin007@gmail.com"

# ============================================================
# UTILITY FUNCTIONS
# ============================================================
def print_header(title):
    """Print formatted header"""
    print("\n" + "="*60)
    print(f"🧬 {title}")
    print("="*60)

# ============================================================
# FILE HANDLING
# ============================================================
class FileHandler:
    """Handles file operations and type detection"""

    @staticmethod
    def detect_file_type(file_name, file_bytes):
        """Detect biological file type (improved heuristics)."""
        name = file_name.lower()

        # Extension-based detection (explicit)
        if name.endswith(".pdb"):
            return "protein_structure"
        if name.endswith((".fasta", ".fa", ".faa", ".fsa", ".fas")):
            return "sequence"
        if name.endswith((".gb", ".gbk", ".genbank")):
            return "genbank"
        if name.endswith((".nwk", ".newick", ".tree")):
            return "phylogenetic_tree"

        # Content-based detection (safer)
        try:
            text = file_bytes.decode('utf-8', errors='ignore').lstrip()
            first_line = text.splitlines()[0] if text else ""

            # FASTA
            if first_line.startswith(">"):
                return "sequence"

            # GenBank (LOCUS or ACCESSION lines)
            if first_line.upper().startswith("LOCUS") or "ACCESSION" in first_line.upper():
                return "genbank"

            # Newick: must contain parentheses and end with semicolon
            if "(" in text and ")" in text and ";" in text:
                # quick sanity check for Newick tokens
                if re.search(r"[^\s\(\),:;]+:[0-9Ee\.\-]+", text) or text.strip().endswith(";"):
                    return "phylogenetic_tree"

            # PDB heuristics (HEADER, ATOM lines)
            if first_line.upper().startswith("HEADER") or "ATOM  " in text or "HETATM" in text:
                return "protein_structure"

        except Exception:
            pass

        return "unknown"

    @staticmethod
    def parse_sequence(bytes_obj):
        """Parse FASTA sequence file"""
        try:
            return list(SeqIO.parse(io.StringIO(bytes_obj.decode(errors='ignore')), "fasta"))
        except Exception as e:
            print(f"❌ Error parsing FASTA: {e}")
            return []

    @staticmethod
    def parse_genbank(bytes_obj):
        """Parse GenBank — returns the first record if multiple are present."""
        stream = io.StringIO(bytes_obj.decode(errors='ignore'))
        # Try reading as single record first, else as iterator
        try:
            rec = SeqIO.read(stream, "genbank")
            return rec
        except Exception:
            # rewind and try iterator
            stream.seek(0)
            records = list(SeqIO.parse(stream, "genbank"))
            if not records:
                raise ValueError("No GenBank records found")
            if len(records) > 1:
                print("⚠️ Multiple GenBank records detected — using the first one for interactive analysis.")
            return records[0]

# ============================================================
# SEQUENCE ANALYSIS
# ============================================================
class SequenceAnalyzer:
    """Analyzes biological sequences"""

    @staticmethod
    def detect_sequence_type(seq):
        """Robustly determine if sequence is DNA, RNA, or Protein (heuristic)."""
        seq_str_raw = str(seq)
        seq_str = seq_str_raw.upper()
        # Quick check for empty
        if not seq_str.strip():
            return "Unknown"

        # Treat U as T for heuristic but preserve detection of RNA if U present without T
        has_u = "U" in seq_str and "T" not in seq_str
        if has_u:
            return "RNA"

        seq_chars = set(re.sub(r'[^A-Z]', '', seq_str))
        if not seq_chars:
            return "Unknown"

        dna_chars = set("ATGCN")
        # If the set is subset of DNA chars -> DNA
        if seq_chars.issubset(dna_chars):
            return "DNA"

        # If many ambiguous letters or B, J, O, X, Z -> treat as Protein
        if any(ch in seq_chars for ch in set("BJOUXZ")):
            return "Protein"

        # If there's >20% letters outside DNA, call Protein
        non_dna_fraction = len([c for c in seq_chars if c not in dna_chars]) / max(1, len(seq_chars))
        if non_dna_fraction > 0.2:
            return "Protein"

        # Fallback: amino-acid letters
        amino_letters = set("ARNDCQEGHILKMFPSTWYVXBZJOU")
        if seq_chars.issubset(amino_letters):
            return "Protein"

        # Default guess
        return "DNA"

    @staticmethod
    def analyze_protein_basic(seq_record):
        """Basic protein analysis"""
        seq = str(seq_record.seq).replace("*", "").replace("-", "")
        analyzer = ProteinAnalysis(seq)

        return {
            "id": seq_record.id,
            "length": len(seq),
            "mw": analyzer.molecular_weight(),
            "pI": analyzer.isoelectric_point(),
            "gravy": analyzer.gravy(),
            "aromaticity": analyzer.aromaticity(),
            "instability": analyzer.instability_index(),
            "aa_comp": analyzer.get_amino_acids_percent(),
            "secondary_structure": analyzer.secondary_structure_fraction()
        }

    @staticmethod
    def detailed_analysis(record):
        """Print detailed sequence analysis"""
        seq_type = SequenceAnalyzer.detect_sequence_type(record.seq)

        print(f"\n{'='*60}")
        print(f"📊 Detailed Analysis: {record.id}")
        print(f"Type: {seq_type}")
        print(f"Length: {len(record.seq)}")

        if seq_type in ["DNA", "RNA"]:
            try:
                gc = gc_fraction(record.seq) * 100
                mw = molecular_weight(record.seq, seq_type)
                print(f"GC Content: {gc:.2f}%")
                print(f"Molecular Weight: {mw:.2f} Da")
            except Exception as e:
                print(f"⚠️ Could not compute GC/MW: {e}")

            # Base counts (safely)
            seq_str = str(record.seq).upper()
            counts = {base: seq_str.count(base) for base in "ATGC"}
            print(f"Base composition: {counts}")

        elif seq_type == "Protein":
            stats = SequenceAnalyzer.analyze_protein_basic(record)
            print(f"Molecular Weight: {stats['mw']:.2f} Da")
            print(f"Isoelectric Point: {stats['pI']:.2f}")
            print(f"Aromaticity: {stats['aromaticity']:.3f}")
            print(f"Instability Index: {stats['instability']:.2f}")
            print(f"GRAVY: {stats['gravy']:.3f}")
            helix, turn, sheet = stats['secondary_structure']
            print(f"Secondary Structure:")
            print(f"  Helix: {helix:.2%}, Turn: {turn:.2%}, Sheet: {sheet:.2%}")

    @staticmethod
    def identify_motifs(seq_record):
        """Identify biological motifs in sequence"""
        seq = str(seq_record.seq)

        motifs_dict = {
            "N-glycosylation site": r"N[^P][ST][^P]",
            "Protein kinase C phosphorylation": r"[ST].[RK]",
            "Casein kinase II phosphorylation": r"[ST]..[DE]",
            "Tyrosine kinase phosphorylation": r"[RK].{2,3}[YFW]",
            "ATP-binding site": r"G.{4}GK[ST]",
            "Zinc finger": r"C.{2,4}C.{12}H.{3,5}H",
            "Nuclear localization signal": r"[KR]{3,5}",
        }

        print(f"\n🔍 Motif Search: {seq_record.id}")
        print("-" * 60)

        found_any = False
        for name, pattern in motifs_dict.items():
            matches = [(m.start(), m.group()) for m in re.finditer(pattern, seq)]
            if matches:
                found_any = True
                print(f"\n✅ {name}: {len(matches)} match(es)")
                for pos, match in matches[:3]:
                    print(f"   Position {pos+1}: {match}")
                if len(matches) > 3:
                    print(f"   ... and {len(matches)-3} more")

        if not found_any:
            print("⚠️  No common motifs detected")

    @staticmethod
    def translate_sequence(record):
        """Translate DNA/RNA to protein"""
        proteins = []

        if hasattr(record, 'features'):
            try:
                cds_features = [f for f in record.features if f.type == "CDS"]
                if cds_features:
                    for f in cds_features:
                        try:
                            cds_seq = f.extract(record.seq)
                            prot = cds_seq.translate(to_stop=True)
                            product = f.qualifiers.get("product", ["Unknown"])[0]
                            proteins.append((product, prot))
                        except Exception:
                            pass
                    return proteins
            except Exception:
                # If features can't be iterated, fallback to full translation
                pass

        try:
            prot = record.seq.translate(to_stop=True)
            proteins.append(("Full sequence translation", prot))
        except Exception as e:
            proteins.append((f"Error: {e}", None))

        return proteins

    @staticmethod
    def align_sequences(records):
        """Perform multiple sequence alignment or pairwise scores"""
        if len(records) < 2:
            print("❌ Need at least 2 sequences for alignment")
            return None

        print(f"\n🔗 Aligning {len(records)} sequences...")

        if len(records) == 2:
            # Pairwise alignment
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5

            alignments = list(aligner.align(records[0].seq, records[1].seq))
            best = alignments[0]

            print(f"Score: {best.score}")
            print(f"\n{str(best)[:500]}")  # Show first 500 chars
            return None
        else:
            # Multiple sequence alignment (simple clustering)
            print("⚠️  For MSA with >2 sequences, use external tools like MUSCLE/ClustalW")
            print("Showing pairwise distances instead...")

            distances = []
            for i in range(len(records)):
                for j in range(i+1, len(records)):
                    aligner = PairwiseAligner()
                    score = list(aligner.align(records[i].seq, records[j].seq))[0].score
                    distances.append((records[i].id, records[j].id, score))

            print("\nPairwise alignment scores:")
            for id1, id2, score in distances:
                print(f"  {id1} <-> {id2}: {score:.1f}")

            return distances

# ============================================================
# PHYLOGENETIC ANALYSIS
# ============================================================
class PhylogeneticAnalyzer:
    """Handles phylogenetic tree construction and visualization"""

    @staticmethod
    def build_tree_from_sequences(records, method="upgma"):
        """Build phylogenetic tree from sequences (works with unaligned sequences)."""
        if len(records) < 3:
            print("❌ Need at least 3 sequences to build a tree")
            return None

        print(f"\n🌳 Building phylogenetic tree using {method.upper()} (fast, heuristic)...")
        try:
            names = [r.id for r in records]
            n = len(records)

            # Compute pairwise distances via normalized pairwise alignment score
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5

            # Create lower-triangular distance matrix expected by Bio.Phylo.TreeConstruction.DistanceMatrix
            matrix = []
            for i in range(n):
                row = []
                seq_i = str(records[i].seq)
                for j in range(i):
                    seq_j = str(records[j].seq)
                    # Align and get best score
                    try:
                        score = list(aligner.align(seq_i, seq_j))[0].score
                    except Exception:
                        score = 0.0
                    # Normalise score by max possible length*match_score to get identity-like fraction
                    max_len = max(len(seq_i), len(seq_j), 1)
                    # Convert to distance: higher score -> smaller distance
                    similarity = score / (2.0 * max_len)  # approx scaling
                    # Clip to [0,1]
                    similarity = max(0.0, min(1.0, similarity))
                    distance = round(1.0 - similarity, 6)
                    row.append(distance)
                matrix.append(row)

            # Build DistanceMatrix and tree
            from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
            dm = DistanceMatrix(names, matrix)
            constructor = DistanceTreeConstructor()
            if method.lower() == "upgma":
                tree = constructor.upgma(dm)
            else:
                tree = constructor.nj(dm)

            # Save and plot
            os.makedirs("results", exist_ok=True)
            out = os.path.join("results", "phylogenetic_tree.nwk")
            Phylo.write(tree, out, "newick")
            print(f"💾 Tree saved to: {out}")

            PhylogeneticAnalyzer.plot_tree(tree_obj=tree)
            return tree

        except Exception as e:
            print(f"❌ Error building tree: {e}")
            return None

    @staticmethod
    def plot_tree(tree_obj=None, newick_text=None):
        """Plot phylogenetic tree"""
        try:
            if newick_text:
                tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".nwk", mode="w")
                tmp.write(newick_text)
                tmp.close()
                tree = Phylo.read(tmp.name, "newick")
                os.unlink(tmp.name)
            else:
                tree = tree_obj

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, axes=ax, do_show=False)
            plt.title("Phylogenetic Tree", fontsize=14, fontweight='bold')
            plt.tight_layout()
            plt.show()

        except Exception as e:
            print(f"❌ Error plotting tree: {e}")

# ============================================================
# STRUCTURE ANALYSIS
# ============================================================
class StructureAnalyzer:
    """Analyzes protein structures"""

    @staticmethod
    def parse_pdb_data(file_bytes, file_name):
        """Parse PDB file"""
        tmpdir = tempfile.mkdtemp()
        fpath = os.path.join(tmpdir, file_name)
        with open(fpath, "wb") as f:
            f.write(file_bytes)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(file_name.split('.')[0], fpath)

        atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atoms.append([
                            chain.id, residue.get_resname(), residue.id[1],
                            atom.get_name(), *atom.get_coord()
                        ])

        df = pd.DataFrame(atoms, columns=["Chain", "Residue", "ResID", "Atom", "X", "Y", "Z"])
        return structure, df, fpath

    @staticmethod
    def visualize_structure_py3dmol(file_bytes):
        """Visualize 3D structure using py3Dmol"""
        print("🧩 Rendering 3D structure...")

        # Read PDB content
        pdb_content = file_bytes.decode('utf-8', errors='ignore')

        # Create py3Dmol view
        view = py3Dmol.view(width=800, height=600)
        view.addModel(pdb_content, 'pdb')

        # Add multiple styles for better visualization
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        # Try to highlight chain A if present
        try:
            view.addStyle({'chain': 'A'}, {'cartoon': {'color': 'blue', 'opacity': 0.8}})
        except Exception:
            pass

        view.zoomTo()
        view.spin(False)  # Disable auto-spin

        print("✅ 3D structure loaded. You can rotate and zoom with mouse.")

        return view

# ============================================================
# DATABASE TOOLS
# ============================================================
class DatabaseTools:
    """Handles database operations"""

    @staticmethod
    def fetch_record(acc_id, db="nucleotide"):
        """Fetch sequence from NCBI"""
        print(f"🔍 Fetching {acc_id} from {db}...")
        rettype = "gb" if db == "nucleotide" else "fasta"
        handle = Entrez.efetch(db=db, id=acc_id, rettype=rettype, retmode="text")
        record = SeqIO.read(handle, rettype)
        handle.close()
        print(f"✅ Retrieved {record.id}")
        return record

    @staticmethod
    def run_blast(seq, record_id, db="nt", blast_program=None, top_hits=5):
        """Run BLAST search with auto-selection of program when possible."""
        outdir = "results"
        os.makedirs(outdir, exist_ok=True)

        # Accept Seq or string
        seq_str = str(seq).strip()

        # Auto choose program if not provided
        if not blast_program:
            # Quick heuristic: if sequence contains amino-acid-only letters -> protein
            s_type = SequenceAnalyzer.detect_sequence_type(seq_str)
            if s_type == "Protein":
                blast_program = "blastp"
                db = "nr" if db == "nt" else db
            else:
                blast_program = "blastn"
                db = "nt"

        print(f"🔬 Running BLAST for {record_id} with program {blast_program} against {db}")
        print("⏳ This may take several minutes...")

        try:
            result_handle = NCBIWWW.qblast(blast_program, db, seq_str)
            outpath = os.path.join(outdir, f"{record_id}_blast.xml")
            with open(outpath, "w") as f:
                f.write(result_handle.read())
            result_handle.close()
            print(f"✅ Results saved: {outpath}")

            with open(outpath) as handle:
                # NCBIXML.read expects a single record — catch if multiple
                try:
                    blast_record = NCBIXML.read(handle)
                except Exception:
                    handle.seek(0)
                    blast_records = list(NCBIXML.parse(handle))
                    blast_record = blast_records[0] if blast_records else None

                if not blast_record:
                    print("⚠️  No BLAST hits returned or parsing failed.")
                    return

                print(f"\n🎯 Top {top_hits} BLAST Hits:")
                print("="*60)
                for i, alignment in enumerate(blast_record.alignments[:top_hits], 1):
                    hsp = alignment.hsps[0]
                    print(f"\n{i}. {alignment.title.splitlines()[0][:100]}")
                    print(f"   Score: {hsp.score:.1f}, E-value: {hsp.expect:.2e}")
                    print(f"   Identities: {hsp.identities}/{hsp.align_length}")
        except Exception as e:
            print(f"❌ BLAST error: {e}")
            print("Tip: Check your internet connection, NCBI usage limits, and that sequence and program types match.")

# ============================================================
# VISUALIZATION
# ============================================================
class Visualizer:
    """Handles all visualization"""

    @staticmethod
    def plot_aa_composition(record):
        """Plot amino acid composition"""
        seq = str(record.seq).replace("*", "").replace("-", "")
        analyzer = ProteinAnalysis(seq)
        comp = analyzer.get_amino_acids_percent()

        df = pd.DataFrame.from_dict(comp, orient="index", columns=["Fraction"])
        df = df.sort_values("Fraction", ascending=True)

        plt.figure(figsize=(10, 6))
        df.plot(kind="barh", legend=False, ax=plt.gca())
        plt.title(f"Amino Acid Composition: {record.id}")
        plt.xlabel("Fraction")
        plt.ylabel("Amino Acid")
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_gc_content(record, window=100):
        """Plot GC content with safety for short sequences."""
        seq = record.seq
        L = len(seq)
        if L < 2 or window < 1:
            print("⚠️  Sequence too short for GC plot")
            return

        gc_values = []
        positions = []

        # Slide or chunk-based—use step = max(1, window//2) to show more points for small sequences
        step = max(1, window // 2)
        for i in range(0, max(1, L - window + 1), step):
            gc = gc_fraction(seq[i:i+window]) * 100
            gc_values.append(gc)
            positions.append(i + window//2)

        if not gc_values:
            print("⚠️  Not enough sequence to compute GC windows with given window size.")
            return

        mean_gc = sum(gc_values)/len(gc_values)
        plt.figure(figsize=(12, 4))
        plt.plot(positions, gc_values, linewidth=2)
        plt.axhline(y=mean_gc, linestyle='--', label=f'Mean: {mean_gc:.1f}%')
        plt.xlabel('Position (bp)')
        plt.ylabel('GC Content (%)')
        plt.title(f'GC Content: {record.id}')
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_structure_coords(atom_df):
        """Plot structure coordinates"""
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        axes[0].scatter(atom_df['X'], atom_df['Y'], s=1, alpha=0.5)
        axes[0].set_xlabel('X')
        axes[0].set_ylabel('Y')
        axes[0].set_title('XY Coordinates')

        axes[1].scatter(atom_df['X'], atom_df['Z'], s=1, alpha=0.5)
        axes[1].set_xlabel('X')
        axes[1].set_ylabel('Z')
        axes[1].set_title('XZ Coordinates')

        axes[2].scatter(atom_df['Y'], atom_df['Z'], s=1, alpha=0.5)
        axes[2].set_xlabel('Y')
        axes[2].set_ylabel('Z')
        axes[2].set_title('YZ Coordinates')

        plt.tight_layout()
        plt.show()

# ============================================================
# MAIN PIPELINE
# ============================================================
class BioinformaticsPipeline:
    """Main pipeline controller"""

    def __init__(self):
        self.current_files = {}

    def run(self):
        """Main execution loop"""
        print_header("BIOINFORMATICS ANALYSIS PIPELINE")

        while True:
            print("\n" + "="*60)
            print("MAIN MENU")
            print("="*60)
            print("1. Upload and analyze files")
            print("2. Fetch from NCBI database")
            print("3. View loaded files")
            print("4. Detailed analysis menu")
            print("5. Run BLAST search")
            print("6. Build phylogenetic tree")
            print("0. Exit")

            choice = input("\n➤ Enter choice: ").strip()

            if choice == "0":
                print("\n✅ Pipeline complete. Goodbye!")
                break
            elif choice == "1":
                self.upload_files()
            elif choice == "2":
                self.fetch_ncbi()
            elif choice == "3":
                self.view_files()
            elif choice == "4":
                self.analysis_menu()
            elif choice == "5":
                self.run_blast_search()
            elif choice == "6":
                self.build_phylogenetic_tree()
            else:
                print("❌ Invalid choice")

    def upload_files(self):
        """Handle file upload"""
        print("\n📤 Upload your files (FASTA, PDB, GenBank, Newick)")
        if COLAB_AVAILABLE:
            uploaded = files.upload()
        else:
            uploaded = _local_upload_prompt()
            if not uploaded:
                # Allow user to type filenames (comma separated) for local use
                fnames = input("Enter local filenames (comma-separated): ").strip().split(",")
                uploaded = {}
                for fn in fnames:
                    fn = fn.strip()
                    if not fn:
                        continue
                    try:
                        with open(fn, "rb") as fh:
                            uploaded[os.path.basename(fn)] = fh.read()
                    except Exception as e:
                        print(f"⚠️ Could not read {fn}: {e}")

        for fname, data in uploaded.items():
            print(f"\n{'='*60}")
            print(f"📁 Processing: {fname}")

            file_type = FileHandler.detect_file_type(fname, data)
            print(f"🧩 Detected type: {file_type}")

            self.current_files[fname] = {
                'data': data,
                'type': file_type
            }

            self.initial_analysis(fname, data, file_type)

    def initial_analysis(self, fname, data, file_type):
        """Perform initial analysis"""
        if file_type == "sequence":
            records = FileHandler.parse_sequence(data)
            print(f"✅ Found {len(records)} sequence(s)")

            seq_type = None
            for record in records:
                seq_type = SequenceAnalyzer.detect_sequence_type(record.seq)
                print(f"\n🧬 {record.id}")
                print(f"   Length: {len(record.seq)}")
                print(f"   Type: {seq_type}")

                if seq_type == "DNA":
                    try:
                        gc = gc_fraction(record.seq) * 100
                        print(f"   GC%: {gc:.2f}")
                    except Exception:
                        pass
                elif seq_type == "Protein":
                    stats = SequenceAnalyzer.analyze_protein_basic(record)
                    print(f"   MW: {stats['mw']:.0f} Da, pI: {stats['pI']:.2f}")

            self.current_files[fname]['records'] = records
            # Save last detected seq_type if multiple; this matches earlier script behavior
            self.current_files[fname]['seq_type'] = seq_type

        elif file_type == "protein_structure":
            try:
                structure, atom_df, fpath = StructureAnalyzer.parse_pdb_data(data, fname)
                print(f"✅ {len(atom_df)} atoms, {atom_df['Chain'].nunique()} chains")
                self.current_files[fname]['atoms'] = atom_df
                self.current_files[fname]['pdb_path'] = fpath
            except Exception as e:
                print(f"❌ Error parsing PDB: {e}")

        elif file_type == "genbank":
            try:
                record = FileHandler.parse_genbank(data)
                print(f"✅ {record.id}: {len(record.seq)} bp, {len(record.features)} features")
                self.current_files[fname]['record'] = record
                self.current_files[fname]['seq_type'] = 'DNA'
            except Exception as e:
                print(f"❌ Error parsing GenBank: {e}")

        elif file_type == "phylogenetic_tree":
            try:
                print("🌳 Phylogenetic tree detected")
                # Plot directly from provided newick text
                PhylogeneticAnalyzer.plot_tree(newick_text=data.decode(errors='ignore'))
                self.current_files[fname]['data'] = data
            except Exception as e:
                print(f"❌ Error handling tree file: {e}")

        else:
            print("⚠️ Unknown file type — stored for manual inspection")
            self.current_files[fname]['data'] = data

    def fetch_ncbi(self):
        """Fetch from NCBI"""
        print("\n🔍 NCBI Fetch")
        acc_ids = input("Accession IDs (comma-separated): ").strip().split(',')
        db = input("Database (nucleotide/protein): ").strip().lower()
        if db not in ("nucleotide", "protein"):
            print("⚠️ Invalid database, defaulting to nucleotide")
            db = "nucleotide"

        for acc_id in acc_ids:
            acc_id = acc_id.strip()
            if not acc_id:
                continue
            try:
                record = DatabaseTools.fetch_record(acc_id, db)
                fname = f"{acc_id}.gb" if db == 'nucleotide' else f"{acc_id}.fasta"
                self.current_files[fname] = {
                    'record': record,
                    'type': 'genbank' if db == 'nucleotide' else 'sequence',
                    'seq_type': 'DNA' if db == 'nucleotide' else 'Protein'
                }
            except Exception as e:
                print(f"❌ Error fetching {acc_id}: {e}")

    def view_files(self):
        """View loaded files"""
        if not self.current_files:
            print("\n⚠️  No files loaded")
            return

        print("\n📚 Loaded Files:")
        for i, (fname, info) in enumerate(self.current_files.items(), 1):
            t = info.get('type', 'unknown')
            print(f"{i}. {fname} ({t})")

    def analysis_menu(self):
        """Analysis submenu with context-aware options"""
        if not self.current_files:
            print("\n⚠️  No files loaded")
            return

        self.view_files()
        idx = input("\nSelect file number: ").strip()

        try:
            fname = list(self.current_files.keys())[int(idx) - 1]
            info = self.current_files[fname]
            file_type = info.get('type', 'unknown')

            print(f"\n{'='*60}")
            print(f"Analysis Options for: {fname}")
            print(f"Type: {file_type}")
            print("="*60)

            # Build context-aware menu
            menu_options = {}
            option_num = 1

            # Sequence file options
            if file_type == "sequence":
                menu_options[str(option_num)] = ("detailed_stats", "Detailed sequence statistics")
                option_num += 1

                if 'seq_type' in info and info['seq_type'] == "Protein":
                    menu_options[str(option_num)] = ("aa_comp", "Amino acid composition plot")
                    option_num += 1
                    menu_options[str(option_num)] = ("motifs", "Motif identification")
                    option_num += 1

                if 'seq_type' in info and info['seq_type'] in ["DNA", "RNA"]:
                    menu_options[str(option_num)] = ("translate", "Translate to protein")
                    option_num += 1
                    menu_options[str(option_num)] = ("gc_plot", "GC content plot")
                    option_num += 1

                if 'records' in info and len(info['records']) >= 2:
                    menu_options[str(option_num)] = ("align", "Sequence alignment")
                    option_num += 1

            # Structure file options
            elif file_type == "protein_structure":
                menu_options[str(option_num)] = ("visualize_3d", "3D structure visualization")
                option_num += 1
                menu_options[str(option_num)] = ("coord_plot", "Coordinate scatter plots")
                option_num += 1
                menu_options[str(option_num)] = ("atom_stats", "Atom statistics")
                option_num += 1

            # GenBank file options
            elif file_type == "genbank":
                menu_options[str(option_num)] = ("detailed_stats", "Detailed sequence statistics")
                option_num += 1
                menu_options[str(option_num)] = ("gc_plot", "GC content plot")
                option_num += 1
                menu_options[str(option_num)] = ("features", "Show features")
                option_num += 1
                menu_options[str(option_num)] = ("translate", "Extract and translate CDS")
                option_num += 1

            # Phylogenetic tree options
            elif file_type == "phylogenetic_tree":
                menu_options[str(option_num)] = ("plot_tree", "Plot phylogenetic tree")
                option_num += 1

            # Display menu
            for key, (action, desc) in menu_options.items():
                print(f"{key}. {desc}")
            print("0. Back")

            choice = input("\n➤ Choice: ").strip()

            if choice == "0":
                return

            if choice in menu_options:
                action = menu_options[choice][0]
                self.execute_analysis(action, fname, info)
            else:
                print("❌ Invalid choice")

        except (ValueError, IndexError):
            print("❌ Invalid selection")

    def execute_analysis(self, action, fname, info):
        """Execute specific analysis action"""

        if action == "detailed_stats":
            if 'records' in info:
                for r in info['records']:
                    SequenceAnalyzer.detailed_analysis(r)
            elif 'record' in info:
                SequenceAnalyzer.detailed_analysis(info['record'])

        elif action == "aa_comp":
            if 'records' in info and info['records']:
                Visualizer.plot_aa_composition(info['records'][0])
            else:
                print("❌ No protein record available for amino acid composition")

        elif action == "motifs":
            if 'records' in info:
                for r in info['records']:
                    SequenceAnalyzer.identify_motifs(r)
            else:
                print("❌ No sequences available for motif search")

        elif action == "translate":
            if 'records' in info:
                for r in info['records']:
                    proteins = SequenceAnalyzer.translate_sequence(r)
                    for desc, prot in proteins:
                        if prot:
                            print(f"\n{desc}:")
                            print(f"{str(prot)[:200]}..." if len(str(prot)) > 200 else str(prot))
            elif 'record' in info:
                proteins = SequenceAnalyzer.translate_sequence(info['record'])
                for desc, prot in proteins:
                    if prot:
                        print(f"\n{desc}:")
                        print(f"{str(prot)[:200]}..." if len(str(prot)) > 200 else str(prot))
            else:
                print("❌ No sequence data to translate")

        elif action == "gc_plot":
            if 'record' in info:
                Visualizer.plot_gc_content(info['record'])
            elif 'records' in info and info['records']:
                Visualizer.plot_gc_content(info['records'][0])
            else:
                print("❌ No sequence data for GC plot")

        elif action == "align":
            if 'records' in info:
                SequenceAnalyzer.align_sequences(info['records'])
            else:
                print("❌ No sequences to align")

        elif action == "visualize_3d":
            if 'data' in info:
                view = StructureAnalyzer.visualize_structure_py3dmol(info['data'])
                try:
                    view.show()
                except Exception:
                    print("⚠️ Could not display interactive view in this environment.")
            else:
                print("❌ No structure data available")

        elif action == "coord_plot":
            if 'atoms' in info:
                Visualizer.plot_structure_coords(info['atoms'])
            else:
                print("❌ No atomic coordinates available")

        elif action == "atom_stats":
            if 'atoms' in info:
                print("\n📊 Atom Statistics:")
                print("="*60)
                try:
                    print(info['atoms'].describe())
                except Exception as e:
                    print(f"⚠️ Could not describe atoms: {e}")
                print("\nChain distribution:")
                print(info['atoms']['Chain'].value_counts())
                print("\nResidue types (top 10):")
                print(info['atoms']['Residue'].value_counts().head(10))
            else:
                print("❌ No atom data available")

        elif action == "features":
            if 'record' in info:
                self.show_genbank_features(info['record'])
            else:
                print("❌ No GenBank record available")

        elif action == "plot_tree":
            if 'data' in info:
                PhylogeneticAnalyzer.plot_tree(newick_text=info['data'].decode(errors='ignore'))
            else:
                print("❌ No tree data available")

    def show_genbank_features(self, record):
        """Display GenBank features"""
        print(f"\n🏷️  Features in {record.id}:")
        print("="*60)

        feature_types = {}
        for f in record.features:
            feature_types[f.type] = feature_types.get(f.type, 0) + 1

        print("\n📊 Feature Summary:")
        for ftype, count in sorted(feature_types.items()):
            print(f"  {ftype}: {count}")

        print("\n📋 Detailed Features (first 20):")
        for i, f in enumerate(record.features[:20], 1):
            print(f"\n{i}. Type: {f.type}")
            print(f"   Location: {f.location}")
            if 'product' in f.qualifiers:
                print(f"   Product: {f.qualifiers['product'][0]}")
            if 'gene' in f.qualifiers:
                print(f"   Gene: {f.qualifiers['gene'][0]}")

    def build_phylogenetic_tree(self):
        """Build phylogenetic tree from sequences"""
        if not self.current_files:
            print("\n⚠️  No files loaded")
            return

        # Find files with multiple sequences
        suitable_files = []
        for fname, info in self.current_files.items():
            if 'records' in info and len(info['records']) >= 3:
                suitable_files.append((fname, len(info['records'])))

        if not suitable_files:
            print("\n⚠️  No suitable files found")
            print("Need at least 3 sequences in a single file to build a tree")
            return

        print("\n🌳 Build Phylogenetic Tree")
        print("="*60)
        print("Available files with multiple sequences:")
        for i, (fname, count) in enumerate(suitable_files, 1):
            print(f"{i}. {fname} ({count} sequences)")

        choice = input("\nSelect file number: ").strip()

        try:
            fname = suitable_files[int(choice) - 1][0]
            info = self.current_files[fname]
            records = info['records']

            print(f"\nSelected: {fname} with {len(records)} sequences")
            print("\nTree construction methods:")
            print("1. UPGMA (Unweighted Pair Group Method with Arithmetic Mean)")
            print("2. NJ (Neighbor Joining)")

            method_choice = input("\nSelect method (1/2): ").strip()
            method = "upgma" if method_choice == "1" else "nj"

            tree = PhylogeneticAnalyzer.build_tree_from_sequences(records, method)

            if tree:
                print("\n✅ Phylogenetic tree successfully created!")
                print("Tree file saved to: results/phylogenetic_tree.nwk")

        except (ValueError, IndexError):
            print("❌ Invalid selection")

    def run_blast_search(self):
        """Run BLAST"""
        if not self.current_files:
            print("\n⚠️  No sequences loaded")
            return

        self.view_files()
        idx = input("\nSelect file for BLAST: ").strip()

        try:
            fname = list(self.current_files.keys())[int(idx) - 1]
            info = self.current_files[fname]

            # Determine sequence to BLAST
            seq = None
            seq_id = None

            if 'records' in info:
                if len(info['records']) > 1:
                    print("\nMultiple sequences found:")
                    for i, r in enumerate(info['records'], 1):
                        print(f"{i}. {r.id}")
                    seq_choice = input("Select sequence number: ").strip()
                    record = info['records'][int(seq_choice) - 1]
                else:
                    record = info['records'][0]
                seq = record.seq
                seq_id = record.id
            elif 'record' in info:
                seq = info['record'].seq
                seq_id = info['record'].id

            if seq:
                print(f"\nBLAST Options:")
                print("1. Nucleotide BLAST (blastn) - for DNA/RNA")
                print("2. Protein BLAST (blastp) - for proteins")
                print("3. Auto (choose based on detected sequence type)")

                blast_choice = input("Select BLAST type (1/2/3): ").strip()
                if blast_choice == "1":
                    db = "nt"
                    prog = "blastn"
                elif blast_choice == "2":
                    db = "nr"
                    prog = "blastp"
                else:
                    db = "nt"
                    prog = None

                DatabaseTools.run_blast(seq, seq_id, db=db, blast_program=prog)
            else:
                print("❌ No sequence found")

        except (ValueError, IndexError):
            print("❌ Invalid selection")

# ============================================================
# RUN PIPELINE
# ============================================================
if __name__ == "__main__":
    pipeline = BioinformaticsPipeline()
    pipeline.run()
