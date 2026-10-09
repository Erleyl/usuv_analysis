#!/usr/bin/env bash
#
# Enhanced USUV Lab Sequence Curation Script
# ==========================================
# 
# This script provides comprehensive curation for local lab sequences including:
# - Metadata validation and standardization
# - Sequence quality assessment
# - Data completeness scoring
# - Integration with existing pipeline
#
# Usage: bash curate_lab_data.sh [OPTIONS]

set -euo pipefail

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
export RUN_DATE=$(date +%Y%m%d)

# Directories
BASE_DIR="$PROJECT_ROOT/usuv_data"
LAB_DATA_DIR="$BASE_DIR/lab_sequences"
CURATED_DIR="$BASE_DIR/curated"
LOG_DIR="$BASE_DIR/logs"

# Create directories
mkdir -p "$LAB_DATA_DIR" "$CURATED_DIR" "$LOG_DIR"

# Logging function
log() {
    printf "\033[1;32m$(date +"[%Y-%m-%d %H:%M:%S]") [CURATE] %s\033[0m\n" "$*" | tee -a "$LOG_DIR/curation_$RUN_DATE.log"
}

error() {
    printf "\033[1;31m$(date +"[%Y-%m-%d %H:%M:%S]") [ERROR] %s\033[0m\n" "$*" | tee -a "$LOG_DIR/curation_$RUN_DATE.log"
}

warning() {
    printf "\033[1;33m$(date +"[%Y-%m-%d %H:%M:%S]") [WARNING] %s\033[0m\n" "$*" | tee -a "$LOG_DIR/curation_$RUN_DATE.log"
}

usage() {
    cat << EOF
=========================================================
   USUV LAB SEQUENCE CURATION PIPELINE
=========================================================

Usage: bash curate_lab_data.sh [OPTIONS]

OPTIONS:
  --input-metadata FILE   Input metadata CSV file (required)
  --input-sequences FILE  Input sequences FASTA file (required)
  --config FILE          Configuration file (JSON, optional)
  --output-dir DIR       Output directory (default: curated/)
  --dry-run             Check dependencies and validate inputs
  --validate-only       Run validation without generating outputs
  --force               Overwrite existing output files
  --help                Show this help message

EXAMPLES:
  # Basic curation
  bash curate_lab_data.sh \\
    --input-metadata lab_metadata.csv \\
    --input-sequences lab_sequences.fasta

  # With custom configuration
  bash curate_lab_data.sh \\
    --input-metadata lab_metadata.csv \\
    --input-sequences lab_sequences.fasta \\
    --config custom_config.json \\
    --output-dir results/

  # Dry run to check setup
  bash curate_lab_data.sh --dry-run

=========================================================
EOF
}

check_dependencies() {
    log "Checking dependencies..."
    
    local missing_deps=()
    
    # Check for required tools
    for tool in python3 pip; do
        if ! command -v "$tool" &> /dev/null; then
            missing_deps+=("$tool")
        fi
    done
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        error "Missing dependencies: ${missing_deps[*]}"
        error "Please install missing tools and try again"
        return 1
    fi
    
    # Check Python packages
    log "Checking Python packages..."
    python3 -c "
import sys
missing = []
try:
    import pandas
except ImportError:
    missing.append('pandas')
try:
    import numpy
except ImportError:
    missing.append('numpy')
try:
    import Bio
except ImportError:
    missing.append('biopython')

if missing:
    print('Missing Python packages:', ', '.join(missing))
    print('Install with: pip install', ' '.join(missing))
    sys.exit(1)
else:
    print('All Python packages available')
"
    
    log "Dependencies check completed ✓"
    return 0
}

validate_inputs() {
    local metadata_file="$1"
    local sequences_file="$2"
    
    log "Validating input files..."
    
    # Check if files exist
    if [[ ! -f "$metadata_file" ]]; then
        error "Metadata file not found: $metadata_file"
        return 1
    fi
    
    if [[ ! -f "$sequences_file" ]]; then
        error "Sequences file not found: $sequences_file"
        return 1
    fi
    
    # Check file formats
    log "Checking metadata file format..."
    if ! python3 -c "import pandas as pd; pd.read_csv('$metadata_file', nrows=1)" &>/dev/null; then
        error "Invalid CSV format in metadata file"
        return 1
    fi
    
    log "Checking sequences file format..."
    if ! python3 -c "
from Bio import SeqIO
try:
    next(SeqIO.parse('$sequences_file', 'fasta'))
    print('FASTA format valid')
except:
    print('Invalid FASTA format')
    exit(1)
" &>/dev/null; then
        error "Invalid FASTA format in sequences file"
        return 1
    fi
    
    log "Input validation completed ✓"
    return 0
}

prepare_workspace() {
    local output_dir="$1"
    local force="$2"
    
    log "Preparing workspace..."
    
    if [[ -d "$output_dir" && "$force" != "true" ]]; then
        warning "Output directory exists: $output_dir"
        read -p "Continue and overwrite? (y/N): " -r
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log "Aborted by user"
            exit 0
        fi
    fi
    
    mkdir -p "$output_dir"
    
    # Copy curation script to workspace
    local curation_script="$output_dir/curation_pipeline.py"
    if [[ ! -f "$curation_script" ]]; then
        if [[ -f "$SCRIPT_DIR/enhanced_curation_pipeline.py" ]]; then
            cp "$SCRIPT_DIR/enhanced_curation_pipeline.py" "$curation_script"
        else
            error "Curation pipeline script not found"
            return 1
        fi
    fi
    
    log "Workspace prepared ✓"
}

generate_summary_stats() {
    local metadata_file="$1"
    local output_dir="$2"
    
    log "Generating summary statistics..."
    
    python3 << EOF
import pandas as pd
import json
from pathlib import Path

# Read metadata
df = pd.read_csv('$metadata_file')

# Calculate summary statistics
stats = {
    'total_records': len(df),
    'unique_species': df['Species_English'].nunique() if 'Species_English' in df.columns else 0,
    'countries': df['Country'].value_counts().to_dict() if 'Country' in df.columns else {},
    'date_range': {
        'start': str(df['DeathOrSampleDate'].min()) if 'DeathOrSampleDate' in df.columns else 'N/A',
        'end': str(df['DeathOrSampleDate'].max()) if 'DeathOrSampleDate' in df.columns else 'N/A'
    },
    'ct_values': {
        'mean': float(df['CT'].mean()) if 'CT' in df.columns and df['CT'].notna().any() else 'N/A',
        'min': float(df['CT'].min()) if 'CT' in df.columns and df['CT'].notna().any() else 'N/A',
        'max': float(df['CT'].max()) if 'CT' in df.columns and df['CT'].notna().any() else 'N/A'
    }
}

# Save to file
output_file = Path('$output_dir') / 'pre_curation_stats.json'
with open(output_file, 'w') as f:
    json.dump(stats, f, indent=2)

print(f"Summary statistics saved to: {output_file}")
print(f"Total records: {stats['total_records']}")
print(f"Unique species: {stats['unique_species']}")
print(f"Countries: {len(stats['countries'])}")
EOF
}

run_curation() {
    local metadata_file="$1"
    local sequences_file="$2"
    local output_dir="$3"
    local config_file="$4"
    
    log "Running curation pipeline..."
    
    local curation_script="$output_dir/curation_pipeline.py"
    local cmd="python3 '$curation_script' --metadata '$metadata_file' --sequences '$sequences_file' --output '$output_dir'"
    
    if [[ -n "$config_file" ]]; then
        cmd+=" --config '$config_file'"
    fi
    
    log "Executing: $cmd"
    eval "$cmd"
    
    if [[ $? -eq 0 ]]; then
        log "Curation completed successfully ✓"
    else
        error "Curation pipeline failed"
        return 1
    fi
}

create_pipeline_integration() {
    local output_dir="$1"
    
    log "Creating pipeline integration scripts..."
    
    # Create integration script for main pipeline
    cat > "$output_dir/integrate_curated_data.sh" << 'EOF'
#!/usr/bin/env bash
#
# Integration script for curated lab data
# Usage: bash integrate_curated_data.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Source directories
CURATED_DIR="$SCRIPT_DIR"
MAIN_DATA_DIR="$PROJECT_ROOT/usuv_data"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log "Integrating curated lab sequences into main pipeline..."

# Copy NCBI-compatible data to main pipeline (these files match NCBI format)
if [[ -f "$CURATED_DIR/ncbi_compatible_sequences.fasta" ]]; then
    cp "$CURATED_DIR/ncbi_compatible_sequences.fasta" "$MAIN_DATA_DIR/lab_sequences_curated.fasta"
    log "✓ NCBI-compatible sequences copied (ready for merge with NCBI data)"
fi

if [[ -f "$CURATED_DIR/ncbi_compatible_metadata.csv" ]]; then
    cp "$CURATED_DIR/ncbi_compatible_metadata.csv" "$MAIN_DATA_DIR/lab_metadata_curated.csv"
    log "✓ NCBI-compatible metadata copied"
fi

# Also preserve original format files
if [[ -f "$CURATED_DIR/high_quality_sequences.fasta" ]]; then
    cp "$CURATED_DIR/high_quality_sequences.fasta" "$MAIN_DATA_DIR/lab_sequences_original.fasta"
    log "✓ Original format sequences preserved"
fi

if [[ -f "$CURATED_DIR/high_quality_metadata.csv" ]]; then
    cp "$CURATED_DIR/high_quality_metadata.csv" "$MAIN_DATA_DIR/lab_metadata_original.csv"
    log "✓ Original format metadata preserved"
fi

# Create a summary file for the merge step
cat > "$MAIN_DATA_DIR/lab_data_summary.txt" << SUMMARY
Lab Data Integration Summary
===========================
Generated on: $(date)

Files created for pipeline integration:
- lab_sequences_curated.fasta    : NCBI-compatible sequences (merge-ready)
- lab_metadata_curated.csv       : NCBI-compatible metadata
- lab_sequences_original.fasta   : Original format sequences (backup)
- lab_metadata_original.csv      : Original format metadata (backup)

Format compatibility:
- NCBI sequences: taxa_name format (e.g., AS202400001.turdus_merula.nld.2024)
- NCBI metadata: Compatible with usuv_master_metadata.csv structure
- Ready for merge step: YES

Next step: Run 'bash main.sh --merge' to combine with NCBI data
SUMMARY

log "Integration completed. Files are ready for merge step."
log "Lab sequences now have NCBI-compatible names like: SeqID.host.country.year"
log "Ready to merge with NCBI data. Run: bash main.sh --merge"
EOF
    
    chmod +x "$output_dir/integrate_curated_data.sh"
    
    # Create quality report script
    cat > "$output_dir/generate_quality_report.sh" << 'EOF'
#!/usr/bin/env bash
#
# Generate quality comparison report
# Usage: bash generate_quality_report.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 << 'PYTHON_EOF'
import pandas as pd
import json
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

script_dir = Path('$SCRIPT_DIR')

# Read curated metadata if available
curated_file = script_dir / 'curated_metadata.csv'
if curated_file.exists():
    df = pd.read_csv(curated_file)
    
    # Generate quality score distribution plot
    if 'QC_SCORE' in df.columns:
        plt.figure(figsize=(10, 6))
        plt.hist(df['QC_SCORE'], bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Quality Score (%)')
        plt.ylabel('Number of Sequences')
        plt.title('Distribution of Sequence Quality Scores')
        plt.savefig(script_dir / 'quality_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Quality distribution plot saved to: {script_dir}/quality_distribution.png")
    
    # Generate species composition plot
    if 'Species_English' in df.columns:
        species_counts = df['Species_English'].value_counts().head(10)
        plt.figure(figsize=(12, 8))
        species_counts.plot(kind='bar')
        plt.xlabel('Species')
        plt.ylabel('Number of Sequences')
        plt.title('Top 10 Species in Dataset')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(script_dir / 'species_composition.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Species composition plot saved to: {script_dir}/species_composition.png")
PYTHON_EOF
EOF
    
    chmod +x "$output_dir/generate_quality_report.sh"
    
    log "Integration scripts created ✓"
}

# Parse command line arguments
METADATA_FILE=""
SEQUENCES_FILE=""
CONFIG_FILE=""
OUTPUT_DIR="$CURATED_DIR"
DRY_RUN=false
VALIDATE_ONLY=false
FORCE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --input-metadata)
            METADATA_FILE="$2"
            shift 2
            ;;
        --input-sequences)
            SEQUENCES_FILE="$2"
            shift 2
            ;;
        --config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --validate-only)
            VALIDATE_ONLY=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Main execution
main() {
    log "Starting USUV lab sequence curation pipeline..."
    
    # Check dependencies
    if ! check_dependencies; then
        exit 1
    fi
    
    # Dry run mode
    if [[ "$DRY_RUN" == true ]]; then
        log "Dry run completed successfully ✓"
        exit 0
    fi
    
    # Validate required arguments
    if [[ -z "$METADATA_FILE" || -z "$SEQUENCES_FILE" ]]; then
        error "Both --input-metadata and --input-sequences are required"
        usage
        exit 1
    fi
    
    # Validate inputs
    if ! validate_inputs "$METADATA_FILE" "$SEQUENCES_FILE"; then
        exit 1
    fi
    
    # Validate only mode
    if [[ "$VALIDATE_ONLY" == true ]]; then
        log "Validation completed successfully ✓"
        exit 0
    fi
    
    # Prepare workspace
    if ! prepare_workspace "$OUTPUT_DIR" "$FORCE"; then
        exit 1
    fi
    
    # Generate pre-curation statistics
    generate_summary_stats "$METADATA_FILE" "$OUTPUT_DIR"
    
    # Run curation pipeline
    if ! run_curation "$METADATA_FILE" "$SEQUENCES_FILE" "$OUTPUT_DIR" "$CONFIG_FILE"; then
        exit 1
    fi
    
    # Create integration scripts
    create_pipeline_integration "$OUTPUT_DIR"
    
    log "=========================================="
    log "Curation pipeline completed successfully!"
    log "=========================================="
    log "Output directory: $OUTPUT_DIR"
    log "Next steps:"
    log "1. Review the curation report: $OUTPUT_DIR/curation_report_*.html"
    log "2. Run integration: bash $OUTPUT_DIR/integrate_curated_data.sh"
    log "3. Continue with main pipeline: bash main.sh --merge"
}

# Run main function
main "$@"
