set -e

#Prints a help message
function print_help() {
  echo "Usage: $0 [options] --in longReads.fasta --out result.fasta --type readsTechnology"
  echo ""
  echo "	Input:"
  echo "	longReads.fasta:               fasta file of long reads to correct, with one sequence per line."
  echo "	result.fasta:                  fasta file where to output the corrected long reads."
  echo "	readsTechnology:               Indicate whether the long reads are from PacBio (--type PB) or Oxford Nanopore (--type ONT)"
  echo ""
  echo "	Options:"
  echo "	--windowSize INT, -l INT:      Size of the windows to process. (default: 500)"
  echo "	--minSupport INT, -s INT:      Minimum support to consider a window for correction. (default: 4)"
  echo "	--maxSupport INT, -S INT:      Maximum number of sequences to include into a window. (default: 1,000)"
  echo "	--maxMSA INT, -M:              Maximum number of sequences to include into the MSA. (default: 150)"
  echo "	--merSize INT, -k INT:         k-mer size for chaining and polishing. (default: 9)"
  echo "	--solid INT, -f INT:           Minimum number of occurrences to consider a k-mer as solid during polishing. (default: 4)"
  echo "	--anchorSupport INT, -c INT:   Minimum number of sequences supporting (Ai) - (Ai+1) to keep the two anchors in the chaining. (default: 8)"
  echo "	--minAnchors INT, -a INT:      Minimum number of anchors in a window to allow consensus computation. (default: 10)"
  echo "	--windowOverlap INT, -o INT:   Overlap size between consecutive windows. (default: 50)"
  echo "	--nproc INT, -j INT:           Number of processes to run in parallel (default: number of cores)."
  echo "  --minimapIndex INT, -m INT:    Split minimap2 index every INT input bases (default: 500M)."
  echo "	--tmpdir STRING, -t STRING:    Path where to store the temporary overlaps file (default: working directory, as Alignments_dateTimeStamp.paf)."
  echo "	--help, -h:                    Print this help message."
  exit 1
}

#Set options to default values
reads=""
proof=""
nproc=$(nproc)
tmpdir="."
out=""
minSupport=4
maxSupport=1000
maxMSA=150
windowSize=500
merSize=9
commonKMers=8
minAnchors=2
solid=4
windowOverlap=50
minimapOption=""
minimapMemory="500M"

#Print help if no argument specified
if [[ "$1" == "" ]]; then
  print_help
fi

#Otions handling
while [[ "$1" != "" ]]; do
  case "$1" in
  "--help" | "-h")
    print_help
    ;;
  "--in")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      reads="$2"
      shift 2
      ;;
    esac
    ;;
  "--out")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      out="$2"
      shift 2
      ;;
    esac
    ;;
  "--type")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *) if [[ "$2" == "PB" ]]; then
      minimapOptions="PB"
      shift 2
    elif [[ "$2" == "ONT" ]]; then
      minimapOptions="ONT"
      shift 2
    else
      echo "Error: $1 must be either PB or ONT"
      exit 1
    fi ;;
    esac
    ;;
  "--minSupport" | "-s")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      minSupport="$2"
      shift 2
      ;;
    esac
    ;;
  "--maxSupport" | "-S")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      maxSupport="$2"
      shift 2
      ;;
    esac
    ;;
  "--maxMSA" | "-M")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      maxMSA="$2"
      shift 2
      ;;
    esac
    ;;
  "--windowSize" | "-l")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      windowSize="$2"
      shift 2
      ;;
    esac
    ;;
  "--merSize" | "-k")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      merSize="$2"
      shift 2
      ;;
    esac
    ;;
  "--anchorSupport" | "-c")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      commonKMers="$2"
      shift 2
      ;;
    esac
    ;;
  "--minAnchors" | "-a")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      minAnchors="$2"
      shift 2
      ;;
    esac
    ;;
  "--solid" | "-f")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      solid="$2"
      shift 2
      ;;
    esac
    ;;
  "--windowOverlap" | "-o")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      windowOverlap="$2"
      shift 2
      ;;
    esac
    ;;
  "--nproc" | "-j")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      nproc="$2"
      shift 2
      ;;
    esac
    ;;
  "--minimapIndex" | "-m")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      minimapMemory="$2"
      shift 2
      ;;
    esac
    ;;
  "--tmpdir" | "-t")
    case "$2" in
    "")
      echo "Error: $1 expects an argument"
      exit 1
      ;;
    *)
      tmpdir="$2"
      shift 2
      ;;
    esac
    ;;
  --)
    shift
    break
    ;;
  *)
    echo "Error: invalid option \"$1\""
    exit 1
    ;;
  esac
done

#Exit if no input or no output files have been specified
if [[ $reads == "" ]]; then
  echo "Error: --in must be specified"
  exit 1
fi
if [[ $out == "" ]]; then
  echo "Error: --out must be specified"
  exit 1
fi
if [[ $minimapOptions == "" ]]; then
  echo "Error: --type must be specified"
  exit 1
fi
mkdir -p $tmpdir

#Remove the output file if it already exists
if [[ -f $out ]]; then
  rm $out
fi
