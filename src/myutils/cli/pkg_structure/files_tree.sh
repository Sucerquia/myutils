#!/bin/bash


print_help() {
echo "
Function that creates a mermaid map of a directory structure.

Usage: myutils files_tree <path to dir to be mapped> [ <ignoring patters> ... ]

The output is the string that you have to locate in an rst file in order to
render your map. In principle, the output can be used in markdown files
changing the corresponding keywords.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== Costumer set up ========================================================
verbose=''
while getopts 'vh' flag;
do
    case "${flag}" in

      v) verbose='true' ;;
      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

# Function to recursively generate Mermaid nodes
generate_tree() {
    local path="$1"
    local parent="$2"
    local indent="$3"
    local node_id="node$((++NODE_COUNTER))"

    # Add the current directory to the Mermaid graph
    if [[ "$parent" == "root" ]]
    then
      # TODO: replace myutils with any arbitrary name and the condition related with the indentation
      echo "${indent}${node_id}[\"myutils\"]"
    else
      echo "${indent}${parent} --> ${node_id}[\"$(basename "$path")\"]"
    fi
    # Add a click action for the node to link to the directory path
    if [[ "$path" == "." ]]
    then
       file=myutils.html
    else
      file=${path//\.\//}
      file=myutils/$file
      file=${file//\//\.}.html
    fi

    link=modules/$file
    echo "${indent}click $node_id \"$link\" _self"

    # Iterate through subdirectories only
    for child in "$path"/*; do
        [ -d "$child" ] || continue # Skip files
        # Skip directories matching ignore patterns
        for pattern in "${IGNORE_PATTERNS[@]}"; do
            if [[ "$(basename "$child")" == $pattern ]]; then
                continue 2
            fi
        done
        generate_tree "$child" "$node_id" "$indent  "
    done
}

# Check if the user provided a directory path
if [ $# -lt 1 ]; then
    echo "Usage: $0 <directory-path> [ignore-pattern1 ignore-pattern2 ...]"
    exit 1
fi

DIR_PATH="$1"
shift
IGNORE_PATTERNS=("$@") # Remaining arguments are ignore patterns
NODE_COUNTER=0

# Start the RST file content
echo ".. mermaid::"
echo "   :align: center"
echo ""
echo "   graph TD"

# Generate the tree
generate_tree "$DIR_PATH" "root" "   "
