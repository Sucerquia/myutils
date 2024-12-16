# Load the molecule in VMD and display it in CPK representation, then save the image as a PNG file.

f { $argc < 2 } {
    puts "Usage: vmd -e create_mol_png.tcl -args input_file output_file"
    exit 1
}

# read input and output files
set input_file [lindex $argv 0]
set output_file [lindex $argv 1]

# Load the xyz file
mol new $input_file type xyz

# Set up CPK representation
mol representation CPK 1.0 0.3 20 20
mol color Name
mol selection all
mol material Opaque
mol addrep top
color Name C gray

# Set rendering parameters
display projection orthographic
display depthcue off
axes location off

color Display Background white

# Render the scene and save it as a PNG
render snapshot $output_file

# Exit VMD
exit