# --- Modified VMD Tcl script to generate an intermolecular distance histogram ---

# Declare that we want to use the global variable 'sel'
global sel

# 2. Define the cutoff (in Å) to search for contacts.
set cutoff 5.0

# 3. Measure contacts: exclude pairs from the same residue.
set pairs [measure contacts $cutoff $sel "not same residue"]

# 4. Define histogram parameters.
set binwidth 0.05
set numBins [expr {int($cutoff / $binwidth) + 1}]
# Initialize a list to hold histogram counts (all bins start at 0).
set hist {}
for {set i 0} {$i < $numBins} {incr i} {
    lappend hist 0
}

# 5. Retrieve coordinates for all selected atoms.
set coords [$sel get {x y z}]

# 6. Loop over each contact pair, compute distance, and update the histogram.
foreach pair $pairs {
    lassign $pair i j
    # Get the coordinates for atoms i and j.
    set r_i [lindex $coords $i]
    set r_j [lindex $coords $j]
    # Compute the Euclidean distance.
    set d [veclength [vecsub $r_i $r_j]]
    # Determine which bin the distance falls into.
    set binIndex [expr {int($d / $binwidth)}]
    if {$binIndex < $numBins} {
         # Increment the count in that bin.
         set current [lindex $hist $binIndex]
         set newVal [expr {$current + 1}]
         set hist [lreplace $hist $binIndex $binIndex $newVal]
    }
}

# 7. Write the histogram data to a file.
set fp [open "distance_histogram.dat" "w"]
puts $fp "# Lower_Angstrom  Upper_Angstrom  Count"
for {set i 0} {$i < $numBins} {incr i} {
    set lower [expr {$i * $binwidth}]
    set upper [expr {($i + 1) * $binwidth}]
    set count [lindex $hist $i]
    puts $fp "$lower $upper $count"
}
close $fp

puts "Histogram data written to distance_histogram.dat"

