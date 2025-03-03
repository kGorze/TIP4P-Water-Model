# Define an atom selection for water. Adjust the selection if your water residue is named differently.
set sel [atomselect top "resname SOL"]

# Define a cutoff distance (in angstroms) that indicates a problematic close contact.
set cutoff 1.0

# Measure contacts: returns a list of pairs of indices for atoms that are closer than cutoff.
set contactPairs [measure contacts $cutoff $sel]

# Print how many contacts were found.
puts "Found [llength $contactPairs] contacts within $cutoff Å:"

# Loop over each pair and print them.
foreach pair $contactPairs {
    puts "$pair"
}

# (Optional) Highlight these atoms:
# You could, for example, change the color of the atoms that are involved in the contacts.
# Here we loop over the contact pairs and set their beta value to 100 for easy visualization.
foreach pair $contactPairs {
    lassign $pair a1 a2
    $sel frame 0
    # Get current beta values (if any) and set them to 100 for atoms a1 and a2.
    $sel frame 0; $sel set beta 0  ;# Clear any existing beta values.
    $sel frame 0; $sel set beta 100 [list $a1 $a2]
}

