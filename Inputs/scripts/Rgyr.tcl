# VMD script to calculate radius of gyration for multiple chains
# and save results to a file

# Define the output file
set outfile [open "radius_of_gyration.dat" w]

# Write header to the file
puts $outfile "Chain\tRadius_of_Gyration(A)"

# Function to calculate radius of gyration
proc calc_rgyr {selection} {
    # Get the center of mass
    set com [measure center $selection weight mass]
    
    # Get coordinates and masses for all atoms in the selection
    set coords [$selection get {x y z}]
    set masses [$selection get mass]
    
    # Calculate radius of gyration
    set total_mass 0
    set rgyr_sum 0
    
    # Loop through each atom in the selection
    foreach coord $coords mass $masses {
        set total_mass [expr $total_mass + $mass]
        
        # Calculate square distance from center of mass
        set dx [expr [lindex $coord 0] - [lindex $com 0]]
        set dy [expr [lindex $coord 1] - [lindex $com 1]]
        set dz [expr [lindex $coord 2] - [lindex $com 2]]
        set dist_sq [expr $dx*$dx + $dy*$dy + $dz*$dz]
        
        # Add to sum
        set rgyr_sum [expr $rgyr_sum + $mass * $dist_sq]
    }
    
    # Calculate and return radius of gyration
    if {$total_mass > 0} {
        return [expr sqrt($rgyr_sum / $total_mass)]
    } else {
        return 0
    }
}

# List of all chain IDs
set chains {A B C D E F G H I J K L M N O P Q R S T U V W X Y Z a b c d e f g h i j k l m n o p q r s t u v w x}

# Loop through each chain
foreach chain $chains {
    # Create a selection for the current chain
    set sel [atomselect top "chain $chain"]
    
    # Check if the selection has any atoms
    if {[$sel num] > 0} {
        # Calculate radius of gyration
        set rgyr [calc_rgyr $sel]
        
        # Write to file with 4 decimal places
        puts $outfile "$chain\t[format "%.4f" $rgyr]"
        
        # Also print to console
        puts "Chain $chain: Radius of Gyration = [format "%.4f" $rgyr] Å"
    } else {
        puts "Chain $chain has no atoms"
        puts $outfile "$chain\t-"
    }
    
    # Delete the selection
    $sel delete
}

# Close the output file
close $outfile

puts "Radius of gyration calculations complete. Results saved to radius_of_gyration.dat"
