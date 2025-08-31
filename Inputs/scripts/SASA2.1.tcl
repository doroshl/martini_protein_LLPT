# VMD script to calculate Solvent-Accessible Surface Area (SASA) for multiple chains
# and save results to a file

# Define the output file
set outfile [open "chain_sasa.dat" w]

# Write header to the file
puts $outfile "Chain\tSASA(Å²)"

# Define calculation parameters
set probe_radius 2.1
set samples_per_atom 500

# Define chain list properly
set chains {}
lappend chains A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
lappend chains a b c d e f g h i j k l m n o p q r s t u v w x

# Loop through each chain
foreach chain $chains {
    # Create a selection for the current chain
    set sel [atomselect top "chain $chain"]
    
    # Check if the selection has any atoms
    if {[$sel num] > 0} {
        # Calculate SASA
        set sasa [measure sasa $probe_radius $sel -samples $samples_per_atom]
        
        # Write to file with 2 decimal places
        puts $outfile "$chain\t[format "%.2f" $sasa]"
        
        # Also print to console
        puts "Chain $chain: SASA = [format "%.2f" $sasa] Å²"
    } else {
        puts "Chain $chain has no atoms"
        puts $outfile "$chain\t-"
    }
    
    # Delete the selection
    $sel delete
}

# Calculate total SASA for comparison
set all_sel [atomselect top "all"]
set total_sasa [measure sasa $probe_radius $all_sel -samples $samples_per_atom]
puts "Total SASA (all chains combined): [format "%.2f" $total_sasa] Å²"
puts $outfile "Total\t[format "%.2f" $total_sasa]"
$all_sel delete

# Close the output file
close $outfile

puts "SASA calculations complete. Results saved to chain_sasa.dat"
