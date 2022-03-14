# Get number of frames, get output files
set nf [molinfo top get numframes]
set outmaxz [open ${outputname_max_z} w]
set outcomz [open ${outputname_com_z} w]
set outminz [open ${outputname_min_z} w]

# Obtain number of drug molecules
set nDrugs [expr {[[atomselect top "segname $drug" frame 0] num] / [[atomselect top "segname $drug and resid 1" frame 0] num]}]

# Loop through each drug
for {set i 1} {$i <= $nDrugs} {incr i} {
    set maxz_list {}
    set comz_list {}
    set minz_list {}

    # Loop through each frame
    for {set f 0} {$f < $nf} {incr f} {

        set drug_mol [atomselect top "segname $drug and resid $i and same residue as sqrt(x^2 + y^2) < 10 and z > $bottomPore and z < $S4" frame $f]
        
        if {[$drug_mol num] != 0} {
            # Get max z for drug molecule in current frame
            lappend maxz_list [tcl::mathfunc::max {*}[$drug_mol get z]]
            # Get z of center of mass
            lappend comz_list [lindex [measure center $drug_mol weight mass] 2]
            # Get min z
            lappend minz_list [tcl::mathfunc::min {*}[$drug_mol get z]]
        } else {
            # Get max z for drug molecule in current frame
            lappend maxz_list NaN
            # Get z of center of mass
            lappend comz_list NaN
            # Get min z
            lappend minz_list NaN
        }
    }

    puts $outmaxz "$i $maxz_list"
    puts $outcomz "$i $comz_list"
    puts $outminz "$i $minz_list"
}

close $outmaxz
close $outcomz
close $outminz