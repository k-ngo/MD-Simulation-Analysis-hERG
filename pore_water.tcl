# get number of frames, open a dat file
set nf [molinfo top get numframes]
set out [open ${outputname} w]    

for {set f 0} {$f < $nf} {incr f} {
	# make the atom selection, divide to get water
	
	# Inside Filter, z range = S4 to S1
	set watselFilter [atomselect top "water and same residue as sqrt(x^2 + y^2) < 2 and z > $S4 and z < $S1" frame $f]
    set watnumFilter [$watselFilter num]
    set watcountFilter [expr {$watnumFilter / 3}]

	# Behind Filter, range = water in range of sequence ALY (helix behind filter) and sequence GFG (filter)
	# PROA
    set watselPROA [atomselect top "water and same residue as within 3.5 of (segname PROA and sequence ALY) and same residue as within 3.5 of (segname PROA and sequence GFG)" frame $f]
    set watnumPROA [$watselPROA num]
    set watcountPROA [expr {$watnumPROA / 3}]
	# PROB
    set watselPROB [atomselect top "water and same residue as within 3.5 of (segname PROB and sequence ALY) and same residue as within 3.5 of (segname PROB and sequence GFG)" frame $f]
    set watnumPROB [$watselPROB num]
    set watcountPROB [expr {$watnumPROB / 3}]
	# PROC
    set watselPROC [atomselect top "water and same residue as within 3.5 of (segname PROC and sequence ALY) and same residue as within 3.5 of (segname PROC and sequence GFG)" frame $f]
    set watnumPROC [$watselPROC num]
    set watcountPROC [expr {$watnumPROC / 3}]
	# PROD
    set watselPROD [atomselect top "water and same residue as within 3.5 of (segname PROD and sequence ALY) and same residue as within 3.5 of (segname PROD and sequence GFG)" frame $f]
    set watnumPROD [$watselPROD num]
    set watcountPROD [expr {$watnumPROD / 3}]

	# Pore, z range = S4 to (S4 - 21 Angstrom)
    set watselPore [atomselect top "water and same residue as sqrt(x^2 + y^2) < 14 and z > $bottomPore and z < $S4" frame $f]
    set watnumPore [$watselPore num]
    set watcountPore [expr {$watnumPore / 3}]
    
    # and write to file
    puts $out "$f,$watcountFilter,$watcountPROA,$watcountPROB,$watcountPROC,$watcountPROD,$watcountPore"}

close $out