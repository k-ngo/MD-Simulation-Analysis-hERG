# get number of frames, open a dat file
set nf [molinfo top get numframes] 

for {set f 0} {$f < $nf} {incr f} {

	# Get location of drug binding residues
	set drug_tyr [atomselect top "sequence YAS and resname TYR" frame $f]
	set drug_phe [atomselect top "sequence IFGNVS and resname PHE" frame $f]
	set drug_ser [atomselect top "sequence IFGNVS and resname SER" frame $f]

    # Output to file center of mass
	puts $drugBindingPoreResLocation "$f,[lindex [measure center $drug_tyr weight mass] 2],[lindex [measure center $drug_phe weight mass] 2],[lindex [measure center $drug_ser weight mass] 2]"}

close $drugBindingPoreResLocation