package require pbctools

source prot_center.tcl

pbc wrap -center com -centersel "segid PROC and helix and backbone" -compound residue -all
pbc wrap -center com -centersel "protein and helix and backbone" -compound residue -all

pbc join res -first 0 -last 0 -sel "protein"

# pbc wrap -center com -centersel "protein and helix and backbone" -compound residue -all
pbc unwrap -sel "protein" -all

source ./prot_center.tcl