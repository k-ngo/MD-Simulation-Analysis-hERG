proc get_backbone_dihedral_indices {resid resname segname mol} {
  
  set ressel [atomselect $mol "segname $segname and resid $resid and name CA"]
  set residue_index [$ressel get residue]
  $ressel delete
  
  set prev_residue_index [expr $residue_index-1]
  set next_residue_index [expr $residue_index+1]
  
  #Nitrogen
  set atomsel [atomselect $mol "segname $segname and residue $residue_index and name \"N\""]
  set Ni [$atomsel get index]
  $atomsel delete

  #CA
  set atomsel [atomselect $mol "segname $segname and residue $residue_index and name \"CA\""]
  set CAi [$atomsel get index]
  $atomsel delete
  
  #C(O)
  set atomsel [atomselect $mol "segname $segname and residue $residue_index and name \"C\""]
  set Ci [$atomsel get index]
  $atomsel delete
  
  

  #if prevous residue exists
  
  if {$prev_residue_index>=0} then {
    
    #previous C(O)
    set atomsel [atomselect $mol "segname $segname and residue $prev_residue_index and name \"C\""]
    set pCi [$atomsel get index]
    $atomsel delete

    
    lappend phi $pCi $Ni $CAi $Ci
  } else  {
    set phi {}
  }


  #next N
  set atomsel [atomselect $mol "segname $segname and residue $next_residue_index and name \"N\""]
  set nNi [$atomsel get index]
  $atomsel delete

  #next CA
  set atomsel [atomselect $mol "segname $segname and residue $next_residue_index and name \"CA\""]
  set nCAi [$atomsel get index]
  $atomsel delete
  

  #if next residue exists
  if {$nNi>0} then {
    lappend psi $Ni $CAi $Ci $nNi
    lappend omega $CAi $Ci $nNi $nCAi
  } else  {
    set psi {}
    set omega {}
  }  
  
  lappend result $phi $psi $omega 
  return $result  
  

}

proc get_sidechain_inidces {resid resname segname mol} {

  lappend result $resid $resname $segname
  lappend result [get_backbone_dihedral_indices $resid $resname $segname $mol]
  
  #puts $result
  switch $resname {
    GLY {lappend result [get_GLY_sidechain_indicies $resid $resname $segname $mol]}
    ALA {lappend result [get_ALA_sidechain_indicies $resid $resname $segname $mol]}
    SER {lappend result [get_SER_sidechain_indicies $resid $resname $segname $mol]}
    CYS {lappend result [get_CYS_sidechain_indicies $resid $resname $segname $mol]}
    VAL {lappend result [get_VAL_sidechain_indicies $resid $resname $segname $mol]}
    THR {lappend result [get_THR_sidechain_indicies $resid $resname $segname $mol]}
    ILE {lappend result [get_ILE_sidechain_indicies $resid $resname $segname $mol]}
    PRO {lappend result [get_PRO_sidechain_indicies $resid $resname $segname $mol]}
    MET {lappend result [get_MET_sidechain_indicies $resid $resname $segname $mol]}
    ASP {lappend result [get_ASP_sidechain_indicies $resid $resname $segname $mol]}
    ASN {lappend result [get_ASN_sidechain_indicies $resid $resname $segname $mol]}
    LEU {lappend result [get_LEU_sidechain_indicies $resid $resname $segname $mol]}
    LYS {lappend result [get_LYS_sidechain_indicies $resid $resname $segname $mol]}
    GLU {lappend result [get_GLU_sidechain_indicies $resid $resname $segname $mol]}
    GLN {lappend result [get_GLN_sidechain_indicies $resid $resname $segname $mol]}
    ARG {lappend result [get_ARG_sidechain_indicies $resid $resname $segname $mol]}
    
    HIS {lappend result [get_HIS_sidechain_indicies $resid $resname $segname $mol]}
    HSE {lappend result [get_HIS_sidechain_indicies $resid $resname $segname $mol]}
    HSP {lappend result [get_HIS_sidechain_indicies $resid $resname $segname $mol]}
    HSD {lappend result [get_HIS_sidechain_indicies $resid $resname $segname $mol]}
    
    PHE {lappend result [get_PHE_sidechain_indicies $resid $resname $segname $mol]}
    TYR {lappend result [get_TYR_sidechain_indicies $resid $resname $segname $mol]}
    TRP {lappend result [get_TRP_sidechain_indicies $resid $resname $segname $mol]}
    UNK {lappend result [get_UNK_sidechain_indicies $resid $resname $segname $mol]}
    
    CMT {lappend result [get_CMT_sidechain_indicies $resid $resname $segname $mol]}
    CMTS {lappend result [get_CMT_sidechain_indicies $resid $resname $segname $mol]}
    } 
    

    return $result
}

proc get_all_sidechain_indices {resinfos segname mol} {
  set result {}
  foreach {resid resname} $resinfos {
     #puts  "$resid $resname $segname $mol"
     lappend result [get_sidechain_inidces $resid $resname $segname $mol] 
  }
  return $result
}




proc print_all_sidechain_dihedrals {indices frame_num output_file} {
  #puts $indices
  foreach {resid resname segname backbone_dihedrals dihedrals} $indices {
    #puts "| $resid $resname $segname $backbone_dihedrals $dihedrals|"
    
    
    set phi_indices [lindex $backbone_dihedrals 0]
    set psi_indices [lindex $backbone_dihedrals 1]
    #set omega_inidces [lindex backbone_dihedrals 2/]
    
    #first phi doesent exist
    if {[llength $phi_indices]} then {              
      set phi [measure dihed $phi_indices]
    } else {
      set phi 999
    }
    
    #lst psi doesent exist
    if {[llength $psi_indices]} then {
      set psi [measure dihed $psi_indices]
    } else {
      set psi 999
    }
    
    set outstring [format "%4.d%4.4s%9.d%6.1f%6.1f" $resid $resname $frame_num $phi $psi]
    #if not empty list
    
    if {[llength $dihedrals]} {
      foreach dihedral $dihedrals {
        set val [measure dihed $dihedral]
        set val_str [format "%6.1f" $val]        
        append outstring $val_str
      } ;#foreach dihedral      
    } ;# if not empty list
    puts $output_file $outstring
  } ;#foreach residue 
}

puts "calculate dihedrals 1.0"
puts "proc print_sidechain_dihedrals {segname residues {mol top} {output_file_name \"stdout\"} {first_frame 0} {last_frame -1} {stride 1} {print_progress 0}}"

proc print_sidechain_dihedrals {segname residues {mol top} {output_file_name "stdout"} {first_frame 0} {last_frame -1} {stride 1} {print_progress 0}} {
  
  #get residue list
  set residues [atomselect $mol "segname $segname and $residues and name CA"]
  set resinfos [join [$residues get {resid resname}]]
  #puts $resinfos
  
  if {($print_progress>0) && ($output_file_name != "stdout")} then {
       set counter [expr [llength $resinfos]/2]
       puts "Getting indices for $counter residues."  
     }
  #get the atom indices of the dihedral (torsion) angles
  set indices [join [get_all_sidechain_indices $resinfos $segname $mol]]
  
  #puts $indices
  if {($print_progress>0) && ($output_file_name != "stdout")} then {
       set counter [expr [llength $resinfos]/2]
       puts "Done."  
     }


  
  #set output
  if {$output_file_name == "stdout"} then {
    set output_file stdout
  } else { 
    set output_file [open $output_file_name w]}
  
  
  #get current mol frame
  set old_frame [molinfo $mol get frame] 
  #don't update display
  display update off
  
  #default last_frame is the number of all frames
  if {$last_frame == -1} then {
    set last_frame [expr [molinfo $mol get numframes] -1]}
  
  set counter 0
  #iterate along through all frames
  for {set frame $first_frame} {$frame <= $last_frame} {incr frame $stride} {
     incr counter
     if {($print_progress>0) && ($counter>=$print_progress) && ($output_file_name != "stdout")} then {
       puts "Working on frame $frame."
       set counter 0
     }
      
     molinfo $mol set frame $frame
     print_all_sidechain_dihedrals $indices $frame $output_file
  }
  
 
  #restore frame
  molinfo $mol set frame $old_frame 
  #restore update
  display update on
  
  if {$output_file_name != "stdout"} then {
  close $output_file}
  $residues delete 
}



