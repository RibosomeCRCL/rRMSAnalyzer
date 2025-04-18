# Check if two samples share the same positions
#   sample_1: First sample to check against sample_2
#   sample_2 : Second sample to check against sample_1
#   sample_1_name : name of sample 1
#   sample_2_name : name of sample 2
#  Returns A boolean indicating if the two samples are identical.
.check_sample_positions <- function(sample_1, sample_2, sample_1_name, sample_2_name) {
  sample_size <- length(sample_1[, "named_position"])
  
  if (length(sample_2[, "named_position"]) == sample_size) { # On compare si sample_2 a le même nombre de positions que sample_1
    if (sum(sample_1["named_position"] == sample_2["named_position"]) == sample_size) # compte le nombre de positions identiques 
      return(TRUE) # Si toutes les positions sont identiques (sum == sample_size) retourne TRUE
    else {
      errored_positions_file_1 <- setdiff(sample_1[["named_position"]], sample_2[["named_position"]]) #Trouve les éléments de sample_1 qui ne sont pas dans sample_2
      errored_positions_file_2 <- setdiff(sample_2[["named_position"]], sample_1[["named_position"]]) #Trouve les éléments de sample_2 qui ne sont pas dans sample_1
      paste("[WARNING] Two samples have differents positions! \n Mismatch at ", sample_1_name, ": ",paste(errored_positions_file_1,collapse = ", "),"\n Mismatch at ",sample_2_name, ": ",paste(errored_positions_file_2,collapse = ", "))
      return(FALSE)
    }
    
  }
  else { # sinon on retourne un message d'erreur
    warning(paste("[WARNING]", sample_1_name, "and", sample_2_name ," have different size!"))
    return(FALSE)
  }
  
}