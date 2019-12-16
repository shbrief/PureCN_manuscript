findMNPs = function(x) {
    # `mnps` is the index of merged nucleotides
    mnps = c()
    
    # MNPs with 3 nucleotides
    for (i in 1:(nrow(x)-2)) {
        if (abs(x$pos[i] - x$pos[i+1]) == 1 & abs(x$pos[i+1] - x$pos[i+2]) == 1) {  
            x$ref[i] = paste0(x$ref[i], x$ref[i+1], x$ref[i+2])
            x$mut[i] = paste0(x$mut[i], x$mut[i+1], x$mut[i+2])
            mnps = c(mnps, x$pos[i+1], x$pos[i+2]) 
        }
    }
    
    # MNPs with 2 nucleotides    
    for (i in 1:(nrow(x)-1)) {
        if (nchar(x$ref[i]) != 1 | x$pos[i] %in% mnps) {
            next
        } else if (abs(x$pos[i] - x$pos[i+1]) == 1) {
            x$ref[i] = paste0(x$ref[i], x$ref[i+1])
            x$mut[i] = paste0(x$mut[i], x$mut[i+1])
            mnps = c(mnps, x$pos[i+1])
        }
    }    

    # if there is MNP longer than 3 nucleotides    
    for (i in 1:(nrow(x)-1)) {
        if (nchar(x$ref[i]) != 1 & nchar(x$ref[i+1]) != 1) {
            mnps = NULL
            print(paste("Check MNP at row", i, "of sample", unique(x$sampleID)))
        }
    }
    
    # remove merged positions
    if (is.null(mnps)) {   
        x = x   
    } else {
        x = x[-which(x$pos %in% mnps),]   
    }
    
    return(x)
}