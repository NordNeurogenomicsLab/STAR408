#AmplifyDNA_2 function, which is a modified by me AmplifyDNA function from the DECIPHER package.
# modifications
#    * Primer hybridization over indels is dissalowed
#    * Primers cannot have any mismatches with the reference. This may not be optimal, in case a small mismatch results in a much shorter product, making it potentially favoured in a PCR reaction. This simplicifaction dramatically increses the function performance though. 

# For further consideration
# * Sequence dissambiguation could be removed

AmplifyDNA_2 <- function (primers, myDNAStringSet, maxProductSize, annealingTemp, 
                          P, ions = 0.2, includePrimers = TRUE, minEfficiency = 0.001, 
                          ...) 
{
  if (is.character(primers)) 
    primers <- DNAStringSet(primers)            # converting primers to BStringSet class
  if (!is(primers, "DNAStringSet"))               # checking for the object class
    stop("primers must be a DNAStringSet.")
  if (any(primers == ""))                         # checking if any of the primers is ""
    stop("primers cannot contain 0 character elements.")
  if (is.character(myDNAStringSet)) 
    myDNAStringSet <- DNAStringSet(myDNAStringSet) # converting DNA template to BStringSet class
  if (!is(myDNAStringSet, "DNAStringSet"))           # checking for the object class
    stop("myDNAStringSet must be a DNAStringSet.") 
  if (length(myDNAStringSet) == 0)                   # checking if DNA string has length 0     
    stop("myDNAStringSet is empty.")
  if (!is.numeric(ions)) 
    stop("ions must be a numeric.")
  if (ions < 0.01 || is.nan(ions)) 
    stop("Sodium equivilent concentration must be at least 0.01M.")
  if (!is.numeric(P))                               # P is the primer concentration  
    stop("P must be a numeric.")
  if (!(P > 0)) 
    stop("P must be greater than zero.")
  if (!is.numeric(annealingTemp)) 
    stop("annealingTemp must be a numeric.")
  if (!is.numeric(maxProductSize)) 
    stop("maxProductSize must be a numeric.")
  if (maxProductSize <= 0) 
    stop("maxProductSize must be greater than zero")
  if (!is.logical(includePrimers)) 
    stop("includePrimers must be a logical.")
  a <- vcountPattern("-", myDNAStringSet)              # checking if there are any gaps in the sequence
  if (any(a > 0)) 
    stop("Gap characters ('-') must be removed before amplification.")
  a <- vcountPattern("+", myDNAStringSet)
  if (any(a > 0)) 
    stop("Mask characters ('+') must be removed before amplification.")   # Checkign for masrking + character
  a <- vcountPattern(".", myDNAStringSet)
  if (any(a > 0)) 
    stop("Unknown characters ('.') must be removed before amplification.") # Checkign for unknown characters
  
  #Expand Ambiguities into All Permutations of a DNAStringSet, Important if IUPAC_CODE_MAP characters are present
  primers <- Disambiguate(primers)           
  l <- unlist(lapply(primers, length))
  l <- rep(seq_along(primers), l)
  primers <- unlist(primers)
  
  # I think this is for handling multiple template sequences? 
  w <- width(myDNAStringSet)
  w <- w[-length(w)]
  myDNAStringSet <- unlist(myDNAStringSet) #Flatten list-like objects from {BiocGenerics}
  if (length(w) > 0) 
    myDNAStringSet <- replaceAt(myDNAStringSet, IRanges(start = cumsum(w) + # ?? test w 2 seqs
                                                          1L, width = 0), paste(rep("-", maxProductSize), 
                                                                                collapse = ""))
  
  # Forward primer position - I made this more stringent
  f <- IRanges()
  fp <- integer()
  for (i in 1:length(primers)) {
    temp <- matchPattern(primers[[i]], myDNAStringSet, max.mismatch = 0, with.indels = FALSE)
    temp <- as(temp, "IRanges")
    f <- c(f, temp)
    fp <- c(fp, rep(i, length(temp)))
  }
  f <- Views(myDNAStringSet, start(f), end(f))
  if (length(f) == 0) 
    return(DNAStringSet())
  
  # Reverse primer position
  r <- IRanges()
  rp <- integer()
  for (i in 1:length(primers)) {
    temp <- matchPattern(reverseComplement(primers[[i]]), 
                         myDNAStringSet, max.mismatch = 0, with.indels = FALSE)
    temp <- as(temp, "IRanges")
    r <- c(r, temp)
    rp <- c(rp, rep(i, length(temp)))
  }
  r <- Views(myDNAStringSet, start(r), end(r))
  if (length(r) == 0) 
    return(DNAStringSet())
  
  ends <- end(r)
  o <- order(ends)
  r <- r[o]
  rp <- rp[o]
  
  # Forward primer efficiency consideration
  targets <- extractAt(myDNAStringSet, at = IRanges(start = start(f), 
                                                    end = end(f)))
  fe <- CalculateEfficiencyPCR(primers[fp], reverseComplement(targets), 
                               annealingTemp, P, ions)
  fw <- which(fe >= minEfficiency)
  if (length(fw) == 0) 
    return(DNAStringSet())
  
  # Reverse primer efficiency consideration
  targets <- extractAt(myDNAStringSet, at = IRanges(start = start(r), 
                                                    end = end(r)))
  re <- CalculateEfficiencyPCR(primers[rp], targets, annealingTemp, 
                               P, ions)
  rw <- which(re >= minEfficiency)
  if (length(rw) == 0) 
    return(DNAStringSet())
  
  # start/end positions for each primer, after efficiency consideration
  sf <- start(f)[fw]
  sr <- start(r)[rw]
  ef <- end(f)[fw]
  er <- end(r)[rw]
  
  
  #Placeholders?
  b <- e <- fs <- rs <- integer(1e+06)
  effs <- numeric(1e+06)
  c <- 0L
  
  for (i in 1:length(sf)) {
    w <- .Call("boundedMatches", sr, as.integer(sf[i] + 
                                                  1L), as.integer(sf[i] + maxProductSize - 1L), PACKAGE = "DECIPHER")
    if (length(w) > 0) {
      r <- (c + 1L):(c + length(w))
      c <- c + length(w)
      if (includePrimers) {
        b[r] <- ef[i] + 1L
        e[r] <- sr[w] - 1L
      }
      else {
        b[r] <- sf[i]
        e[r] <- er[w]
      }
      effs[r] <- sqrt(fe[fw[i]] * re[rw[w]])
      fs[r] <- fp[fw[i]]
      rs[r] <- rp[rw[w]]
    }
  }
  
  # This extracts the sequence strings. Full seq can be retrieved with as.character(v)
  if (c > 0) {
    length(b) <- length(e) <- length(effs) <- c
    o <- order(effs, decreasing = TRUE)
    if (includePrimers) {
      extra <- ifelse(b > e, b - e - 1L, 0)
      e <- ifelse(b > e, b - 1L, e)
      v <- Views(myDNAStringSet, start = b[o], end = e[o])
      v <- as(v, "DNAStringSet")
      v <- xscat(substring(primers[fs[o]], 1L, width(primers[fs[o]]) - 
                             extra[o]), v, reverseComplement(primers[rs[o]]))
      #names(v) <- paste(round(100 * effs[o], 1), "% (", 
      #    l[fs[o]], " x ", l[rs[o]], ")", sep = "")
      
      #I'm using "XStringViews" instead of "DNAStringSet" to preserve coordinates
      v_kc <- Views(myDNAStringSet, start = b[o] - width(primers[1]), end = e[o] + width(primers[2]))
      names(v_kc) <- paste(round(100 * effs[o], 1), "% (", 
                           l[fs[o]], " x ", l[rs[o]], ")", sep = "")
      
    }
    else {
      v <- Views(myDNAStringSet, start = b[o], end = e[o])
      #,names = paste(round(100 * effs[o], 1), "% (", 
      #      l[fs[o]], " x ", l[rs[o]], ")", 
      #      sep = ""))
      v <- as(v, "DNAStringSet")
      
      #I'm using "XStringViews" instead of "DNAStringSet" to preserve coordinates
      v_kc <- Views(myDNAStringSet, start = b[o] - width(primers[1]), end = e[o] + width(primers[2]))
      names(v_kc) <- paste(round(100 * effs[o], 1), "% (", 
                           l[fs[o]], " x ", l[rs[o]], ")", sep = "")
      
      
    }
    w <- which(width(v) > maxProductSize)
    if (length(w) > 0) 
      v <- v[-w]
    v_kv <- v_kc
    
  }
  else {
    #v <- DNAStringSet()
    v_kv <- DNAStringSet()
  }
  return(v_kv)
}
