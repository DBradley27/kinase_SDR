## Code written by Dr. Omar Wagih

#' MATCH-TM Implemntation + Log weights score
#' @author omarwagih
#' @references pubmed:12824369
require(parallel)


#' Namespaces
AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
DNA = c('A', 'C', 'G', 'T')
RNA = c('A', 'C', 'G', 'U')

#' Priors for amino acids and dna
AA_PRIORS  =   list(human = c(A=0.070, R=0.056, N=0.036, D=0.048,C=0.023,
                              Q=0.047, E=0.071, G=0.066, H=0.026, I=0.044,
                              L=0.100, K=0.058, M=0.021, F=0.037, P=0.063,
                              S=0.083, T=0.053, W=0.012, Y=0.027, V=0.060),
                    
                    yeast = c(A=0.055, R=0.045, N=0.061, D=0.058, C=0.013,
                              Q=0.039, E=0.064, G=0.05, H=0.022, I=0.066,
                              L=0.096, K=0.073, M=0.021, F=0.045, P=0.044,
                              S=0.091, T=0.059, W=0.01, Y=0.034, V=0.056))

DNA_PRIORS =   list(human = c(A=0.293, C=0.207, G=0.200, T=0.300),
                    yeast = c(A=0.313, C=0.187, G=0.171, T=0.329),
                    ecoli = c(A=0.247, C=0.260, G=0.257, T=0.236))


#' Compute information content (equation from http://en.wikipedia.org/wiki/Sequence_logo)
#' Different from that used in MATCH algorithm
#' 
#' @param pwm position weight matrix containing relation freqiences
#' @param N number of letters (20 for AA, 4 for DNA/RNA)
#' @param Nseqs number of sequences in the alignment
computeBits <- function(pwm, N=4, Nseqs){
  H_i = - apply(pwm, 2, function(col) sum(col * log2(col), na.rm=T))
  e_n = (1/logb(2)) * (N-1)/(2*Nseqs) 
  R_i = log2(N) - (H_i  + e_n)
  # Set any negatives to 0
  R_i = pmax(R_i, 0)
  return(R_i)
}

getPriors <- function(priors='eq', seq.type='DNA', N=4){
  namespace = get(seq.type)
  # If priors is equiprobable
  if(priors[1] == 'eq'){
    priors = rep(1/N, N)
    names(priors) = namespace
  }else if(is.character(priors[1])){
    organism = priors[1]
    priors.list = get(sprintf('%s_PRIORS', seq.type))
    if(! organism %in% names(priors.list)) stop('Could not find organism in priors!')
    priors = priors.list[[organism]]
  }
  
  # Ensure namespace and priors have same length
  if(length(priors) != N) stop('Priors and namespace must have same length')
  
  return(priors)
}

findNamespace <- function(seq.type, sp){
  
  if(seq.type == 'auto'){
    dat = setdiff(intersect(sp, AA), c(DNA, RNA))
    if(length(dat) > 0) return('AA')
    if('U' %in% sp) return('RNA')
    return('DNA')
  }
  return(seq.type)
  
}


# Get first consecutive x highest conservation 
coreIndicies <- function(ic, core=5){
  # Get possible start indices
  starts = 1:(length(ic)-core+1)
  # Get indicies of length core
  cons_index = lapply(starts, function(s)  s:(s+core-1))
  # Compute conservation for these streches 
  cons = sapply(cons_index, function(ind) sum(ic[ind]) )
  # Find best set
  p = cons_index[[which.max(cons)]]
  
  p
}



letterMatrix <- function(input){
  # Ensure kmers are the same length characters 
  seq.len = sapply(input, nchar)
  num.pos = seq.len[1]
  if(! all(seq.len == num.pos)) stop('Unequal length of sequences')
  
  # Construct matrix of letters
  split = unlist( sapply(input, function(seq){strsplit(seq, '')}) )
  
  m = t( matrix(split, seq.len, length(split)/num.pos) )
  m
}


#' Construct position weight matrix
#' 
#' Makes a position weight matrix given aligned sequences.
#'
#' @param input Aligned sequences all of the same length OR PFM as matrix, where rows are letters and columns are positions
#' @param pseudocount Pseudocount factor. Final pseudocount is pseudocount / number of of letters in namespace
#' @param relative.freq TRUE if each column should be divided by the sum
#' @param seq.type seq.type of sequences 'AA', 'DNA'. Default is autodetect.
#' @param log.bg If true, relative frequencies will be converted to weights using -log of freq/bgfreq
#' @param priors Named character vector containing priors of letters. If set to organism name, it'll be picked up automatically: "human", "yeast". If priors is set to "eq", then equiprobable priors are used (i.e. prior for each letter is 1/total number of letters)
#' @keywords pwm construct
#' @export
#' @examples
#' # No examples
makePWM <- function(input, pseudocount=1, relative.freq=T, seq.type='auto', log.bg=F, priors='eq', core=5){
  if(!all(is.character(input)) & !is.matrix(input) ) stop('Input must be char vector or a matrix!')
  no.pfm = all(is.character(input))
  
  # We dont need pseudocounts if we're not logging!
  if(!log.bg) pseudocount = 0
  
  # Ensure seq.type is correct
  if(!seq.type %in% c('AA', 'DNA', 'auto')) stop('seq.type must be AA, DNA, or auto')
  
  if(no.pfm){
    num.pos = nchar(input[1])
    m = letterMatrix(input)
    split = as.character(m)
    nseqs = length(input)
  }else{
    num.pos = ncol(input)
    split = rownames(input)
    nseqs = max(apply(input, 2, sum))
    m = input
  }
  
  # Chose correct namespace
  seq.type = findNamespace(seq.type, split)
  namespace = get(seq.type)
  
  # Number of letters in ns
  N = length(namespace)
  
  # Get priors
  my.priors = getPriors(priors, seq.type, N)
  
  # Match priors to namespace 
  bg.prob = my.priors[match(namespace, names(my.priors))]
  # Construct PWM
  pwm.matrix = apply(m, 2, function(pos.data){
    
    t = pos.data
    # Get frequencies 
    if(no.pfm) t = table(pos.data)
    
    # Match to aa
    ind = match(namespace, names(t))
    
    # Create column
    col = t[ind]
    col[is.na(col)] = 0
    names(col) = namespace
    
    # Do pseudocounts
    #if(log.bg) col = col + (pseudocount / N) 
    #if(log.bg) col = col + bg.prob
    if(log.bg) col = col + (sqrt(nseqs) * bg.prob)
    
    # Do relative frequencies
    if(relative.freq) col = col / sum(col)
    
    col
  })
  
  # Set attributes
  attr(pwm.matrix, 'pseudocount') = pseudocount
  
  # Information vector as computed by MATCH paper
  match.ic = apply(pwm.matrix, 2, function(col) sum(col * logb(N * col), na.rm=T))
  attr(pwm.matrix, 'match.ic') = match.ic
  attr(pwm.matrix, 'bits') = computeBits(pwm.matrix, N, nseqs)
  
  # If log.bg, compute weights as log(fij/bi)
  if(log.bg) pwm.matrix = apply(pwm.matrix, 2, function(col) log2(col / bg.prob))
  #if(log.bg) pwm.matrix = apply(pwm.matrix, 2, function(col) log2(col / 0.25))
  
  # Assign AA names to rows/pos col
  rownames(pwm.matrix) = namespace
  colnames(pwm.matrix) = 1:num.pos
  attr(pwm.matrix, 'log.bg') = log.bg
  attr(pwm.matrix, 'seq.type') = seq.type
  
  core = pmin(core, ncol(pwm.matrix))
  if(!is.na(core))
    attr(pwm.matrix, 'core') = coreIndicies(match.ic, core)
  
  
  return(pwm.matrix)
}

#' Get weight/probability for each amino acid in a sequence 
#' 
#' Gets weight/probability for the amino acid at each position of the sequence
#' as an array.
#'
#' @param seqs One or more sequences to be processed
#' @param pwm Position weight matrix
#'  
#' @keywords pwm mss match tfbs
#' @examples
#' # No Examples
scoreArray <- function(seqs, pwm){
  # Split sequence
  sp = strsplit(seqs, '')
  
  seq.lens = sapply(sp, length)
  seq.len = seq.lens[1]
  if(any(seq.lens != seq.len)) stop('Input sequences must be same length')
  
  # Iterate through sequences
  dat = lapply(sp, function(seq){
    # Match sequence to the PWM
    mat = matrix(c(match(seq, rownames(pwm)), 1:seq.len), seq.len, 2)
    prob.vector = pwm[ mat ]
    prob.vector
  })
  names(dat) = seqs
  return(dat)
}

#' Given a position weight matrix, find the best matching sequence
#' 
#' Finds the amino acid at each position of the PWM with the highest occurence.
#' Used in matrix similarity score calculation.  
#'
#' @param pwm Position weight matrix
#' @keywords pwm best
#' @examples
#' # No Examples
bestSequence <- function(pwm){
  b = rownames(pwm)[apply(pwm, 2, which.max)]
  return(paste(b, collapse=''))
}

#' Given a position weight matrix, find the worst matching sequence
#' 
#' Finds the amino acid at each position of the PWM with the lowest occurence.
#' Used in matrix similarity score calculation.  
#'
#' @param pwm Position weight matrix
#' @keywords pwm worst
#' @examples
#' # No Examples
worstSequence <- function(pwm){
  w = rownames(pwm)[apply(pwm, 2, which.min)]
  return(paste(w, collapse=''))
}

#' Compute classical log score 
#' 
#' Computes sum of weights of a pwm for a givet set of sequences.
#'
#' @param seqs Sequences to be scored
#' @param pwm Position weight matrix
#' @param na.rm Remove NA scores?
#' @param ignore.central If true, ignores central residue weight
#'  
#' @keywords pwm log
#' @examples
#' # No Examples
logScore <- function(seqs, pwm, na.rm=F, ignore.central=F, relative=T){
  
  seqs = c(bestSequence(pwm), worstSequence(pwm), seqs)
  if(ignore.central){
    # Central residue index
    central.ind = ceiling(ncol(pwm)/2)
    scores = sapply(scoreArray(seqs, pwm), function(sa) sum(sa[-central.ind], na.rm=T))
  }else{
    scores = sapply(scoreArray(seqs, pwm), function(sa) sum(sa, na.rm=T))
  }
  
  max.score = scores[1]
  min.score = scores[2]
  scores = scores[-(1:2)]
  if(relative){
    # (site_score - min_matrix_score) / (max_matrix_score - min_matrix_score)
    scores = (scores - min.score)/(max.score - min.score)
  }
  
  # Remove NA if requested
  if(na.rm) scores = scores[!is.na(scores)]
  return(scores)
}


## a useful function: rev() for strings
strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

revComp <- function(seqs, rev=T, comp=T){
  if(nchar(seqs[1]) == 1) warning('Reverse compliment wont work on split sequence!')
  cseqs = seqs
  if(comp) cseqs = chartr("ATGCatgc","TACGtacg",seqs)
  if(rev) cseqs = strReverse(cseqs)
  cseqs
}

#' Compute matrix similarity score as described in MATCH algorithm
#' 
#' Computes score of a PWM with a k-mer.
#' There are two methods, 
#' "mss" gives a score from 0-1, as described in [PMID: 12824369]
#' "log" gives the sum of log weights
#'
#' @param seqs Sequences to be scored
#' @param pwm Position weight matrix
#' @param kinase.pwm TRUE if PWM is that of a kinase (special case). If you have DNA sequences, this will be automatically set to FALSE.
#' @param na.rm Remove NA scores?
#' @param ignore.central Ignore central residue, for things like PTMs where you want to ignore the weight of the modified site
#' @param both.strands If you are scoring DNA sequences, the method will score both the sequence and the reverse compliment (reverse strand). If you want the forward strand only (i.e. the input), set to FALSE. Default is TRUE.
#'  
#' @keywords pwm mss match tfbs log pfm
#' @examples
#' # No Examples
matchScore <- function(seqs, pwm, kinase.pwm=T, na.rm=F, ignore.central=F, both.strands=F){
  
  if(attr(pwm, 'seq.type') != 'AA') kinase.pwm = F
  
  # Must ignore central for kinase pwms
  if(kinase.pwm) ignore.central = T
  
  # Central residue index
  central.ind = ceiling(ncol(pwm)/2)
  
  # Cannot have log weights for mss
  if(attr(pwm, 'log.bg')) 
      stop('Cannot have log weights for MSS, please reconstruct the pwm with log.bg=FALSE')
    
    
  # Best/worst sequence match
  oa = scoreArray(bestSequence(pwm), pwm)[[1]]
  wa = scoreArray(worstSequence(pwm), pwm)[[1]]
    
  # MATCH information content
  IC = attr(pwm, 'match.ic')
    
  # Get vector of frequencies 
  score.arr = scoreArray(seqs, pwm)
    
  mssHelper <- function(sa, wa, oa, IC){
    curr.score  = sum( IC * sa, na.rm=T ) 
    opt.score   = sum( IC * oa, na.rm=T )
    worst.score = sum( IC * wa, na.rm=T )
    score.final = ( (curr.score - worst.score) / (opt.score - worst.score) )
    score.final
  }
    
  # Get core indicies
  core.ind = attr(pwm, 'core')
  if(ignore.central) core.ind = setdiff(core.ind, central.ind)
    
  scores = lapply(score.arr, function(sa){
    na = is.na(sa)
    na[central.ind] = ignore.central
    keep = !na
    
    mss.score = mssHelper(sa[keep], wa[keep], oa[keep], IC[keep])
    css.score = mssHelper(sa[core.ind], wa[core.ind], oa[core.ind], IC[core.ind])
    
    c(mss=mss.score, css=css.score)
  })
    
  scores.mss = sapply(scores, function(i) i[[1]])
  scores.css = sapply(scores, function(i) i[[2]])
  
  
  # Only score sequences which have a central residue S/T or Y depending on the PWM
  central.res = '*'
  if(kinase.pwm){
    kinase.seq.type = names(which.max(pwm[,ceiling(ncol(pwm)/2)]))
    kinase.seq.type = ifelse(grepl('S|T', kinase.seq.type), 'S|T', 'Y')
    central.res = kinase.seq.type
  }
  
  # Set scores to NA for kinase stuff
  if(central.res != '*'){
    keep = grepl(central.res, substr(seqs, central.ind, central.ind))
    scores.mss[!keep] = NA
    scores.css[!keep] = NA
  }
  
  # Remove NA if requested
  if(na.rm){
    scores.mss = scores.mss[!is.na(scores.mss)]
    scores.css = scores.css[!is.na(scores.css)]
  }
  
  scores = list(mss=scores.mss, css=scores.css)
  # Do reverse strands
  if(attr(pwm, 'seq.type') == 'DNA' & both.strands){
    rc = revComp(seqs)
    li.scores = list('+' = scores, 
                     '-' = matchScore(rc, pwm, kinase.pwm, na.rm, ignore.central, both.strands = F))
    return(li.scores)
  }
  
  return(scores)
}

