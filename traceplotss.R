
#create trace plots from alex's code 
createTracePlots <- function(trace, model,genome,numMixtures,samples,mixture.labels,samples.percent.keep=1)
{
  for (i in 1:numMixtures)
  {
    plot(trace, what = "Mutation", mixture = i)
    plot(trace, what = "Selection", mixture = i)
    
    plot(model, genome, samples = samples*samples.percent.keep, mixture = i,main = mixture.labels[i])
  }
  plot(trace, what="AcceptanceRatio")
}


if (with.phi && !is.null(obs.phi))
{
  genome <- initializeGenomeObject(file=input,codon_table=codon.table,match.expression.by.id = FALSE,observed.expression.file=obs.phi)
} else{
  if (dev)
  {
    genome <- initializeGenomeObject(file=input,codon_table=codon.table)
  } else {
    genome <- initializeGenomeObject(file=input)
  }
}

