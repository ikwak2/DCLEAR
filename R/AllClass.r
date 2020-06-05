
setOldClass('phyDat')
setOldClass('phylo')
setOldClass('igraph')

setClass(
	'lineage_tree',
	representation(
		x = 'phyDat',
		graph = 'igraph',
		outcome_prob = 'numeric',
		alphabets = 'character',
		division = 'numeric',
		n_samples = 'numeric',
		n_targets = 'numeric',
		deletion = 'logical',
		dropout_prob = 'numeric'
	)
)

setClass(
	'kmer_summary',
	representation(
		df = 'data.frame',
    alphabets = 'character',
    kmers = 'character',
    outcome_prob = 'numeric',
    sequence_length = 'numeric',
    division = 'numeric',
    reps = 'numeric',
    max_distance = 'numeric',
    mutation_prob = 'numeric',
    k = 'numeric',
		dropout_prob = 'numeric'
	)
)
