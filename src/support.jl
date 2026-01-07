"""
	get_chunks(n, tasks_per_thread=2)

Returns chunk size and an iterator of chunks for parallelizing work over `n` items.
Each chunk is a range of indices. The chunk size is determined by dividing `n` by
`tasks_per_thread * Threads.nthreads()` to ensure good load balancing.
"""
function get_chunks(n::Int, tasks_per_thread::Int=2)
	chunk_size = max(1, n ÷ (tasks_per_thread * Threads.nthreads()))
	data_chunks = Iterators.partition(1:n, chunk_size)
	return chunk_size, data_chunks
end

"""
	getindX(supportX,val) 

Returns index of `val` within `supportX`. Requires `supportX` to be sorted. 
"""
function getindX(supportX,val) 
	searchsorted(supportX,val)[1]
end

"""
	estimateΠ(X,supportX)

Returns estimate of transition matrix `Π`, using data `X` and support `supportX`.
"""
function estimateΠ(X,supportX)
	nperiod = size(X,1) 
	nfirm = size(X,2) 

	Π = zeros(Int,length(supportX),length(supportX))
	Π0 = zeros(Int,length(supportX))

	for t in 2:nperiod, f in 1:nfirm
			_ix = getindX(supportX,X[t-1,f])
			Π[_ix,getindX(supportX,X[t,f])] += 1 
			Π0[_ix]+=1
	end

	return Π ./ Π0 
end



