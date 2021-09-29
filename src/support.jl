"""
	getrange(n)

For a thread, returns respective range within 1:`n` to iterate over. 
"""
@inline function getrange(n)
	tid = Threads.threadid()
	nt = Threads.nthreads()
	d , r = divrem(n, nt)
	from = (tid - 1) * d + min(r, tid - 1) + 1
	to = from + d - 1 + (tid ≤ r ? 1 : 0)
	from:to
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



