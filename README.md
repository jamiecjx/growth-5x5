To replicate results, the following commands should be run.
The number of threads running can be adjusted: the below commands are specified for 64 threads.

Lemma 1 is contained in the commands
```
julia -t 64 lemma1.jl 215//100 4
julia -t 64 lemma1.jl 220//100 4
```
Which output two `.jld2` files consisting of a single variable which is a `Vector{MMatrix{3, 3, Interval{Float64}}}`, with lengths 1882 and 19511 respectively.
Approximate runtime is 2 hours with 64 threads for the 2.15 case.

Lemma 1 is contained in the commands

```
julia -t 64 lemma2.jl 215//100 220//100 lemma1_2.15_4_matlist.jld2
julia -t 64 lemma2.jl 220//100 210//100 lemma1_2.2_4_matlist.jld2
```
which output two `.csv` files which consist of pairs `(index, value)`. The resulting file should always have `value = 0` which completes the proof.
Approximate runtime is 7 hours with 64 threads for the 2.15 case.
