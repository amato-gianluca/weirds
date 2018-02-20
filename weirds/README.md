# Primitive Weirds and Abundant Numbers

Counts and lists of primitive weird numbers. With the exception of the `experiments.txt` file, all weird numbers are provided in both a sequence of factors and a linear index sequence, i.e., an index sequence where `None` means that the last factor has been repeated. The `*` aside a weird means that its abundance is *small*, namely, less than the last prime factor.


* `experiments.txt`: a file with some random experiments on weird numbers in the form of a screen log. Please note that since experiments have been conducted with different versions of our software, the format of the output is not consistent.
* `sfweirds-list-by-bigomega.txt`: a screenlog of commands which produce all square-free weird numbers for a given number of factors (from 2 to 5).
* `weirds-list-by-bigomega.txt`: a screenlog of commands which produce all weird numbers for a given number of factors, counted with multiplicity (from 2 to 5).
* `weirds-nosfo-list-by-bigomega.txt`: a screenlog of commands which produce all weird numbers with a non square-free odd part, for a given number of factors, counted with multiplicity (from 3 to 7).