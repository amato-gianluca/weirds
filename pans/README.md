# Primitive Weirds and Abundant Numbers

Counts and lists of primitive abundant numbers. File with extension `.xz` are compressed in the [Xz format](https://en.wikipedia.org/wiki/Xz).

1. `count-opan-bigomega<n>.txt`: number of odd PANs with *n* prime factors, counted with multiplicities. Each line of the file has the form `[p1, ..., pn] c min max found` meaning that there are `m` odd PANs whose factorization start with `p1 ... pn`. `min` and `max` are the minimum and maximum values respectively of these PANs, and `found` is `True` when there is at least an abundant number, not necessarily primitive abundant, starting with [p1,...., pn].
2. `count-pan-bigomega<n>.txt`: like 1, but for all PANs, both odd and even.
3. `count-sfopan-bigomega<n>.txt`: similar to 1, but for square-free odd PANs. Each line has the form `[p1, ..., pn] c`, where `c` is the number of odd square-free PANs starting with `p1 ... pn`.
4. `count-sfpan-bigomega<n>.txt`: like 3, but for square-free PANs, both odd and even.
5. `pan-bigomega<n>.txt`: list of all PANs with *n* prime factors, counted with multiplicities. Each line of the file is a pair `(ps, psidx)` where *ps* is a sequence of primes whose product is the number, and `psidx` is the index sequence in *linear form* (i.e., a `None` means that the same prime is repeated once again).
6. `pan-nsfo-bigomega<n>.txt`: like 5, but the PANs listed are those with a non square-free odd part.
7. `sfpan-bigomega<n>.txt`: like 5,  but the PANs listed are square-free.

Moreover, the folllowing are the result of in-progress experiments:

1. `count-opan-omega<n>.txt`: number of odd PANs with *n* distinct prime factors.
