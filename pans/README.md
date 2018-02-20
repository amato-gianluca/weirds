# Primitive Weirds and Abundant Numbers

Counts and lists of primitive abundant numbers. File with extension `.xz` are compressed in the [Xz format](https://en.wikipedia.org/wiki/Xz).

1. `count-opan-bigomega<n>.txt`: number of odd PANs with *n* prime factos, counted with multiplicities. Each line of the file has the form `[p1, ..., pn] m` meaning that there are `m` odd PANs whose factorization starts with `p1 ... pn`.
2. `count-pan-bigomega<n>.txt`: like 1, but for all PANs, both odd and even.
3. `count-sfopan-omega<n>.txt`: like 1, but for square-free odd PANs.
4. `count-sfopan-omega<n>.txt`: like 1, but for square-free PANs, both odd and even.
5. `pan-bigomega<n>.txt`: list of all PANs with *n* prime factors, counted with multiplicities. Each line of the file is a pair `(ps, psidx)` where *ps* is a sequence of primes whose product is the number, and `psidx` is the index sequence in *linear form* (i.e., a `None` means that the same prime is repeated once again).
