## General
This is an implementation of SLAMP-FSS based on https://eprint.iacr.org/2024/1394.pdf. There are implementations for both Python and C++.
### Getting irreducible polynomials
The C++ version is unable to figure out what polynomial to use, so this needs to be given to it. To get a working polynomial, you can use Python's galois library.
```
import galois
a = galois.GF(2**80)
a.properties
```
This will output something like
```
Galois Field:
  name: GF(2^80)
  characteristic: 2
  degree: 80
  order: 1208925819614629174706176
  irreducible_poly: x^80 + x^48 + x^46 + x^42 + x^41 + x^38 + x^33 + x^32 + x^31 + x^29 + x^26 + x^25 + x^24 + x^22 + x^21 + x^20 + x^17 + x^15 + x^13 + x^11 + x^9 + x^8 + x^6 + x^5 + x^4 + x^2 + 1
  is_primitive_poly: True
  primitive_element: x
```
Take the irreducible_poly value from here and put it into `poly2vec.py`. That will output a correctly formatted C++ vector that represents the same polynomial, this should then be copied to `main.cpp` into the beginning of main.