

def main():
    
    inp = "x^80 + x^48 + x^46 + x^42 + x^41 + x^38 + x^33 + x^32 + x^31 + x^29 + x^26 + x^25 + x^24 + x^22 + x^21 + x^20 + x^17 + x^15 + x^13 + x^11 + x^9 + x^8 + x^6 + x^5 + x^4 + x^2 + 1"
    #inp = "x^40 + x^23 + x^21 + x^18 + x^16 + x^15 + x^13 + x^12 + x^8 + x^5 + x^3 + x + 1"
    #inp = "x^32 + x^15 + x^9 + x^7 + x^4 + x^3 + 1"
    #inp  = "x^31 + x^3 + 1"
    #inp = "x^30 + x^17 + x^16 + x^13 + x^11 + x^7 + x^5 + x^3 + x^2 + x + 1"
    #inp = "x^10 + x^6 + x^5 + x^3 + x^2 + x + 1"
    #inp = "x^6 + x^4 + x^3 + x + 1"

    inp = inp.split(" + ")
    size = int(inp[0].split("^")[1]) + 1
    pol = [0 for i in range(size)]
    for el in inp:
        if el == "1":
            pol[-1] = 1
        elif el == "x":
            pol[-2] = 1
        else:
            v = int(el.split("^")[1])
            pol[size - 1 - v] = 1
    print("{" + ", ".join([str(x) for x in pol]) + "}")

if __name__ == '__main__':
    main()
