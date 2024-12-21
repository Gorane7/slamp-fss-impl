import random
import time

import galois

import hashlib

import pickle

VERBOSE = 1


def sample_fv(v, GF):
    res = []
    for i in range(v):
        res.append(GF(random.randint(0, 1)))
    return res

def sample_fv2k(v, k, GF):
    res = []
    for i in range(v):
        res.append(sample_f2k(k, GF))
    return res

def sample_f2k(k, GF):
    return GF(random.randint(1, 2**k - 1))

def alive(start, fun):
    for a, b in fun:
        if a.startswith(start):
            return True
    return False


def prg(base_seed, inp, v, k, GF):
    random.seed(base_seed + int(inp))
    x = sample_fv(v, GF)
    s = sample_f2k(k, GF)
    #s = GF(random.randint(0, 1))  # This didn't find a solution, i think this is not correct
    return x, s

def disp_matrix(A, msg):
    print(msg)
    for row in A:
        print(f"{[int(x) for x in row]}")
    print()


def find_rank(A, B, GF):
    #disp_matrix(A, "Starting to find rank of matrix")
    y = 0
    x = 0
    # I can convert matrices to row echelon form now!?
    while True:
        if y >= len(A) or x >= len(A[0]):
            break
        base_val = A[y][x]
        #print(f"Coordinates are row {y}, column {x}")
        if base_val == GF(0):
            # Find row to switch
            found = False
            for ty in range(y + 1, len(A)):
                if A[ty][x] != GF(0):
                    A[ty], A[y] = A[y], A[ty]
                    B[ty], B[y] = B[y], B[ty]
                    #disp_matrix(A, "After row swap")
                    found = True
                    break
            if not found:
                x += 1
            continue
        # Go through all bottom rows and subtract
        for ty in range(y + 1, len(A)):
            if A[ty][x] != GF(0):
                B[ty] -= B[y]
                for tx in range(x, len(A[0])):
                    A[ty][tx] -= A[y][tx]
        #disp_matrix(A, "After making column down zero")
        x += 1
        y += 1
    c = 0
    for row in A:
        all_zero = True
        for el in row:
            if el != GF(0):
                all_zero = False
        if not all_zero:
            c += 1
    #print(f"Rank of matrix is: {c}")
    return c


def free_remaining(vals, row, GF):
    free = []
    for i, (val, el) in enumerate(zip(vals, row)):
        if val is None:
            if el == GF(1):
                free.append(i)
    return free


def sample_solution(A, B, k, GF):
    if VERBOSE >= 1:
        print(f"B: {[int(x) for x in B]}")
        disp_matrix(A, "Sampling from")
    vals = [None for i in range(len(A[0]))]
    y = len(A) - 1
    while y >= 0:
        if VERBOSE >= 2:
            print(vals)
        free = free_remaining(vals, A[y], GF)
        if len(free) > 1:
            # Take random
            vals[free[-1]] = sample_f2k(k, GF)
            continue
        # Calculate
        e = B[y]
        for val, el in zip(vals, A[y]):
            if val is not None and el == GF(1):
                e = e - val
        vals[free[0]] = e
        y -= 1
    vals = [(val if val is not None else sample_f2k(k, GF)) for val in vals]
    return vals


def solve_equations(linear_equations, k, GF):
    v = len(linear_equations[0][0])
    if VERBOSE >= 2:
        print(f"Need to solve, length of d (same as v) is {v}")
        print(f"Have {len(linear_equations)} rows in linear equation")
        print(f"Each row has length {len(linear_equations[0][0])}, plus singular element ({linear_equations[0][1]})")
    B = []
    for a, b in linear_equations:
        B.append(b)
    if VERBOSE >= 2:
        print(f"B: {[int(b) for b in B]}")
    A = []
    for a, b in linear_equations:
        A.append(a)
    if VERBOSE >= 2:
        for row in A:
            print(f"Row: {[int(el) for el in row]}")

    rank = find_rank(A, B, GF)
    if VERBOSE >= 2:
        print(f"Rank is: {rank}")
        disp_matrix(A, "After gaussian eliminiation")
    if rank < len(linear_equations):
        print(f"Failed to solve")
        #exit()
        return False
    sol = sample_solution(A, B, k, GF)
    #print(f"Sampled solution: {[x for x in sol]}")
    return sol

def combine_parts(xa, xb):
    return [a + b for a, b in zip(xa, xb)]


def gen_stuff(tf, v, k, n, GF):
    base_seed = int(time.time()) + random.randint(0, 100000)
    #tf = {x[0]: x[1] for x in tf}
    # F^v_2k: Vector of length v, where elements are from a field of length 2^k (so basically bit strings of length k)
    # t: number of points in multi-point function
    
    # Step 1
    R_s = []
    R0 = set()
    R0.add('')
    R_s.append(R0)
    
    X_e_1 = sample_fv(v, GF)
    s_e_1 = sample_f2k(k, GF)

    X_e_2 = sample_fv(v, GF)
    s_e_2 = sample_f2k(k, GF)
    
    party1_X = {
        "": X_e_1
    }
    party2_X = {
        "": X_e_2
    }
    party1_s = {
        "": s_e_1
    }
    party2_s = {
        "": s_e_2
    }
    if VERBOSE >= 2:
        print([int(a) for a in X_e_1], s_e_1)
        print([int(a) for a in X_e_2], s_e_2)
        print()
    
    # Step 2
    D_s = []
    w_0s = []
    w_1s = []
    for i in range(1, n + 1):
        
        # Substep a
        w_i_0 = sample_f2k(k, GF)
        w_i_1 = w_i_0
        while w_i_0 == w_i_1:
            w_i_1 = sample_f2k(k, GF)
        if VERBOSE >= 2:
            print("wi0_1:", int(w_i_0), ",", int(w_i_1))
        w_0s.append(w_i_0)
        w_1s.append(w_i_1)
        
        # Substep b
        R_i = set()  # Alive
        R_i_prim = set()  # Dead
        R_i_roof = set()  # Both children alive
        for r in R_s[-1]:
            r0 = r + "0"
            r1 = r + "1"
            numalive = 0
            for rr in [r0, r1]:
                if alive(rr, tf):
                    numalive += 1
                    R_i.add(rr)
                else:
                    R_i_prim.add(rr)
            if numalive == 2:
                R_i_roof.add(r)
        if VERBOSE >= 2:
            print(f"Alive: {R_i}")
            print(f"Dead: {R_i_prim}")
            print(f"Both children alive: {R_i_roof}")
        
        linear_equations = []

        # Substep c
        for dead_r in R_i_prim:
            wib = w_i_0 if dead_r[-1] == "0" else w_i_1
            X_r = combine_parts(party1_X[dead_r[:-1]], party2_X[dead_r[:-1]])
            s_r = party1_s[dead_r[:-1]] + party2_s[dead_r[:-1]]
            linear_equations.append((X_r, s_r * wib))
        
        # Substep d
        for alive_r in R_i_roof:
            wr2 = w_i_0
            while wr2 == w_i_0 or wr2 == w_i_1:
                wr2 = sample_f2k(k, GF)
            X_r = combine_parts(party1_X[alive_r], party2_X[alive_r])
            s_r = party1_s[alive_r] + party2_s[alive_r]
            linear_equations.append((X_r, s_r * wr2))
        
        # Substep e-g (fail if not full rank, also fail with some probability)
        # TODO: Add failure with calculated probability
        d_i = solve_equations(linear_equations, k, GF)
        if not d_i:
            return 0
        if VERBOSE >= 2:
            print()
            print("Equations")
            [print([int(a) for a in x[0]], ' = ', x[1]) for x in linear_equations]
            print("Answers")
            print([int(a) for a in d_i])
            print()
        D_s.append(d_i)
        
        # Substep h: Final
        next_values_1 = []
        z1 = {}
        for r in R_i:
            x = party1_X[r[:-1]]
            total = GF(0)
            assert len(x) == len(d_i)

            # Apply solutions of linear system
            for a, b in zip(x, d_i):
                total += a * b
            total += party1_s[r[:-1]] * (w_i_0 if r[-1] == "0" else w_i_1)
            val = prg(base_seed, total, v, k, GF)
            if i == n:
                z1[r] = total
            next_values_1.append((r, val[0], val[1]))

        next_values_2 = []
        z2 = {}
        for r in R_i:
            x = party2_X[r[:-1]]
            total = GF(0)
            assert len(x) == len(d_i)
            for a, b in zip(x, d_i):
                total += a * b
            total += party2_s[r[:-1]] * (w_i_0 if r[-1] == "0" else w_i_1)
            val = prg(base_seed, total, v, k, GF)
            if i == n:
                z2[r] = total
            next_values_2.append((r, val[0], val[1]))
        
        if VERBOSE >= 2:
            print(next_values_1)
            print(next_values_2)
        
        # Extra substep, save data
        # Currently this is a bit ugly, for the last step we should have a new type of thing called z, we just save it as X and s remains None
        for r, X, s in next_values_1:
            party1_X[r] = X
            party1_s[r] = s
        for r, X, s in next_values_2:
            party2_X[r] = X
            party2_s[r] = s
        R_s.append(R_i)
        if VERBOSE >= 2:
            print()
    return party1_X, party1_s, party2_X, party2_s, w_0s, w_1s, D_s, base_seed, z1, z2


def evaluate(sk, inp, n, v, k, GF, base_seed):
    names = ["Xe", "se", "w0", "w1", "d", "g"]
    X0 = sk[0]
    s0 = sk[1]
    X_s = [X0]
    s_s = [s0]
    if VERBOSE >= 3:
        print(f"Starting evaluate for {inp}")
        print(f"X0: {[int(x) for x in X0]}, s0: {s0}")
    for i in range(1, n + 1):
        if i >= len(inp) + 1:
            break
        z_i = GF(0)
        for a, b in zip(X_s[-1], sk[4][i - 1]):
            z_i += a * b
        if VERBOSE >= 3:
            print(f"pre z: {z_i}")
        z_i += s_s[-1] * (sk[2][i - 1] if inp[i - 1] == "0" else sk[3][i - 1])
        if VERBOSE >= 3:
            print(f"z: {z_i}")

        #if i < n:
        X_i, s_i = prg(base_seed, z_i, v, k, GF)
        if VERBOSE >= 3:
            print(f"X after PRG: {[int(x) for x in X_i]}")
            print(f"s after PRG: {s_i}")
            print()
        X_s.append(X_i)
        s_s.append(s_i)
    z = GF(0)
    z_i = to_bits(z_i, k)
    for a, b in zip(z_i, sk[5]):
        z += GF(int(a)) * b
    if VERBOSE >= 3:
        print(f"Final z: {z}")
        print(f"Ended evaluate for {inp}")
        print()
    return int(z)


def tree_eval(inp_s, n, sk1, sk2, v, k, GF, base_seed):
    if len(inp_s) > n:
        return
    p1 = evaluate(sk1, inp_s, n, v, k, GF, base_seed)
    p2 = evaluate(sk2, inp_s, n, v, k, GF, base_seed)
    print(f"p1({inp_s}) = {p1}")
    print(f"p2({inp_s}) = {p2}")
    print(f"Combined: {p1 ^ p2}")
    print()
    tree_eval(inp_s + "0", n, sk1, sk2, v, k, GF, base_seed)
    tree_eval(inp_s + "1", n, sk1, sk2, v, k, GF, base_seed)


def define_plus(GF, am):
    all_ans = ""
    for a in range(am):
        for b in range(am):
            c = GF(a) + GF(b)
            #print(f"{a} + {b} = {int(c)}")
            all_ans += str(int(c))
    print(all_ans)
    print(hashlib.md5(all_ans.encode('utf-8')).hexdigest())
    print()

def to_bits(z, k):
    return bin(int(z))[2:].zfill(k)


def de_rand(Z1, Z2, tf, k, GF):
    z1 = []
    z2 = []
    B = []
    for source, sink in tf:
        z1.append(to_bits(Z1[source], k))
        z2.append(to_bits(Z2[source], k))
        B.append(sink)
    Z = []
    for a, b in zip(z1, z2):
        Z.append([int(q) for q in bin(int(a, 2) ^ int(b, 2))[2:].zfill(k)])
    linear_system = []
    for row, b in zip(Z, B):
        linear_system.append(([GF(x) for x in row], GF(b)))
    return solve_equations(linear_system, k, GF)


if __name__ == '__main__':
    times = []
    gen_stuff_time = 0
    gen_stuff_times = 0

    k = 6
    n = 4  # Number of bits in input
    t = 3  # Number of points
    v = 8
    tf = []
    used = set()
    for i in range(t):
        ss = "".join([random.choice("01") for q in range(n)])
        while ss in used:
            ss = "".join([random.choice("01") for q in range(n)])
        val = random.randint(1, 2**k - 1)
        tf.append((ss, val))
        used.add(ss)

    
    GF = galois.GF(2**k)  # Initializing this takes like a second
    start = time.time()

    if VERBOSE >= 1:
        print("Target function")
        [print(f"f({a}) = {b}") for a, b in tf]
        print()

    while True:
        while True:
            gen_stuff_time -= time.time()

            data = gen_stuff(tf, v, k, n, GF)

            gen_stuff_times += 1
            gen_stuff_time += time.time()
            if data:
                break
            print("FAILED TO GEN STUFF")
        # party1_X, party1_s, party2_X, party2_s, w_0s, w_1s, D_s, base_seed, z1, z2

        base_seed = data[7]
        if g := de_rand(data[8], data[9], tf, k, GF):
            break

        print("FAILED TO SOLVE WITH GENED STUFF")
        print(tf)
        if VERBOSE >= 2:
            print("Unable to find solution, restarting process")
    
    sk1 = [data[0][""], data[1][""], data[4], data[5], data[6], g]
    sk2 = [data[2][""], data[3][""], data[4], data[5], data[6], g]
    
    end = time.time()
    
    names = ["Xe", "se", "w0", "w1", "d", "g"]
    if VERBOSE >= 1:
        print(f"g: {[int(a) for a in g]}")
        print()
        print("SK1")
        for a, name in zip(sk1, names):
            if isinstance(a, list):
                if isinstance(a[0], list):
                    a = [[int(c) for c in b] for b in a]
                else:
                    a = [int(b) for b in a]
            print(f"{name}: {a}")
        print()
        print("SK2")
        for a, name in zip(sk2, names):
            if isinstance(a, list):
                if isinstance(a[0], list):
                    a = [[int(c) for c in b] for b in a]
                else:
                    a = [int(b) for b in a]
            print(f"{name}: {a}")
    if VERBOSE >= 1:
        print(f"Time taken: {round(end - start, 2)} s")
    times.append(end - start)
    avg = sum(times) / len(times)
    print(f"On average took: {round(avg, 2)} s")
    print(f"Gen stuff time: {round(gen_stuff_time, 2)} s")
    print(f"Did gen stuff {gen_stuff_times} times")
    print(f"All times: {times}")
    
    
    
