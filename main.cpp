#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <optional>
#include <chrono>
#include <bitset>

using namespace std;

const int n = 128;
const int k = 80;

struct Data {
    map<string, vector<int>> party1_X;
    map<string, vector<int>> party2_X;
    map<string, bitset<k>> party1_s;
    map<string, bitset<k>> party2_s;
    vector<bitset<k>> w_0s;
    vector<bitset<k>> w_1s;
    vector<vector<bitset<k>>> D_s;
    int base_seed;
    map<string, bitset<k>> z1;
    map<string, bitset<k>> z2;
};

vector<int> sample_fv(int v, mt19937& gen) {
    vector<int> res;
    uniform_int_distribution<> dis(0, 1);
    for (int i = 0; i < v; i++) {
        res.push_back(dis(gen));
    }
    return res;
}

bitset<k> sample_f2k(mt19937& gen) {
    uniform_int_distribution<> dis(0, 1);
    bitset<k> val;
    for (int i = 0; i < k; i++) {
        val[i] = dis(gen);
    }
    return val;
}

pair<vector<int>, bitset<k>> prg(int base_seed, bitset<k> inp, int v, vector<int> monic) {
    // TODO: Use all bits from inp
    //mt19937 gen(base_seed + (int)(inp.to_ulong()));
    //mt19937 gen(inp);
    string seed = inp.to_string() + to_string(base_seed);
    seed_seq true_seed(seed.begin(), seed.end());
    mt19937 gen(true_seed);
    vector<int> x = sample_fv(v, gen);
    bitset<k> s = sample_f2k(gen);
    return make_pair(x, s);
}

bool alive(string rr, vector<bitset<n>> tf_in, int t, int diff) {
    for (int i = 0; i < t; i++) {
        bitset<n> val = tf_in[i];
        bool works = true;
        for (int j = 0; j < rr.size(); j++) {
            if (rr[j] == '0' && val[n - 1 - j] == 1 || rr[j] == '1' && val[n - 1 - j] == 0) {
                works = false;
                break;
            }
        }
        if (works) {
            return true;
        }
    }
    return false;
}

vector<int> combine_parts(vector<int> a, vector<int> b) {
    vector<int> res;
    for (int i = 0; i < a.size(); i++) {
        res.push_back(a[i] ^ b[i]);
    }
    return res;
}

vector<int> to_pol(bitset<k> a) {
    vector<int> res;
    bool active = false;
    for (int i = k; i >= 0; i--) {
        if (active) {
            res.push_back(a[i]);
        } else {
            if (a[i]) {
                res.push_back(a[i]);
                active = true;
            }
        }
    }
    //reverse(res.begin(), res.end());
    return res;
}

bitset<k> my_prod(bitset<k> a, bitset<k> b, vector<int> monic) {
    if (a.none() && b.none()) {
        return 0;
    }
    //cout << "My prod called\n";
    //cout << "a: " << a << endl;
    //cout << "b: " << b << endl;
    vector<int> a_pol = to_pol(a);
    vector<int> b_pol = to_pol(b);
    int la = a_pol.size();
    int lb = b_pol.size();
    vector<int> new_pol(la + lb - 1, 0);
    for (int i = 0; i < la; i++) {
        if (!a_pol[i]) {
            continue;
        }
        for (int j = 0; j < lb; j++) {
            if (!b_pol[j]) {
                continue;
            }
            int pa = la - i - 1;
            int pb = lb - j - 1;
            int p = pa + pb;
            int ti = la + lb - p - 2;
            new_pol[ti] = (new_pol[ti] + 1) % 2;
        }
    }
    while (new_pol.size() >= monic.size()) {
        int diff = new_pol.size() - monic.size();
        vector<int> multed;
        for (int el : monic) {
            multed.push_back(el);
        }
        for (int i = 0; i < diff; i++) {
            multed.push_back(0);
        }
        for (int i = 0; i < multed.size(); i++) {
            new_pol[i] = (new_pol[i] + multed[i]) % 2;
        }
        while (new_pol[0] == 0) {
            new_pol.erase(new_pol.begin());
            if (new_pol.size() < monic.size()) {
                break;
            }
        }
    }
    bitset<k> prod;
    reverse(new_pol.begin(), new_pol.end());
    for (int i = 0; i < new_pol.size(); i++) {
        prod[i] = new_pol[i];
    }
    return prod;
}

int find_rank(vector<vector<int>>& A, vector<bitset<k>>& B, vector<int> monic) {
    //cout << "Starting find rank\n";
    int x = 0;
    int y = 0;
    while (true) {
        if (y >= A.size() || x >= A[0].size()) {
            break;
        }
        int base_val = A[y][x];
        if (base_val == 0) {
            bool found = false;
            for (int ty = y + 1; ty < A.size(); ty++) {
                if (A[ty][x] != 0) {
                    vector<int> tmpvec = A[y];
                    A[y] = A[ty];
                    A[ty] = tmpvec;
                    bitset<k> tmp = B[y];
                    B[y] = B[ty];
                    B[ty] = tmp;
                    found = true;
                    break;
                }
            }
            if (!found) {
                x += 1;
            }
            continue;
        }
        for (int ty = y + 1; ty < A.size(); ty++) {
            if (A[ty][x] != 0) {
                B[ty] = B[ty] ^ B[y];
                for (int tx = x; tx < A[0].size(); tx++) {
                    A[ty][tx] = A[ty][tx] ^ A[y][tx];
                }
            }
        }
        x += 1;
        y += 1;
    }
    int c = 0;
    for (vector<int> row : A) {
        bool all_zero = true;
        for (int el : row) {
            if (el != 0) {
                all_zero = false;
                break;
            }
        }
        if (!all_zero) {
            c += 1;
        }
    }
    //cout << "Finished find rank\n";
    return c;
}

vector<int> free_remaining(vector<optional<bitset<k>>> vals, vector<int> row, vector<int> monic) {
    vector<int> free;
    for (int i = 0; i < vals.size(); i++) {
        if (!vals[i].has_value()) {
            if (row[i] == 1) {
                free.push_back(i);
            }
        }
    }
    return free;
}

vector<bitset<k>> sample_solution(vector<vector<int>> A, vector<bitset<k>> B, vector<int> monic, mt19937 gen) {
    //cout << "Starting sampling\n";
    vector<optional<bitset<k>>> vals(A[0].size(), nullopt);
    int y = A.size() - 1;
    while (y >= 0) {
        vector<int> free = free_remaining(vals, A[y], monic);
        if (free.size() > 1) {
            vals[free.back()] = sample_f2k(gen);
            continue;
        }
        bitset<k> e = B[y];
        for (int i = 0; i < vals.size(); i++) {
            optional<bitset<k>> val = vals[i];
            int el = A[y][i];
            if (val.has_value() && el == 1) {
                e = e ^ val.value();
            }
        }
        vals[free[0]] = e;
        y -= 1;
    }
    for (int i = 0; i < vals.size(); i++) {
        if (!vals[i].has_value()) {
            vals[i] = sample_f2k(gen);
        }
    }
    vector<bitset<k>> retvals;
    for (int i = 0; i < vals.size(); i++) {
        retvals.push_back(vals[i].value());
    }
    //cout << "Ended sampling\n";
    return retvals;
}

optional<vector<bitset<k>>> solve_equations(vector<pair<vector<int>, bitset<k>>> linear_equations, vector<int> monic, mt19937 gen) {
    //cout << "Starting solving of equations\n";

    vector<vector<int>> A;
    vector<bitset<k>> B;
    for (pair<vector<int>, bitset<k>> el : linear_equations) {
        A.push_back(el.first);
        B.push_back(el.second);
    }

    //for (vector<int> row : A) {
    //    for (int el : row) {
    //        cout << el << " ";
    //    }
    //    cout << endl;
    //}

    int rank = find_rank(A, B, monic);

    //cout << endl;
    //for (vector<int> row : A) {
    //    for (int el : row) {
    //        cout << el << " ";
    //    }
    //    cout << endl;
    //}
    if (rank < linear_equations.size()) {
        return nullopt;
    }
    vector<bitset<k>> sol = sample_solution(A, B, monic, gen);
    return sol;
}

optional<Data> gen_stuff(vector<bitset<n>> tf_in, vector<int> tf_out, int v, int n, int t, vector<int> monic) {
    Data ret;
    random_device rd;
    mt19937 gen(rd());
    int base_seed = rd() + time(NULL);
    ret.base_seed = base_seed;

    // Step 1
    vector<set<string>> R_s;
    set<string> R0;
    R0.insert("");
    R_s.push_back(R0);

    vector<int> X_e_1 = sample_fv(v, gen);
    bitset<k> s_e_1 = sample_f2k(gen);

    vector<int> X_e_2 = sample_fv(v, gen);
    bitset<k> s_e_2 = sample_f2k(gen);

    map<string, vector<int>> party1_X;
    party1_X[""] = X_e_1;
    map<string, vector<int>> party2_X;
    party2_X[""] = X_e_2;
    map<string, bitset<k>> party1_s;
    party1_s[""] = s_e_1;
    map<string, bitset<k>> party2_s;
    party2_s[""] = s_e_2;


    // Step 2
    vector<vector<bitset<k>>> D_s;
    vector<bitset<k>> w_0s;
    vector<bitset<k>> w_1s;
    for (int i = 1; i < n + 1; i++) {
        // Substep a
        bitset<k> w_i_0 = sample_f2k(gen);
        bitset<k> w_i_1 = w_i_0;
        while (w_i_0 == w_i_1) {
            w_i_1 = sample_f2k(gen);
        }
        w_0s.push_back(w_i_0);
        w_1s.push_back(w_i_1);

        // Substep b
        set<string> R_i;       // Alive
        set<string> R_i_prim;  // Dead
        set<string> R_i_roof;  // Both children alive
        for (string r : R_s.back()) {
            
            string r0 = r + "0";
            string r1 = r + "1";
            //cout << "String: '" << r << "' to '" << r0 << "' and '" << r1 << "'" << endl;
            int numalive = 0;
            for (string rr : {r0, r1}) {
                if (alive(rr, tf_in, t, n - i)) {
                    cout << rr << " is alive" << endl;
                    numalive++;
                    R_i.insert(rr);
                } else {
                    cout << rr << " is dead" << endl;
                    R_i_prim.insert(rr);
                }
            }
            if (numalive == 2) {
                R_i_roof.insert(r);
            }
            cout << endl;
        }

        vector<pair<vector<int>, bitset<k>>> linear_equations;

        // Substep c
        for (string dead_r : R_i_prim) {
            bitset<k> wib = dead_r.back() == '0' ? w_i_0 : w_i_1;
            string par = dead_r.substr(0, dead_r.size() - 1);
            vector<int> X_r = combine_parts(party1_X[par], party2_X[par]);
            bitset<k> s_r = party1_s[par] ^ party2_s[par];
            bitset<k> prod = my_prod(s_r, wib, monic);
            linear_equations.push_back(make_pair(X_r, prod));
        }

        // Substep d
        for (string alive_r : R_i_roof) {
            bitset<k> wr2 = w_i_0;
            while (wr2 == w_i_0 || wr2 == w_i_1) {
                wr2 = sample_f2k(gen);
            }
            vector<int> X_r = combine_parts(party1_X[alive_r], party2_X[alive_r]);
            bitset<k> s_r = party1_s[alive_r] ^ party2_s[alive_r];
            bitset<k> prod = my_prod(s_r, wr2, monic);
            linear_equations.push_back(make_pair(X_r, prod));
        }

        // Substep e-g
        // TODO: Add failure with calculated probability
        optional<vector<bitset<k>>> d_i = solve_equations(linear_equations, monic, gen);
        if (!d_i.has_value()) {
            cout << "Failed to solve equations\n";
            return nullopt;
        }
        D_s.push_back(d_i.value());


        // Substep h: Final
        map<string, bitset<k>> z1;
        for (string r : R_i) {
            string par = r.substr(0, r.size() - 1);
            vector<int> x = party1_X[par];
            bitset<k> total;
            for (int j = 0; j < x.size(); j++) {
                total = total ^ my_prod(x[j], d_i.value()[j], monic);
            }
            total = total ^ my_prod(party1_s[par], r.back() == '0' ? w_i_0 : w_i_1, monic);
            pair<vector<int>, bitset<k>> val = prg(base_seed, total, v, monic);
            if (i == n) {
                z1[r] = total;
            }
            party1_X[r] = val.first;
            party1_s[r] = val.second;
        }


        map<string, bitset<k>> z2;
        for (string r : R_i) {
            string par = r.substr(0, r.size() - 1);
            vector<int> x = party2_X[par];
            bitset<k> total;
            for (int j = 0; j < x.size(); j++) {
                total = total ^ my_prod(x[j], d_i.value()[j], monic);
            }
            total = total ^ my_prod(party2_s[par], r.back() == '0' ? w_i_0 : w_i_1, monic);
            pair<vector<int>, bitset<k>> val = prg(base_seed, total, v, monic);
            if (i == n) {
                z2[r] = total;
            }
            party2_X[r] = val.first;
            party2_s[r] = val.second;
        }
        ret.z1 = z1;
        ret.z2 = z2;
        R_s.push_back(R_i);
    }

    ret.party1_X = party1_X;
    ret.party2_X = party2_X;
    ret.party1_s = party1_s;
    ret.party2_s = party2_s;
    ret.w_0s = w_0s;
    ret.w_1s = w_1s;
    ret.D_s = D_s;


    for (int i = 0; i < t; i++) {
        cout << tf_in[i] << " -> " << tf_out[i] << endl;
    }
    return ret;
}

string to_bitstring(int val, int n) {
    char chars[n];
    for (int i = 0; i < n; i++) {
        chars[n - i - 1] = (val % 2) ? '1' : '0';
        val = val / 2;
    }
    string str(chars, n);
    return str;
}

optional<vector<bitset<k>>> de_rand(map<string, bitset<k>> Z1, map<string, bitset<k>> Z2, vector<bitset<n>> tf_in, vector<int> tf_out, int n, int t, vector<int> monic, mt19937 gen) {
    vector<bitset<k>> z1;
    vector<bitset<k>> z2;
    vector<int> B;
    for (int i = 0; i < t; i++) {
        string bitstring = tf_in[i].to_string();
        z1.push_back(Z1[bitstring]);
        z2.push_back(Z2[bitstring]);
        B.push_back(tf_out[i]);
        //cout << bitstring << endl;
    }
    vector<vector<int>> Z;
    for (int i = 0; i < z1.size(); i++) {
        bitset<k> a = z1[i];
        bitset<k> b = z2[i];
        vector<int> res;
        for (int j = 0; j < k; j++) {
            res.push_back((a[j] % 2) ^ (b[j] % 2));
            //a = a / 2;
            //b = b / 2;
        }
        reverse(res.begin(), res.end());
        Z.push_back(res);
    }
    vector<pair<vector<int>, bitset<k>>> linear_system;
    for (int i = 0; i < Z.size(); i++) {
        linear_system.push_back(make_pair(Z[i], B[i]));
    }
    return solve_equations(linear_system, monic, gen);
}

int main() {
    int t = 20;
    int v = 128;

    vector<int> monic = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1};
    //vector<int> monic = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1};
    //vector<int> monic = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1};
    //vector<int> monic = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1};
    //vector<int> monic = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1};
    //vector<int> monic = {1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1};
    
    vector<bitset<n>> tf_in;
    vector<int> tf_out;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, (1 << 30) - 1);
    uniform_int_distribution<> binary_dis(0, 1);

    // Generate multi-point function
    for (int i = 0; i < t; i++) {
        bitset<n> gened;
        for (int j = 0; j < n; j++) {
            gened[j] = binary_dis(gen);
        }
        tf_in.push_back(gened);
        tf_out.push_back(dis(gen));
    }

    chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
    Data data;
    vector<bitset<k>> g;
    while (true) {
        while (true) {
            optional<Data> val = gen_stuff(tf_in, tf_out, v, n, t, monic);
            if (val.has_value()) {
                data = val.value();
                break;
            }
            cout << "Retrying\n";
        }
        mt19937 gen2(rd());
        optional<vector<bitset<k>>> maybe_g = de_rand(data.z1, data.z2, tf_in, tf_out, n, t, monic, gen2);
        if (maybe_g.has_value()) {
            g = maybe_g.value();
            cout << "Derandomized\n";
            break;
        }
    }
    chrono::high_resolution_clock::time_point end = chrono::high_resolution_clock::now();

    chrono::duration<double, milli> duration = chrono::duration<double, milli>(end - start);
    cout << "Execution time: " << duration.count() << " ms" << endl;
    
    return 0;
}
