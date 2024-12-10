#include <RcppArmadillo.h>
#include <unordered_set>
#define ARMA_64BIT_WORD

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Union-Find structure
struct UnionFind {
    std::vector<int> parent;
    std::vector<int> size;
    UnionFind(int n) : parent(n), size(n,1) {
        for (int i = 0; i < n; i++) parent[i] = i;
    }
    int find_set(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    void union_set(int a, int b) {
        a = find_set(a);
        b = find_set(b);
        if (a != b) {
            if (size[a] < size[b]) std::swap(a,b);
            parent[b] = a;
            size[a] += size[b];
        }
    }
};

// Fenwicks tree for prefix sums
struct Fenw {
    int n;
    std::vector<double> fenw;
    Fenw(int n) : n(n), fenw(n+1,0.0) {}
    void update(int i, double delta) {
        for (; i<=n; i+=(i & (-i))) fenw[i]+=delta;
    }
    double prefix_sum(int i) {
        double s=0.0;
        for (; i>0; i-=(i & (-i))) s+=fenw[i];
        return s;
    }
    double total_sum() {
        return prefix_sum(n);
    }
    int find(double value) {
        int pos=0;
        int bit_mask=1 << (int)std::floor(std::log2(n));
        for (; bit_mask>0; bit_mask>>=1) {
            int next_pos=pos+bit_mask;
            if (next_pos<=n && fenw[next_pos]<value) {
                value-=fenw[next_pos];
                pos=next_pos;
            }
        }
        return pos+1;
    }
};

// [[Rcpp::export]]
List rec_coarsen_impl(const S4& W, int T, 
                      Nullable<NumericVector> phi_ = R_NilValue,
                      bool deterministic_first_edge=false,
                      bool verbose=false, 
                      int seed = NA_INTEGER) {
    IntegerVector dims = W.slot("Dim");
    int N = dims[0];

    IntegerVector Wi = W.slot("i");
    IntegerVector Wp = W.slot("p");
    NumericVector Wx = W.slot("x");

    // Convert to triplet edges i<j
    std::vector<int> Wj(Wi.size());
    for (int col=0; col<N; col++) {
        int start_idx = Wp[col];
        int end_idx = Wp[col+1];
        for (int idx=start_idx; idx<end_idx; idx++) {
            Wj[idx] = col;
        }
    }

    std::vector<int> Ei;
    std::vector<int> Ej;
    std::vector<double> Ew;
    for (int k=0; k<(int)Wi.size(); k++) {
        int r = Wi[k]+1;
        int c = Wj[k]+1;
        if (r<c) {
            Ei.push_back(r);
            Ej.push_back(c);
            Ew.push_back(Wx[k]);
        }
    }

    int M=(int)Ei.size();
    std::vector<double> phi(M);
    if (phi_.isNull()) {
        for (int i=0; i<M; i++) phi[i]=Ew[i];
    } else {
        NumericVector phi_in(phi_);
        if ((int)phi_in.size()!=M) stop("phi length must match number of edges");
        for (int i=0; i<M; i++) phi[i]=phi_in[i];
    }

    if (!R_IsNA(seed)) {
        Rcpp::RNGScope scope;
        srand((unsigned)seed);
    }

    std::vector<std::vector<int>> vertex_edge_map(N+1);
    for (int e=0; e<M; e++) {
        vertex_edge_map[Ei[e]].push_back(e);
        vertex_edge_map[Ej[e]].push_back(e);
    }

    UnionFind uf(N+1);
    std::vector<bool> cand(M,true);

    Fenw fenw(M);
    for (int e=0; e<M; e++) {
        fenw.update(e+1, phi[e]);
    }
    double total_phi = fenw.total_sum();

    for (int iter=1; iter<=T; iter++) {
        if (total_phi<=0) break;

        int chosen_edge_id=-1;
        if (deterministic_first_edge && iter==1) {
            for (int e=0; e<M; e++) {
                if (cand[e] && phi[e]>0) {
                    chosen_edge_id=e;
                    break;
                }
            }
            if (chosen_edge_id<0) {
                continue;
            }
        } else {
            double r = R::runif(0,1)*total_phi;
            int fenw_idx = fenw.find(r);
            chosen_edge_id = fenw_idx-1;
            if (!cand[chosen_edge_id] || phi[chosen_edge_id]<=0) {
                // Should not happen if fenw is maintained correctly,
                // but in worst case just continue
                continue;
            }
        }

        int ei=Ei[chosen_edge_id];
        int ej=Ej[chosen_edge_id];

        int pi=uf.find_set(ei);
        int pj=uf.find_set(ej);

        if (pi!=pj) {
            uf.union_set(pi,pj);
            // Remove neighborhood
            std::unordered_set<int> neigh;
            for (auto x: vertex_edge_map[ei]) if (cand[x]) neigh.insert(x);
            for (auto x: vertex_edge_map[ej]) if (cand[x]) neigh.insert(x);

            double removed_pot=0.0;
            for (auto e: neigh) {
                if (cand[e]) {
                    cand[e]=false;
                    fenw.update(e+1, -phi[e]);
                    removed_pot += phi[e];
                    phi[e]=0.0;
                }
            }
            total_phi -= removed_pot;
            if (total_phi<1e-14) total_phi=0.0;
        } else {
            // same set
            if (cand[chosen_edge_id]) {
                cand[chosen_edge_id]=false;
                fenw.update(chosen_edge_id+1, -phi[chosen_edge_id]);
                total_phi -= phi[chosen_edge_id];
                phi[chosen_edge_id]=0.0;
            }
        }
    }

    // final compression
    for (int v=1; v<=N; v++) uf.find_set(v);

    std::unordered_map<int,int> new_id_map;
    int idx=0;
    for (int v=1; v<=N; v++) {
        int p=uf.find_set(v);
        if (new_id_map.find(p)==new_id_map.end()) {
            idx++;
            new_id_map[p]=idx;
        }
    }

    IntegerVector mapping(N);
    for (int v=1; v<=N; v++) {
        mapping[v-1]=new_id_map[uf.find_set(v)];
    }

    int n=idx;

    std::vector<std::vector<int>> groups(n+1);
    for (int v=1; v<=N; v++) {
        groups[mapping[v-1]].push_back(v);
    }

    int total_count=0;
    for (int i=1; i<=n; i++) total_count+=(int)groups[i].size();

    IntegerVector row_idx(total_count);
    IntegerVector col_idx(total_count);
    NumericVector x_vals(total_count);

    {
        int pos=0;
        for (int ci=1; ci<=n; ci++) {
            int sz=(int)groups[ci].size();
            double val=1.0/std::sqrt((double)sz);
            for (int k=0; k<sz; k++) {
                row_idx[pos]=ci-1; // 0-based
                col_idx[pos]=groups[ci][k]-1;
                x_vals[pos]=val;
                pos++;
            }
        }
    }

    // Build C as a dgCMatrix
    // We have triplets (row_idx,col_idx,x_vals)
    // Convert to armadillo sp_mat
    {
        arma::umat locs(2,total_count);
        arma::vec vals_arma(total_count);
        for (int i=0; i<total_count; i++) {
            locs(0,i)= (arma::uword)row_idx[i];
            locs(1,i)= (arma::uword)col_idx[i];
            vals_arma[i]=x_vals[i];
        }
        arma::sp_mat C_sp(locs, vals_arma, (arma::uword)n, (arma::uword)N);

        // Convert W to sp_mat
        // We already have W in dgCMatrix form
        // Let's convert W to arma::sp_mat
        // We know W is dgCMatrix: i, p, x
        IntegerVector Wi2 = W.slot("i");
        IntegerVector Wp2 = W.slot("p");
        NumericVector Wx2 = W.slot("x");
        arma::uword W_nnz = (arma::uword)Wx2.size();
        arma::umat W_locs(2,W_nnz);
        arma::vec W_vals(W_nnz);
        {
            for (arma::uword c=0; c<(arma::uword)N; c++) {
                for (int idx=Wp2[c]; idx<Wp2[c+1]; idx++) {
                    W_locs(0,(arma::uword)idx)=(arma::uword)Wi2[idx];
                    W_locs(1,(arma::uword)idx)=(arma::uword)c;
                    W_vals[(arma::uword)idx]=Wx2[idx];
                }
            }
        }
        arma::sp_mat W_sp(W_locs, W_vals, (arma::uword)N,(arma::uword)N);

        // Compute W_c = C W C^T
        arma::sp_mat CW = C_sp*W_sp;
        arma::sp_mat W_c_sp = CW*C_sp.t();

        // Compute L_c = C L C^T
        // L = D - W
        // compute d (row sums of W)
        arma::vec d(N,arma::fill::zeros);
        for (arma::uword c=0; c<(arma::uword)N; c++) {
            for (arma::uword idx=Wp2[c]; idx<(arma::uword)Wp2[c+1]; idx++) {
                d[(arma::uword)Wi2[idx]] += Wx2[idx];
            }
        }
        // D - W
        // L = D - W
        // D is diagonal
        // Just form L_sp = D - W_sp
        arma::sp_mat D_sp = arma::sp_mat(arma::diagmat(d));
        arma::sp_mat L_sp = D_sp - W_sp;
        arma::sp_mat L_c_sp = C_sp*L_sp*C_sp.t();

        // Convert back to dgCMatrix
        auto to_dgC = [&](const arma::sp_mat &X)->S4 {
            S4 out("dgCMatrix");
            int nnz=(int)X.n_nonzero;
            // i, p, x
            IntegerVector iX(nnz);
            IntegerVector pX(X.n_cols+1);
            NumericVector xX(nnz);
            // X.col_ptrs and X.row_indices are zero-based
            for (int c=0; c<(int)X.n_cols; c++) {
                pX[c]=(int)X.col_ptrs[c];
            }
            pX[X.n_cols]=(int)X.col_ptrs[X.n_cols]; // last one
            for (int k=0; k<nnz; k++) {
                iX[k]=(int)X.row_indices[k];
                xX[k]=X.values[k];
            }
            out.slot("i")=iX;
            out.slot("p")=pX;
            out.slot("x")=xX;
            out.slot("Dim")=IntegerVector::create((int)X.n_rows,(int)X.n_cols);
            return out;
        };

        S4 C_out("dgCMatrix");
        {
            // C_sp also needs conversion
            int c_nnz=(int)C_sp.n_nonzero;
            IntegerVector iC(c_nnz);
            IntegerVector pC((int)C_sp.n_cols+1);
            NumericVector xC(c_nnz);
            for (int c=0; c<(int)C_sp.n_cols; c++) {
                pC[c]=(int)C_sp.col_ptrs[c];
            }
            pC[C_sp.n_cols]=(int)C_sp.col_ptrs[C_sp.n_cols];
            for (int k=0; k<c_nnz; k++) {
                iC[k]=(int)C_sp.row_indices[k];
                xC[k]=C_sp.values[k];
            }
            C_out.slot("i")=iC;
            C_out.slot("p")=pC;
            C_out.slot("x")=xC;
            C_out.slot("Dim")=IntegerVector::create(n,N);
        }

        S4 W_c = to_dgC(W_c_sp);
        S4 L_c = to_dgC(L_c_sp);

        return List::create(
            _["C"]=C_out,
            _["W_c"]=W_c,
            _["L_c"]=L_c,
            _["mapping"]=mapping
        );
    }
}