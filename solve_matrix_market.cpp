#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "umfpack.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// Code from UMFPACK: compute the residual, r = Ax-b and return max_norm(r)
static double resid(int32_t n,
                    const double *x,
                    const double *b,
                    double *r,
                    int32_t Ap[],
                    int32_t Ai[],
                    double Ax[]) {
    int32_t i, j, p;
    double norm;
    for (i = 0; i < n; i++) {
        r[i] = -b[i];
    }
    for (j = 0; j < n; j++) {
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            i = Ai[p];
            r[i] += Ax[p] * x[j];
        }
    }
    norm = 0.;
    for (i = 0; i < n; i++) {
        norm = MAX(ABS(r[i]), norm);
    }
    return norm;
}

struct CooMatrix {
    int32_t m;
    int32_t nnz;
    int32_t nnz_max;
    bool symmetric;
    bool triangular;
    std::vector<int32_t> indices_i;
    std::vector<int32_t> indices_j;
    std::vector<double> values_aij;

    inline static std::unique_ptr<CooMatrix> make_new(int32_t m, int32_t nnz_max) {
        return std::unique_ptr<CooMatrix>{new CooMatrix{
            m,
            0,
            nnz_max,
            false,
            false,
            std::vector<int32_t>(nnz_max, 0),
            std::vector<int32_t>(nnz_max, 0),
            std::vector<double>(nnz_max, 0.0),
        }};
    }

    void put(int32_t i, int32_t j, double aij) {
        if (i < 0 || i >= this->m) {
            throw "CooMatrix::put: index of row is outside range";
        }
        if (j < 0 || j >= this->m) {
            throw "CooMatrix::put: index of column is outside range";
        }
        if (this->nnz >= this->nnz_max) {
            throw "CooMatrix::put: max number of items has been exceeded";
        }
        this->indices_i[this->nnz] = i;
        this->indices_j[this->nnz] = j;
        this->values_aij[this->nnz] = aij;
        this->nnz++;
    }

    void mat_vec_mul(std::vector<double> &v, double alpha, const std::vector<double> &u) {
        if (v.size() != static_cast<size_t>(this->m)) {
            throw "sp_mat_vec_mul: size of v must be equal to the dimension of a";
        }
        if (u.size() != static_cast<size_t>(this->m)) {
            throw "sp_mat_vec_mul: size of u must be equal to the dimension of a";
        }
        // v = alpha * A * u
        std::fill(v.begin(), v.end(), 0.0);
        for (size_t k = 0; k < static_cast<size_t>(this->nnz); k++) {
            auto i = this->indices_i[k];
            auto j = this->indices_j[k];
            auto aij = this->values_aij[k];
            v[i] += alpha * aij * u[j];
            if (this->triangular) {
                if (i != j) {
                    v[j] += alpha * aij * u[i];
                }
            }
        }
    }
};

std::unique_ptr<CooMatrix> read_matrix_market(const std::string &filename) {
    FILE *f = fopen(filename.c_str(), "r");
    if (f == NULL) {
        std::cout << "filename = '" << filename << "'" << std::endl;
        throw "read_matrix_market: cannot open file";
    }

    const int line_max = 500;
    char line[line_max];

    if (fgets(line, line_max, f) == NULL) {
        fclose(f);
        throw "read_matrix_market: cannot read any line in the file";
    }

    char mm[24], opt[24], fmt[24], kind[24], sym[24];
    int nread = sscanf(line, "%24s %24s %24s %24s %24s", mm, opt, fmt, kind, sym);
    if (nread != 5) {
        fclose(f);
        throw "read_matrix_market: number of tokens in the header is incorrect";
    }
    if (strncmp(mm, "%%MatrixMarket", 14) != 0) {
        fclose(f);
        throw "read_matrix_market: header must start with %%MatrixMarket";
    }
    if (strncmp(opt, "matrix", 6) != 0) {
        fclose(f);
        throw "read_matrix_market: option must be \"matrix\"";
    }
    if (strncmp(fmt, "coordinate", 10) != 0) {
        fclose(f);
        throw "read_matrix_market: type must be \"coordinate\"";
    }
    if (strncmp(kind, "real", 4) != 0) {
        fclose(f);
        throw "read_matrix_market: number kind must be \"real\"";
    }
    if (strncmp(sym, "general", 7) != 0 && strncmp(sym, "symmetric", 9) != 0) {
        fclose(f);
        throw "read_matrix_market: matrix must be \"general\" or \"symmetric\"";
    }

    std::unique_ptr<CooMatrix> coo;
    bool symmetric = strncmp(sym, "symmetric", 9) == 0;

    bool initialized = false;
    size_t m, n, nnz_max, i, j;
    double x;

    while (fgets(line, line_max, f) != NULL) {
        if (!initialized) {
            if (line[0] == '%') {
                continue;
            }
            nread = sscanf(line, "%zu %zu %zu", &m, &n, &nnz_max);
            if (nread != 3) {
                fclose(f);
                throw "read_matrix_market: cannot parse the dimensions (m,n,nnz)";
            }
            if (symmetric) {
                // umfpack_di_triplet_to_col requires both sides of the diagonal, and MatrixMarket
                // is lower triangular. So, let's allocate more space than needed.
                nnz_max *= 2;
            }
            coo = CooMatrix::make_new(m, nnz_max);
            coo->symmetric = symmetric;
            coo->triangular = false; // umfpack_di_triplet_to_col requires full matrix
            initialized = true;
        } else {
            nread = sscanf(line, "%zu %zu %lg", &i, &j, &x);
            if (nread != 3) {
                fclose(f);
                throw "read_matrix_market: cannot parse the values (i,j,x)";
            }
            coo->put(i - 1, j - 1, x);
            if (symmetric) {
                // umfpack_di_triplet_to_col requires full matrix
                coo->put(j - 1, i - 1, x);
            }
        }
    }

    fclose(f);
    return coo;
}

int main(int argc, char **argv) try {
    // first argument: matrix name
    auto matrix = std::string("bfwb62.mtx");
    if (argc > 1) {
        std::string arg1(argv[1]);
        matrix = arg1;
    }
    std::cout << "... FILE: " << matrix << std::endl;

    // second argument: verbose mode
    bool verbose = false;
    if (argc > 2) {
        std::string arg2(argv[2]);
        verbose = arg2 == "1";
    }

    // third argument: enforce unsymmetric strategy
    bool enforce_unsymmetric = false;
    if (argc > 3) {
        std::string arg3(argv[3]);
        enforce_unsymmetric = arg3 == "1";
    }

    // message
    if (enforce_unsymmetric) {
        printf("... ENFORCING UNSYMMETRIC STRATEGY\n");
    }

    // read COO
    auto coo = read_matrix_market(matrix);
    printf("... SUCCESS: matrix loaded\n");

    // define the solution vector (x) and right-hand side (rhs)
    auto rhs = std::vector<double>(coo->m, 1.0); // filled with ones
    auto x = std::vector<double>(coo->m, 0.0);

    // default control parameters
    double info[UMFPACK_INFO];
    double control[UMFPACK_CONTROL];
    umfpack_di_defaults(control);
    control[UMFPACK_PRL] = 2.0;

    // enforce unsymmetric strategy
    if (enforce_unsymmetric) {
        control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    }

    // allocate CSC arrays
    int32_t *ap = (int32_t *)malloc((coo->m + 1) * sizeof(int32_t));
    if (ap == NULL) {
        throw "cannot allocate ap array";
    }
    int32_t *ai = (int32_t *)malloc(coo->nnz * sizeof(int32_t));
    if (ai == NULL) {
        free(ap);
        throw "cannot allocate ai array";
    }
    double *ax = (double *)malloc(coo->nnz * sizeof(double));
    if (ax == NULL) {
        free(ai);
        free(ap);
        throw "cannot allocate ax array";
    }
    printf("... SUCCESS: CSC arrays allocated\n");

    // convert COO to CSC
    auto status = umfpack_di_triplet_to_col(coo->m,
                                            coo->m,
                                            coo->nnz,
                                            coo->indices_i.data(),
                                            coo->indices_j.data(),
                                            coo->values_aij.data(),
                                            ap,
                                            ai,
                                            ax,
                                            NULL);
    if (status < 0) {
        umfpack_di_report_status(control, status);
        throw "umfpack_di_triplet_to_col failed";
    }
    printf("... SUCCESS: COO converted to CSC\n");

    // symbolic factorization
    void *symbolic;
    status = umfpack_di_symbolic(coo->m,
                                 coo->m,
                                 ap,
                                 ai,
                                 ax,
                                 &symbolic,
                                 control,
                                 info);
    if (status < 0) {
        umfpack_di_report_info(control, info);
        umfpack_di_report_status(control, status);
        throw "umfpack_di_symbolic failed";
    }
    printf("... SUCCESS: symbolic factorization completed\n");

    // numeric factorization
    void *numeric;
    status = umfpack_di_numeric(ap,
                                ai,
                                ax,
                                symbolic,
                                &numeric,
                                control,
                                info);
    if (status < 0) {
        umfpack_di_report_info(control, info);
        umfpack_di_report_status(control, status);
        umfpack_di_free_symbolic(&symbolic);
        throw "umfpack_di_numeric failed";
    }
    printf("... SUCCESS: numeric factorization completed\n");

    // compute the solution A * x = rhs
    status = umfpack_di_solve(UMFPACK_A,
                              ap,
                              ai,
                              ax,
                              x.data(),
                              rhs.data(),
                              numeric,
                              control,
                              info);
    if (status < 0) {
        umfpack_di_free_symbolic(&symbolic);
        umfpack_di_free_numeric(&numeric);
        throw "umfpack_di_solve failed";
    }
    printf("... SUCCESS: solution calculated\n");

    if (verbose) {
        umfpack_di_report_info(control, info);
    }

    // print residual
    auto r = std::vector<double>(coo->m, 0.0);
    double max_norm_r = resid(coo->m, x.data(), rhs.data(), r.data(), ap, ai, ax);
    printf("... max_norm of residual: %g\n", max_norm_r);

    // check the solution
    const double TOLERANCE = 1e-10;
    auto rhs_new = std::vector<double>(coo->m, 0.0);
    coo->mat_vec_mul(rhs_new, 1.0, x);
    for (size_t k = 0; k < static_cast<size_t>(coo->m); k++) {
        auto diff = fabs(rhs[k] - rhs_new[k]);
        if (diff > TOLERANCE) {
            std::cout << "... ERROR: diff[" << k << "] = " << diff << " is too high" << std::endl;
            return 1;
        }
    }
    printf("... SUCCESS: numerical solution is within tolerance\n");
    return 0;

} catch (std::exception &e) {
    std::cout << e.what() << std::endl;
} catch (std::string &msg) {
    std::cout << msg.c_str() << std::endl;
} catch (const char *msg) {
    std::cout << msg << std::endl;
} catch (...) {
    std::cout << "some unknown exception happened" << std::endl;
}
