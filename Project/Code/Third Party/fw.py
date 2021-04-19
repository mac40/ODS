import numpy as np
from scipy import sparse

# .. for plotting ..
import pylab as plt
# .. to generate a synthetic dataset ..
from sklearn import datasets

n_samples, n_features = 10000, 10000
A, b = datasets.make_regression(n_samples, n_features)

def FW(alpha, max_iter=100, tol=1e-8, callback=None):
    # .. initial estimate, could be any feasible point ..
    x_t = sparse.dok_matrix((n_features, 1))

    # .. some quantities can be precomputed ..
    Atb = A.T.dot(b)
    for it in range(max_iter):
        # .. compute gradient. Slightly more involved than usual because ..
        # .. of the use of sparse matrices ..
        Ax = x_t.T.dot(A.T).ravel()
        grad = (A.T.dot(Ax) - Atb)
        print(grad)
        # .. the LMO results in a vector that is zero everywhere except for ..
        # .. a single index. Of this vector we only store its index and magnitude ..
        idx_oracle = np.argmax(np.abs(grad))
        print(idx_oracle)
        mag_oracle = alpha * np.sign(-grad[idx_oracle])
        print(alpha * np.sign(-grad[idx_oracle]))
        d_t = -x_t.copy()
        print(f"dt ->{d_t}")
        d_t[idx_oracle] += mag_oracle
        print(f"dt_idx ->{d_t[idx_oracle]}")
        g_t = - d_t.T.dot(grad).ravel()
        if g_t <= tol:
            break
        q_t = A[:, idx_oracle] * mag_oracle - Ax
        step_size = min(q_t.dot(b - Ax) / q_t.dot(q_t), 1.)
        x_t += step_size * d_t
        if callback is not None:
            callback(g_t)
    return x_t

# .. plot evolution of FW gap ..
trace = []
def callback(g_t):
    trace.append(g_t)

sol = FW(.5 * n_features, callback=callback)
plt.plot(trace / trace[0], lw=3)
plt.yscale('log')
plt.xlabel('Number of iterations')
plt.ylabel('Relative FW gap')
plt.title('FW on a Lasso problem')
plt.xlim((0, 100))
plt.grid()
plt.show()

density = np.mean(sol.toarray().ravel() != 0)
print('Density of solution: %s%%' % (density * 100))