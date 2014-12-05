import numpy as np
import cvxpy as cvx


def solve_sdp(G):
    """ Solves the SDP problem:
    
    maximize    1' * s
    subject to  0 <= s <= 1
                [G, G - diag(s); G - diag(s), G] >= 0
    """
    G = np.asarray(G)
    assert G.ndim == 2 and G.shape[0] == G.shape[1]
    p = G.shape[0]
    
    # Define the problem.
    s = cvx.Variable(p)
    objective = cvx.Maximize(cvx.sum_entries(s))
    constraints = [
        0 <= s, s <= 1,
        2*G - cvx.diag(s) == cvx.Semidef(p),
    ]
    prob = cvx.Problem(objective, constraints)
    
    # Solve the problem.
    prob.solve()
    assert prob.status == cvx.OPTIMAL
    
    # Return as array, not as matrix.
    return np.asarray(s.value).flatten()


if __name__ == '__main__':
    import json
    import sys
    
    G = json.load(sys.stdin)
    s = solve_sdp(G)
    json.dump(s.tolist(), sys.stdout)