import numpy as np


def pairwise_distances(a, b, p=None):
    """
    Compute the (batched) pairwise distance matrix between a and b which both have size [m, n, d] or [n, d]. The result is a tensor of size [m, n, n] (or [n, n]) whose entry [m, i, j] contains the distance_tensor between a[m, i, :] and b[m, j, :].

    Args:
      a : A tensor containing m batches of n points of dimension d. i.e. of size (m, n, d)
      b : A tensor containing m batches of n points of dimension d. i.e. of size (m, n, d)
      p : Norm to use for the distance_tensor

    Returns:
      M : A (m, n, n)-shaped array containing the pairwise distance_tensor between each pair of inputs in a batch.
    """

    squeezed = False
    if len(a.shape) == 2 and len(b.shape) == 2:
        a = a[np.newaxis, :, :]
        b = b[np.newaxis, :, :]
        squeezed = True

    if len(a.shape) != 3:
        raise ValueError("Invalid shape for a. Must be [m, n, d] or [n, d] but got", a.shape)
    if len(b.shape) != 3:
        raise ValueError("Invalid shape for a. Must be [m, n, d] or [n, d] but got", b.shape)

    ret = np.linalg.norm(a[:, :, np.newaxis, :] - b[:, np.newaxis, :, :], axis=-1, ord=p)
    # ret = np.power(np.abs(a[:, :, np.newaxis, :] - b[:, np.newaxis, :, :]), p).sum(3)
    if squeezed:
        ret = np.squeeze(ret)

    return ret


def sinkhorn(a, b, M, eps, max_iters=100, stop_thresh=1e-3):
    """
    Compute the (batched) Sinkhorn correspondences between two dirac delta distributions, U, and V.
    This implementation is numerically stable with float32.

    Args:
      a : A m-sized minibatch of weights for each dirac in the first distribution, U. i.e. shape = (m, n)
      b : A m-sized minibatch of weights for each dirac in the second distribution, V. i.e. shape = (m, n)
      M : A minibatch of n-by-n tensors storing the distance between each pair of diracs in U and V. i.e. shape = (m, n, n) and each i.e. M[k, i, j] = ||u[k,_i] - v[k, j]||
      eps : The reciprocal of the sinkhorn regularization parameter
      max_iters : The maximum number of Sinkhorn iterations
      stop_thresh : Stop if the change in iterates is below this value

    Returns:
      P : An (m, n, n)-shaped array of correspondences between distributions U and V
    """
    # a and b are tensors of size [nb, m] and [nb, n]
    # M is a tensor of size [nb, m, n]

    M = np.squeeze(M)
    a = np.squeeze(a)
    b = np.squeeze(b)
    squeezed = False

    if len(M.shape) == 2 and len(a.shape) == 1 and len(b.shape) == 1:
        M = M[np.newaxis, :, :]
        a = a[np.newaxis, :]
        b = b[np.newaxis, :]
        squeezed = True
    elif len(M.shape) == 2 and len(a.shape) != 1:
        raise ValueError("Invalid shape for a %s, expected [m,] where m is the number of samples in a and "
                         "M has shape [m, n]" % str(a.shape))
    elif len(M.shape) == 2 and len(b.shape) != 1:
        raise ValueError("Invalid shape for a %s, expected [m,] where n is the number of samples in a and "
                         "M has shape [m, n]" % str(b.shape))

    if len(M.shape) != 3:
        raise ValueError("Got unexpected shape for M %s, should be [nb, m, n] where nb is batch size, and "
                         "m and n are the number of samples in the two input measures." % str(M.shape))
    elif len(M.shape) == 3 and len(a.shape) != 2:
        raise ValueError("Invalid shape for a %s, expected [nb, m]  where nb is batch size, m is the number of samples "
                         "in a and M has shape [nb, m, n]" % str(a.shape))
    elif len(M.shape) == 3 and len(b.shape) != 2:
        raise ValueError("Invalid shape for a %s, expected [nb, m]  where nb is batch size, m is the number of samples "
                         "in a and M has shape [nb, m, n]" % str(b.shape))

    nb = M.shape[0]
    m = M.shape[1]
    n = M.shape[2]

    if a.dtype != b.dtype or a.dtype != M.dtype:
        raise ValueError("Tensors a, b, and M must have the same dtype got: dtype(a) = %s, dtype(b) = %s, dtype(M) = %s"
                         % (str(a.dtype), str(b.dtype), str(M.dtype)))
    if a.shape != (nb, m):
        raise ValueError("Got unexpected shape for tensor a (%s). Expected [nb, m] where M has shape [nb, m, n]." %
                         str(a.shape))
    if b.shape != (nb, n):
        raise ValueError("Got unexpected shape for tensor b (%s). Expected [nb, n] where M has shape [nb, m, n]." %
                         str(b.shape))

    u = np.zeros_like(a)
    v = np.zeros_like(b)

    M_t = np.transpose(M, axes=(0, 2, 1))

    def stabilized_log_sum_exp(x):
        max_x = x.max(2)
        x = x - max_x[:, :, np.newaxis]
        ret = np.log(np.sum(np.exp(x), axis=2)) + max_x
        return ret

    for current_iter in range(max_iters):
        u_prev = u
        v_prev = v

        summand_u = (-M + np.expand_dims(v, 1)) / eps
        u = eps * (np.log(a) - stabilized_log_sum_exp(summand_u))

        summand_v = (-M_t + np.expand_dims(u, 1)) / eps
        v = eps * (np.log(b) - stabilized_log_sum_exp(summand_v))

        err_u = np.sum(np.abs(u_prev-u), axis=1).max()
        err_v = np.sum(np.abs(v_prev-v), axis=1).max()

        if err_u < stop_thresh and err_v < stop_thresh:
            break

    log_P = (-M + np.expand_dims(u, 2) + np.expand_dims(v, 1)) / eps

    P = np.exp(log_P)

    if squeezed:
        P = np.squeeze(P)

    return P


def earth_movers_distance(p, q, p_norm=2, eps=1e-4, max_iters=100, stop_thresh=1e-3):
    """
    Compute the (batched) Sinkhorn correspondences between two dirac delta distributions, U, and V.
    This implementation is numerically stable with float32.

    Args:
      p : An (n, d)-shaped array of d-dimensional points
      b : An (m, d)-shaped array of d-dimensional points
      p_norm : Which norm to use. Must be one of {non-zero int, inf, -inf, ‘fro’, ‘nuc’} (default is 2),
      eps : The reciprocal of the sinkhorn regularization parameter (default 1e-4)
      max_iters : The maximum number of Sinkhorn iterations
      stop_thresh : Stop if the change in iterates is below this value

    Returns:
      emd : The earth mover's distance between point clouds p and q
      P : An (n, m)-shaped array of correspondences between point clouds p and q
    """

    M = pairwise_distances(p, q, p_norm)
    a = np.ones(p.shape[0]) / p.shape[0]
    b = np.ones(q.shape[0]) / q.shape[0]
    P = sinkhorn(a, b, M, eps, max_iters, stop_thresh)

    return (P * M).sum(), P
