# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False

import numpy as np
from libc.math cimport sqrt
from libc.stdint cimport uint64_t


def decode(const long long[:, ::1] histories, int n):
    cdef Py_ssize_t m = histories.shape[0], W = (n + 63) // 64
    kids_array = np.empty((m, n - 1, 2), dtype=np.int32)
    words_array = np.zeros((m, 2*n - 1, W), dtype=np.uint64)
    cdef int[:, :, ::1] kids = kids_array
    cdef uint64_t[:, :, ::1] words = words_array
    cdef int[::1] active = np.empty(n, dtype=np.int32)
    cdef Py_ssize_t t, i, s, a, b, size
    cdef long long c
    for t in range(m):
        for i in range(n):
            active[i] = i
            words[t, i, i >> 6] = (<uint64_t>1) << (i & 63)
        size = n
        for i in range(n - 1):
            c = histories[t, i]
            b = <Py_ssize_t>((1 + sqrt(1 + 8.0*c))/2)
            while b*(b - 1)//2 > c:
                b -= 1
            while (b + 1)*b//2 <= c:
                b += 1
            a = c - b*(b - 1)//2
            kids[t, i, 0] = active[a]
            kids[t, i, 1] = active[b]
            for s in range(W):
                words[t, n + i, s] = words[t, active[a], s] | words[t, active[b], s]
            for s in range(b, size - 1):
                active[s] = active[s + 1]
            for s in range(a, size - 2):
                active[s] = active[s + 1]
            active[size - 2] = n + i
            size -= 1
    return kids_array, words_array


def rf_matrix(const long long[:, ::1] ids):
    cdef Py_ssize_t m = ids.shape[0], k = ids.shape[1], i, j, p, q, shared
    d_array = np.zeros((m, m))
    cdef double[:, ::1] d = d_array
    for i in range(m):
        for j in range(i + 1, m):
            p = q = shared = 0
            while p < k and q < k:
                if ids[i, p] == ids[j, q]:
                    shared += 1
                    p += 1
                    q += 1
                elif ids[i, p] < ids[j, q]:
                    p += 1
                else:
                    q += 1
            d[i, j] = d[j, i] = 2*(k - shared)
    return d_array


cdef long long findpath(const int[:, ::1] kt, const uint64_t[:, ::1] wt, const int[:, ::1] kr, const uint64_t[:, ::1] wr,
                        int n, int[:, ::1] kids, uint64_t[:, ::1] words, int[::1] parent, int[::1] node_at,
                        int[::1] rank_of) noexcept nogil:
    cdef Py_ssize_t W = wt.shape[1], i, s
    cdef int k, r, x, y, u, v, w, z, other
    cdef bint hit
    cdef long long moves = 0
    for i in range(2*n - 1):
        for s in range(W):
            words[i, s] = wt[i, s]
    for i in range(n):
        rank_of[i] = -1
    for i in range(n - 1):
        kids[i, 0] = kt[i, 0]
        kids[i, 1] = kt[i, 1]
        parent[kt[i, 0]] = n + i
        parent[kt[i, 1]] = n + i
        node_at[i] = n + i
        rank_of[n + i] = i
    for k in range(n - 2):
        x = kr[k, 0]
        y = kr[k, 1]
        if x >= n:
            x = node_at[x - n]
        if y >= n:
            y = node_at[y - n]
        while x != y:
            if rank_of[x] < rank_of[y]:
                x = parent[x]
            else:
                y = parent[y]
        r = rank_of[x]
        while r > k:
            u = node_at[r - 1]
            v = node_at[r]
            if parent[u] == v:
                w = kids[v - n, 0] if kids[v - n, 1] == u else kids[v - n, 1]
                z = kids[u - n, 0]
                other = kids[u - n, 1]
                hit = False
                for s in range(W):
                    if words[z, s] & wr[n + k, s]:
                        hit = True
                        break
                if not hit:
                    z, other = other, z
                kids[u - n, 0] = z
                kids[u - n, 1] = w
                parent[w] = u
                kids[v - n, 0] = u
                kids[v - n, 1] = other
                parent[other] = v
                for s in range(W):
                    words[u, s] = words[z, s] | words[w, s]
            else:
                node_at[r - 1] = v
                node_at[r] = u
                rank_of[v] = r - 1
                rank_of[u] = r
            r -= 1
            moves += 1
    return moves


def rnni_matrix(const int[:, :, ::1] kids, const uint64_t[:, :, ::1] words):
    """FindPath, Collienne & Gavryushkin 2021."""
    cdef Py_ssize_t m = kids.shape[0], i, j
    cdef int n = kids.shape[1] + 1
    d_array = np.zeros((m, m))
    cdef double[:, ::1] d = d_array
    cdef int[:, ::1] work_kids = np.empty((n - 1, 2), dtype=np.int32)
    cdef uint64_t[:, ::1] work_words = np.empty((2*n - 1, words.shape[2]), dtype=np.uint64)
    cdef int[::1] parent = np.empty(2*n - 1, dtype=np.int32)
    cdef int[::1] node_at = np.empty(n - 1, dtype=np.int32)
    cdef int[::1] rank_of = np.empty(2*n - 1, dtype=np.int32)
    for i in range(m):
        for j in range(i + 1, m):
            d[i, j] = d[j, i] = findpath(kids[i], words[i], kids[j], words[j], n, work_kids, work_words, parent,
                                         node_at, rank_of)
    return d_array


cdef long long triplet(const int[:, ::1] ka, const int[:, ::1] kb, int n, int[:, ::1] inter, int[::1] sa,
                       int[::1] sb) noexcept nogil:
    cdef Py_ssize_t x, y, i, j
    cdef int a0, a1, b0, b1
    cdef long long p, agree = 0
    for x in range(n):
        sa[x] = 1
        sb[x] = 1
    for i in range(n - 1):
        sa[n + i] = sa[ka[i, 0]] + sa[ka[i, 1]]
        sb[n + i] = sb[kb[i, 0]] + sb[kb[i, 1]]
    for x in range(n):
        for y in range(n):
            inter[x, y] = 0
        inter[x, x] = 1
        for j in range(n - 1):
            inter[x, n + j] = inter[x, kb[j, 0]] + inter[x, kb[j, 1]]
    for i in range(n - 1):
        a0 = ka[i, 0]
        a1 = ka[i, 1]
        for y in range(2*n - 1):
            inter[n + i, y] = inter[a0, y] + inter[a1, y]
    for i in range(n - 1):
        a0 = ka[i, 0]
        a1 = ka[i, 1]
        for j in range(n - 1):
            b0 = kb[j, 0]
            b1 = kb[j, 1]
            p = <long long>inter[a0, b0]*inter[a1, b1] + <long long>inter[a0, b1]*inter[a1, b0]
            if p:
                agree += p*(n - sa[n + i] - sb[n + j] + inter[n + i, n + j])
    return <long long>n*(n - 1)*(n - 2)//6 - agree


def triplet_matrix(const int[:, :, ::1] kids):
    """Brodal et al. 2013."""
    cdef Py_ssize_t m = kids.shape[0], i, j
    cdef int n = kids.shape[1] + 1
    d_array = np.zeros((m, m))
    cdef double[:, ::1] d = d_array
    cdef int[:, ::1] inter = np.empty((2*n - 1, 2*n - 1), dtype=np.int32)
    cdef int[::1] sa = np.empty(2*n - 1, dtype=np.int32)
    cdef int[::1] sb = np.empty(2*n - 1, dtype=np.int32)
    for i in range(m):
        for j in range(i + 1, m):
            d[i, j] = d[j, i] = triplet(kids[i], kids[j], n, inter, sa, sb)
    return d_array
