//
//  History.cpp
//  SINGER
//

#include "ARG.hpp"
#include "Data_reader.hpp"

static int draw_from(const double *w, int n) {
    double total = 0;
    int last = 0;
    for (int i = 0; i < n; i++) {
        total += max(w[i], 0.0);
        if (w[i] > 0) {
            last = i;
        }
    }
    double u = uniform_random()*total;
    int i = 0;
    while (i < last and (w[i] <= 0 or u >= w[i])) {
        u -= max(w[i], 0.0);
        i++;
    }
    return i;
}

static double draw_rate(double shape, double rate) {
    double nu = gamma_distribution<double>(shape, 1/rate)(random_engine);
    if (!(nu > 0)) {
        cerr << "Error: a site rate drawn from Gamma(" << shape << ", " << rate << ") is " << nu << " in double precision (alpha too small for the continuous gamma). " << endl;
        exit(1);
    }
    return nu;
}

static long double log_stay(const Site_model &m, const Join_case &c, int b, double uT) {
    if (m.continuous_gamma()) {
        return m.log_rate_integrated_stay(b, m.alpha - m.rate(b, b)*uT, c.length);
    }
    return (long double) m.rate(b, b)*c.length;
}

static bool is_root_row(const Mutation_row &r) {
    return r.branch.upper_node->index == -1;
}

int ARG::exception_base(int bin, double pos) const {
    const vector<pair<double, int>> &v = exception_sites[bin];
    auto it = lower_bound(v.begin(), v.end(), make_pair(pos, -1));
    return it != v.end() and it->first == pos ? it->second : -1;
}

void ARG::column_below(const Flat_tree &flat, const Tree_order &t, double pos, const vector<char> *forced, vector<int> &below) {
    int nb = (int) flat.parents.size();
    below.assign(nb, 0);
    int site = forced == nullptr and any_missing ? Node::site_index(pos) : -1;
    for (int k = 0; k < nb; k++) {
        Node *l = flat.parents[k].first;
        if (l->is_sample) {
            below[k] = forced != nullptr ? !(*forced)[l->index] : !(any_missing and l->has_missing() and l->missing_at(site, pos));
        }
        if (t.up[k] >= 0) {
            below[t.up[k]] += below[k];
        }
    }
}

unsigned ARG::query_mask(Node *query, const vector<int> &below, double pos, int base) {
    if (query->is_sample) {
        if (any_missing and query->has_missing() and query->missing_at(Node::site_index(pos), pos)) {
            return 0;
        }
        return base >= 0 ? 1u << base : (unsigned) query->get_state(pos);
    }
    return below[query->flat_index] > 0 ? (unsigned) query->get_state(pos) : 0u;
}

Join_case ARG::join_case(const Flat_tree &flat, const Tree_order &t, const vector<int> &below, const Branch &b, double time, Node *query, bool query_data, double unit) {
    Join_case c;
    c.query = query_data;
    bool top_branch = b.upper_node == root.get();
    double t_d = query->time;
    Node *top = flat.parents[t.top].first;
    if (top_branch) {
        c.kind = ROOT;
    } else if (below[b.lower_node->flat_index] > 0) {
        c.kind = E1;
    } else {
        c.kind = E2;
    }
    if (below[t.top] == 0) {
        c.length = unit*((top_branch ? time : top->time) - t_d);
        return c;
    }
    if (c.kind == ROOT) {
        c.anchor = t.top;
        c.length = query_data ? unit*(2*time - top->time - t_d) : unit*(time - top->time);
    } else if (c.kind == E1) {
        c.anchor = b.lower_node->flat_index;
        c.length = unit*(time - t_d);
    } else {
        int w = t.up[b.lower_node->flat_index];
        while (below[w] == 0) {
            w = t.up[w];
        }
        c.anchor = w;
        c.length = unit*(flat.parents[w].first->time - t_d);
    }
    return c;
}

int ARG::join_state(const Join_case &c, const Flat_tree &flat, const vector<Mutation_row> &rows, const Branch &b, double time, double pos) {
    if (c.kind == E2) {
        return node_state(flat.parents[c.anchor].first, pos);
    }
    if (c.kind == ROOT) {
        for (const Mutation_row &r : rows) {
            if (is_root_row(r)) {
                return r.state;
            }
        }
    }
    return branch_state(rows, b, time, pos);
}

int ARG::branch_state(const vector<Mutation_row> &rows, const Branch &b, double time, double pos) {
    int s = node_state(b.upper_node, pos);
    double best = numeric_limits<double>::infinity();
    for (const Mutation_row &r : rows) {
        if (r.branch == b and r.time > time and r.time < best) {
            best = r.time;
            s = r.state;
        }
    }
    return s;
}

double ARG::node_weights(const Branch &b, double time, double t_d, int xu, int xl, unsigned mask, double unit, double *w) const {
    double U[16], L[16], Q[16], total = 0;
    model.transition(unit*(b.upper_node->time - time), U);
    model.transition(unit*(time - b.lower_node->time), L);
    model.transition(unit*(time - t_d), Q);
    for (int s = 0; s < 4; s++) {
        w[s] = U[4*xu + s]*L[4*s + xl]*masked_sum(Q + 4*s, mask);
        total += w[s];
    }
    return total;
}

bool ARG::all_same_state(const Flat_tree &flat, const vector<int> &below, double pos, int base) {
    for (int k = 0; k < (int) flat.parents.size(); k++) {
        Node *n = flat.parents[k].first;
        if (below[k] > 0 and !n->is_sample and node_state(n, pos) != base) {
            return false;
        }
    }
    return true;
}

bool ARG::only_replaced_changes(const vector<Mutation_row> &rows, int base, const Join_case &cn, const Join_case &co, const Branch &bn, const Branch &bo, double pos) {
    bool old_in = co.kind == E1 and !(bo == bn), new_in = cn.kind == E1;
    for (const Mutation_row &r : rows) {
        if (is_root_row(r)) {
            if (r.state != base) {
                return false;
            }
        } else if (!((old_in and r.branch == bo) or (new_in and r.branch == bn))) {
            return false;
        }
    }
    return (!old_in or branch_state(rows, bo, bo.lower_node->time, pos) == base) and (!new_in or branch_state(rows, bn, bn.lower_node->time, pos) == base);
}

const Joint_terms &ARG::joint_terms(const Join_case &cn, const Join_case &co, const Branch &bn, double tn, const Branch &bo, double to, double T, double unit) {
    double fwd[4], rev[4];
    auto lengths = [&](const Join_case &c, const Branch &b, double time, double *l) {
        if (c.kind == E1) {
            l[0] = unit*(b.upper_node->time - time);
            l[1] = unit*(time - b.lower_node->time);
            l[2] = c.length;
        } else {
            l[0] = c.length;
            l[1] = l[2] = 0;
        }
    };
    lengths(cn, bn, tn, fwd);
    lengths(co, bo, to, rev);
    bool same = bo == bn;
    double whole_n = unit*(bn.upper_node->time - bn.lower_node->time);
    fwd[3] = co.kind == E1 and !same ? unit*(bo.upper_node->time - bo.lower_node->time) : 0.0;
    rev[3] = cn.kind == E1 and !same ? whole_n : 0.0;
    double rest = T - fwd[3] - (cn.kind == E1 ? whole_n : 0.0);
    array<double, 10> key = {fwd[0], fwd[1], fwd[2], fwd[3], rev[0], rev[1], rev[2], rev[3], rest, (double) model.stamp};
    auto [it, fresh] = joint_cache.try_emplace(key);
    Joint_terms &j = it->second;
    if (!fresh) {
        return j;
    }
    copy(fwd, fwd + 4, j.fwd);
    j.T = rest;
    bool gamma = model.continuous_gamma();
    long double al = model.alpha, Lf = (long double) fwd[0] + fwd[1] + fwd[2] + fwd[3];
    for (int b = 0; b < 4; b++) {
        long double q = -(long double) model.rate(b, b);
        j.A[b] = (double) (al + q*rest);
        j.log_c[b] = logl(model.pi[b]) + (gamma ? al*(logl(al) - logl(al + q*rest)) : -q*rest);
        j.F[b] = model.node_expectation(b, j.A[b], fwd[0], fwd[1], fwd[2], fwd[3]);
        j.R[b] = model.node_expectation(b, j.A[b], rev[0], rev[1], rev[2], rev[3]);
        j.Fc[b] = gamma ? expl(-al*log1pl(q*Lf/j.A[b])) : expl(-q*Lf);
    }
    return j;
}

long double ARG::joint_log_mix(const Joint_terms &j, const long double *v, int base) const {
    if (base >= 0) {
        return logl(v[base]);
    }
    long double top = *max_element(j.log_c, j.log_c + 4), total = 0;
    for (int b = 0; b < 4; b++) {
        total += expl(j.log_c[b] - top)*v[b];
    }
    return top + logl(total);
}

double ARG::join_value(const Join_case &c, int s, unsigned mask) const {
    if (!c.query) {
        return 1.0;
    }
    if (c.anchor >= 0) {
        return model.transition_dot(c.length, s, leaf_message((int) mask).data());
    }
    double v = 0;
    for (int a = 0; a < 4; a++) {
        if (mask >> a & 1) {
            v += model.pi[a];
        }
    }
    return v;
}

double ARG::pattern_entry(int b, double uT, double l) const {
    if (model.continuous_gamma()) {
        return model.rate_integrated_entry(b, model.alpha - model.rate(b, b)*uT, l);
    }
    return model.transition_entry(l, b, b);
}

double ARG::group_entry(const Join_case &c, int b) const {
    return c.query ? model.pi[b] : 1.0;
}

double ARG::group_constant(const Join_case &c, int b, double uT) const {
    double stay = model.continuous_gamma() ? exp((double) log_stay(model, c, b, uT)) : exp(model.rate(b, b)*c.length);
    if (!c.query) {
        return c.kind == ROOT and c.anchor >= 0 ? stay : 1.0;
    }
    return c.anchor < 0 ? model.pi[b]*stay : stay;
}

double ARG::data_length(const Flat_tree &flat, const Tree_order &t, const vector<int> &below) const {
    double T = 0;
    for (int k = 0; k < (int) flat.parents.size(); k++) {
        if (t.up[k] >= 0 and below[k] > 0) {
            T += flat.lengths[k];
        }
    }
    return T;
}

void ARG::pattern_weights(const Flat_tree &flat, const Tree_order &t, const vector<int> &below, double unit, array<double, 4> &w) const {
    double T = data_length(flat, t, below);
    double x[4], top = -numeric_limits<double>::infinity();
    for (int b = 0; b < 4; b++) {
        if (model.continuous_gamma()) {
            x[b] = model.alpha*log(model.alpha/(model.alpha - model.rate(b, b)*(unit*T)));
        } else {
            x[b] = model.rate(b, b)*unit*T;
        }
        x[b] += below[t.top] > 0 ? log(model.pi[b]) : 0.0;
        top = max(top, x[b]);
    }
    for (int b = 0; b < 4; b++) {
        w[b] = exp(x[b] - top);
    }
}

void ARG::pattern_groups(int bin, vector<Pattern_group> &out) {
    const Bin_gaps &g = gaps_of(bin);
    int n = 0;
    for (Node *s : sample_nodes) {
        n = max(n, s->index + 1);
    }
    double x0 = coordinates[bin], x1 = coordinates[bin + 1];
    out.resize(1);
    Pattern_group &m = out[0];
    m.forced.assign(n, 0);
    m.ranges.clear();
    m.count = g.all_bases;
    m.unknown = g.all_count;
    vector<char> forced(n, 0);
    int gaps = 0;
    size_t used = 1;
    double prev = x0;
    auto segment = [&](double lo, double hi, double count, const array<double, 4> &bases) {
        if (hi <= lo) {
            return;
        }
        if (gaps == 0) {
            out[0].ranges.push_back({lo, hi});
            return;
        }
        if (count <= 0) {
            return;
        }
        if (used == out.size()) {
            out.emplace_back();
        }
        Pattern_group &p = out[used++];
        p.forced = forced;
        p.ranges.assign(1, {lo, hi});
        p.count = bases;
        p.unknown = count;
        for (int b = 0; b < 4; b++) {
            out[0].count[b] -= bases[b];
        }
        out[0].unknown -= count;
    };
    for (size_t j = 0; j < g.sample.size(); j++) {
        if (g.at[j] > prev) {
            segment(prev, g.at[j], g.before[j], g.bases[j]);
            prev = g.at[j];
        }
        forced[g.sample[j]] = g.d[j] < 0;
        gaps -= g.d[j];
    }
    segment(prev, x1, g.tail, g.tail_bases);
    out.resize(used);
    for (auto &[x, b] : exception_sites[bin]) {
        size_t j = 1;
        while (j < out.size() and !(out[j].ranges[0].first <= x and x < out[j].ranges[0].second)) {
            j++;
        }
        Pattern_group &p = j < out.size() ? out[j] : out[0];
        p.count[b] -= 1;
        p.unknown -= 1;
    }
}

void ARG::history_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Node *query, Site_messages &s) {
    vector<int> &below = history_below;
    column_below(flat, t, pos, nullptr, below);
    int base = exception_base(get_index(pos), pos);
    unsigned mask = query_mask(query, below, pos, base);
    s.silent = mask == 0;
    if (s.silent) {
        return;
    }
    if (base >= 0 and all_same_state(flat, below, pos, base)) {
        pattern_fill(flat, t, unit, below, query, s);
        s.base = base;
        return;
    }
    int nb = (int) flat.parents.size();
    s.base = -1;
    s.parts.clear();
    s.unit = unit;
    s.top = t.top;
    s.nu = leaf_message((int) mask);
    s.kind.assign(nb, E2);
    s.upper_state.assign(nb, -1);
    s.lower_state.assign(nb, -1);
    s.fixed.assign(nb, 0.0);
    if (below[t.top] == 0) {
        double v = 0;
        for (int a = 0; a < 4; a++) {
            v += (mask >> a & 1) ? model.pi[a] : 0.0;
        }
        fill(s.fixed.begin(), s.fixed.end(), v);
        return;
    }
    s.rate = column_rate(pos);
    double u = unit*s.rate;
    vector<double> lowest(nb, numeric_limits<double>::infinity());
    int top_state = -1;
    for (const Mutation_row &r : mutations.at(pos)) {
        int k = r.branch.lower_node->flat_index;
        if (is_root_row(r)) {
            top_state = r.state;
        } else if (r.time < lowest[k]) {
            lowest[k] = r.time;
            s.lower_state[k] = (signed char) r.state;
        }
    }
    vector<char> &main = history_main;
    main.assign(nb, 0);
    main[t.top] = 1;
    for (int k = nb - 1; k >= 0; k--) {
        if (k != t.top and t.up[k] >= 0 and main[t.up[k]]) {
            main[k] = 1;
        }
    }
    double t_d = query->time;
    for (int k = 0; k < nb; k++) {
        if (!main[k]) {
            continue;
        }
        if (k == t.top) {
            s.kind[k] = ROOT;
            s.upper_state[k] = (signed char) top_state;
        } else if (below[k] > 0) {
            s.kind[k] = E1;
            s.upper_state[k] = (signed char) node_state(flat.parents[k].second, pos);
            if (s.lower_state[k] < 0) {
                s.lower_state[k] = s.upper_state[k];
            }
            s.fixed[k] = model.transition_entry(u*flat.lengths[k], s.upper_state[k], s.lower_state[k]);
        } else {
            int w = t.up[k];
            while (below[w] == 0) {
                w = t.up[w];
            }
            Node *n = flat.parents[w].first;
            s.fixed[k] = model.transition_dot(u*(n->time - t_d), node_state(n, pos), s.nu.data());
        }
    }
}

void ARG::history_pattern(const Flat_tree &flat, const Tree_order &t, double unit, const vector<char> &forced, Node *query, Site_messages &s) {
    vector<int> &below = history_below;
    column_below(flat, t, 0.0, &forced, below);
    s.silent = query->is_sample ? (bool) forced[query->index] : below[query->flat_index] == 0;
    if (s.silent) {
        return;
    }
    pattern_fill(flat, t, unit, below, query, s);
}

void ARG::pattern_fill(const Flat_tree &flat, const Tree_order &t, double unit, const vector<int> &below, Node *query, Site_messages &s) {
    int nb = (int) flat.parents.size();
    s.parts.clear();
    s.unit = unit;
    s.top = t.top;
    s.kind.assign(nb, E2);
    s.fixed.assign(4*nb, 0.0);
    pattern_weights(flat, t, below, unit, s.weight);
    s.T = unit*data_length(flat, t, below);
    if (below[t.top] == 0) {
        for (int k = 0; k < nb; k++) {
            copy(model.pi.begin(), model.pi.end(), s.fixed.begin() + 4*k);
        }
        return;
    }
    vector<char> &main = history_main;
    main.assign(nb, 0);
    main[t.top] = 1;
    for (int k = nb - 1; k >= 0; k--) {
        if (k != t.top and t.up[k] >= 0 and main[t.up[k]]) {
            main[k] = 1;
        }
    }
    double t_d = query->time, al = model.alpha;
    bool gamma = model.continuous_gamma(), unknown = column_base.empty();
    for (int k = 0; k < nb; k++) {
        if (!main[k]) {
            continue;
        }
        if (k == t.top) {
            s.kind[k] = ROOT;
        } else if (below[k] > 0) {
            s.kind[k] = E1;
            double l0 = unit*flat.lengths[k], total = 0;
            for (int b = 0; b < 4; b++) {
                double q = -model.rate(b, b), D = gamma ? model.rate_integrated_entry(b, al + q*(s.T - l0), l0) : model.transition_entry(l0, b, b);
                if (!unknown) {
                    s.fixed[4*k + b] = D;
                    continue;
                }
                double w = s.weight[b]*(gamma ? exp(al*(log(al + q*s.T) - log(al + q*(s.T - l0)))) : exp(q*l0));
                s.fixed[4*k + b] = w;
                total += w*D;
            }
            for (int b = 0; unknown and b < 4; b++) {
                s.fixed[4*k + b] /= total;
            }
        } else {
            int w = t.up[k];
            while (below[w] == 0) {
                w = t.up[w];
            }
            double L = unit*(flat.parents[w].first->time - t_d);
            for (int b = 0; b < 4; b++) {
                s.fixed[4*k + b] = pattern_entry(b, s.T, L);
            }
        }
    }
}

vector<Mutation_row> &ARG::journal_rows(double pos) {
    if (row_journaled.emplace(pos, 1).second) {
        auto it = mutations.find(pos);
        bool existed = it != mutations.end();
        row_journal.push_back({pos, existed, existed ? it->second : vector<Mutation_row>{}});
    }
    return mutations[pos];
}

void ARG::journal_state(Node *n, double pos, int s) {
    history_journal.push_back({n, pos, n->get_state(pos)});
    if (s < 0) {
        n->write_state(pos, 0);
    } else {
        n->set_site(pos, 1 << s);
    }
}

void ARG::journal_rate(double pos, double v, bool keep) {
    auto it = site_rate.find(pos);
    bool existed = it != site_rate.end();
    rate_journal.push_back({pos, existed, existed ? it->second : 0.0});
    if (keep) {
        site_rate[pos] = v;
    } else {
        site_rate.erase(it);
    }
}

vector<pair<double, int>> &ARG::journal_exceptions(int bin) {
    if (!exception_journaled[bin]) {
        exception_journaled[bin] = 1;
        exception_journal.push_back({bin, exception_sites[bin]});
    }
    return exception_sites[bin];
}

void ARG::restore_histories() {
    for (auto it = history_journal.rbegin(); it != history_journal.rend(); ++it) {
        auto &[n, pos, s] = *it;
        if (s == 0) {
            n->write_state(pos, 0);
        } else {
            n->set_site(pos, s);
        }
    }
    for (auto &[pos, existed, rows] : row_journal) {
        if (existed) {
            mutations[pos] = move(rows);
        } else {
            mutations.erase(pos);
        }
    }
    for (auto &[bin, list] : exception_journal) {
        exception_sites[bin] = move(list);
    }
    for (auto it = rate_journal.rbegin(); it != rate_journal.rend(); ++it) {
        auto &[pos, existed, v] = *it;
        if (existed) {
            site_rate[pos] = v;
        } else {
            site_rate.erase(pos);
        }
    }
    keep_histories();
}

void ARG::keep_histories() {
    for (auto &e : exception_journal) {
        exception_journaled[e.first] = 0;
    }
    exception_journal.clear();
    rate_journal.clear();
    history_journal.clear();
    row_journal.clear();
    row_journaled.clear();
}

void ARG::detach_histories() {
    joint_cache.clear();
    Flat_tree flat;
    Tree_order t;
    flat.assign(start_tree);
    vector<int> &below = history_below;
    auto recomb_it = recombinations.upper_bound(start);
    auto mut_it = mutation_sites.lower_bound(start);
    auto on = removed_branches.begin(), oj = joining_branches.begin();
    Branch ro, bo;
    bool stale = true;
    vector<pair<double, int>> listed;
    for (int i = get_index(start); i < get_index(end); i++) {
        double x = coordinates[i];
        while (on != removed_branches.end() and on->first <= x) {
            ro = on->second;
            ++on;
        }
        while (oj != joining_branches.end() and oj->first <= x) {
            bo = oj->second;
            ++oj;
        }
        if (recomb_it->first == x) {
            flat.forward_update(recomb_it->second);
            ++recomb_it;
            stale = true;
        }
        if (stale) {
            t.assign(flat, root.get());
            stale = false;
        }
        double unit = penalty*thetas[i]/(coordinates[i + 1] - x);
        Node *d = ro.lower_node, *p = ro.upper_node;
        double to = p->time;
        listed.clear();
        while (*mut_it < coordinates[i + 1]) {
            listed.push_back({*mut_it, -1});
            ++mut_it;
        }
        listed.insert(listed.end(), exception_sites[i].begin(), exception_sites[i].end());
        for (auto &[pos, base] : listed) {
            column_below(flat, t, pos, nullptr, below);
            unsigned mask = query_mask(d, below, pos, base);
            Join_case c = join_case(flat, t, below, bo, to, d, mask != 0, unit);
            if (!c.query and (c.kind == E2 or c.anchor < 0)) {
                continue;
            }
            vector<Mutation_row> &rows = journal_rows(pos);
            Node *lo = bo.lower_node, *up = bo.upper_node;
            int xp = node_state(p, pos);
            vector<Branch> drop;
            if (c.query) {
                drop.push_back(Branch(d, p));
            }
            if (c.kind == E1 and c.anchor >= 0) {
                Branch below_p(lo, p), above_p(p, up);
                for (Mutation_row &r : rows) {
                    if (r.branch == below_p or r.branch == above_p) {
                        r.branch = bo;
                    }
                }
            } else if (c.kind == ROOT and c.anchor >= 0) {
                Branch lp(lo, p);
                int xl = xp;
                double lowest = numeric_limits<double>::infinity();
                for (const Mutation_row &r : rows) {
                    if (r.branch == lp and r.time < lowest) {
                        lowest = r.time;
                        xl = r.state;
                    }
                }
                drop.push_back(lp);
                rows.erase(remove_if(rows.begin(), rows.end(), is_root_row), rows.end());
                rows.push_back({Branch(lo, root.get()), xl, 1, lo->time});
            } else {
                vector<Node *> chain;
                if (up != root.get()) {
                    Node *top = c.anchor >= 0 ? flat.parents[c.anchor].first : flat.parents[t.top].first;
                    for (Node *v = up; ; v = flat.parents[v->flat_index].second) {
                        chain.push_back(v);
                        if (v == top) {
                            break;
                        }
                    }
                    reverse(chain.begin(), chain.end());
                }
                chain.push_back(p);
                chain.push_back(d);
                size_t n = chain.size();
                for (size_t j = 0; j + 2 < n; j++) {
                    drop.push_back(Branch(chain[j + 1], chain[j]));
                }
                for (size_t j = c.anchor < 0 ? 0 : 1; j + 2 < n; j++) {
                    journal_state(chain[j], pos, -1);
                }
                if (c.anchor < 0) {
                    rows.erase(remove_if(rows.begin(), rows.end(), is_root_row), rows.end());
                }
            }
            rows.erase(remove_if(rows.begin(), rows.end(), [&](const Mutation_row &r) { return find(drop.begin(), drop.end(), r.branch) != drop.end(); }), rows.end());
            bool back = base >= 0 and none_of(rows.begin(), rows.end(), [](const Mutation_row &r) { return !is_root_row(r); });
            bool subtree = !d->is_sample and below[d->flat_index] > 0;
            if (model.continuous_gamma() and (back or (rows.empty() and !subtree))) {
                journal_rate(pos, 0.0, false);
            }
            if (back) {
                rows.clear();
                for (int k = 0; k < (int) flat.parents.size(); k++) {
                    Node *n = flat.parents[k].first;
                    if (below[k] > 0 and !n->is_sample) {
                        journal_state(n, pos, -1);
                    }
                }
                vector<pair<double, int>> &ex = journal_exceptions(i);
                ex.erase(lower_bound(ex.begin(), ex.end(), make_pair(pos, -1)));
            }
            if (rows.empty()) {
                mutations.erase(pos);
            }
        }
    }
}

void ARG::attach_histories(map<double, Branch> &joins, map<double, Branch> &added) {
    Flat_tree flat;
    Tree_order t;
    flat.assign(start_tree);
    vector<int> &below = history_below;
    vector<Pattern_group> &groups = history_groups;
    auto recomb_it = recombinations.upper_bound(start);
    auto mut_it = mutation_sites.lower_bound(start);
    auto nn = added.begin(), nj = joins.begin(), on = removed_branches.begin(), oj = joining_branches.begin();
    Branch rn, bn, ro, bo;
    bool stale = true, unknown = column_base.empty();
    vector<pair<double, int>> path, path_u, path_l, listed, replaced;
    vector<Node *> chain;
    vector<double> chosen;
    double unit = 0, unit_bin = 0, tn = 0, to = 0, pos = 0;
    Node *d = nullptr, *pn = nullptr;
    int bin = 0;
    auto place_on = [&](vector<Mutation_row> &rows, const Branch &b) {
        for (auto &[tau, next] : path) {
            rows.push_back({b, next, 0, b.upper_node->time - tau/unit});
        }
    };
    auto drop_branch = [&](vector<Mutation_row> &rows, const Branch &b) {
        rows.erase(remove_if(rows.begin(), rows.end(), [&](const Mutation_row &r) { return r.branch == b; }), rows.end());
    };
    auto place_down = [&](vector<Mutation_row> &rows, int s, size_t from, double offset, double h0) {
        int cur = s;
        size_t j = 0;
        for (size_t i = from; i < path.size(); i++) {
            double h = h0 - (path[i].first - offset)/unit;
            while (h < chain[j + 1]->time) {
                j++;
                journal_state(chain[j], pos, cur);
            }
            rows.push_back({Branch(chain[j + 1], chain[j]), path[i].second, 0, h});
            cur = path[i].second;
        }
        while (j + 2 < chain.size()) {
            j++;
            journal_state(chain[j], pos, cur);
        }
    };
    auto place_up = [&](vector<Mutation_row> &rows, int s, size_t count) {
        Node *l = bn.lower_node;
        int cur = s;
        for (size_t i = 0; i < count; i++) {
            rows.push_back({Branch(l, pn), cur, 0, l->time + path[i].first/unit});
            cur = path[i].second;
        }
        return cur;
    };
    auto chain_from = [&](Node *top) {
        chain.clear();
        if (top != pn) {
            for (Node *v = bn.upper_node; ; v = flat.parents[v->flat_index].second) {
                chain.push_back(v);
                if (v == top) {
                    break;
                }
            }
            reverse(chain.begin(), chain.end());
        }
        chain.push_back(pn);
        chain.push_back(d);
    };
    auto drop_root_row = [&](vector<Mutation_row> &rows) {
        rows.erase(remove_if(rows.begin(), rows.end(), is_root_row), rows.end());
    };
    auto split = [&](vector<Mutation_row> &rows) {
        for (Mutation_row &r : rows) {
            if (r.branch == bn) {
                r.branch = r.time > tn ? Branch(pn, bn.upper_node) : Branch(bn.lower_node, pn);
            }
        }
    };
    auto join_root = [&](vector<Mutation_row> &rows, int s) {
        Node *l = bn.lower_node;
        double up = unit*(tn - l->time);
        size_t m = 0;
        while (m < path.size() and path[m].first < up) {
            m++;
        }
        drop_root_row(rows);
        int x = place_up(rows, s, m);
        journal_state(pn, pos, x);
        chain_from(pn);
        place_down(rows, x, m, up, tn);
        rows.push_back({Branch(pn, root.get()), x, 1, tn});
    };
    auto join_top = [&](vector<Mutation_row> &rows, int x) {
        Node *top = bn.upper_node == root.get() ? pn : flat.parents[t.top].first;
        chain_from(top);
        journal_state(top, pos, x);
        place_down(rows, x, 0, 0.0, top->time);
        rows.push_back({Branch(top, root.get()), x, 1, top->time});
    };
    auto free_root = [&](vector<Mutation_row> &rows, int s, int end, double L, bool change) {
        model.sample_path(s, end, L, model.transition_entry(L, s, end), path, change);
        drop_root_row(rows);
        int x = place_up(rows, s, path.size());
        journal_state(pn, pos, x);
        rows.push_back({Branch(pn, root.get()), x, 1, tn});
    };
    auto mark_states = [&](int b) {
        for (int k = 0; k < (int) flat.parents.size(); k++) {
            Node *n = flat.parents[k].first;
            if (below[k] > 0 and !n->is_sample) {
                journal_state(n, pos, b);
            }
        }
    };
    auto add_exception = [&](int b) {
        vector<pair<double, int>> &ex = journal_exceptions(bin);
        ex.insert(lower_bound(ex.begin(), ex.end(), make_pair(pos, -1)), {pos, b});
    };
    auto exception = [&](int b, const Join_case &c, bool change) {
        vector<Mutation_row> &rows = journal_rows(pos);
        mark_states(b);
        double L = c.length, stay = exp(model.rate(b, b)*L), w[4];
        Node *top = flat.parents[t.top].first;
        if (!c.query or c.anchor < 0) {
            for (int s = 0; s < 4; s++) {
                w[s] = model.transition_entry(L, b, s) - (s == b ? stay : 0.0);
            }
            int x = draw_from(w, 4);
            if (!c.query) {
                free_root(rows, b, x, L, x == b);
            } else {
                model.sample_path(x, b, L, model.transition_entry(L, x, b), path, x == b);
                join_top(rows, x);
            }
        } else {
            if (change) {
                model.sample_path(b, b, L, 0.0, path, true);
            } else {
                path.clear();
            }
            if (c.kind == ROOT) {
                join_root(rows, b);
            } else {
                chain_from(flat.parents[c.anchor].first);
                place_down(rows, b, 0, 0.0, chain[0]->time);
                rows.push_back({Branch(top, root.get()), b, 1, top->time});
            }
        }
        add_exception(b);
    };
    auto one_entry_exception = [&](int b, bool data, double uT) {
        if (model.continuous_gamma()) {
            double A = model.alpha - model.rate(b, b)*uT;
            double nu = model.sample_change_rate(b, A, join_case(flat, t, below, bn, tn, d, data, unit_bin).length);
            journal_rate(pos, nu, true);
            unit = unit_bin*nu;
        }
        exception(b, join_case(flat, t, below, bn, tn, d, data, unit), true);
        unit = unit_bin;
    };
    auto node_place = [&](vector<Mutation_row> &rows, int s) {
        Node *ln = bn.lower_node, *un = bn.upper_node;
        for (auto &[tau, next] : path_u) {
            rows.push_back({Branch(pn, un), next, 0, un->time - tau/unit});
        }
        for (auto &[tau, next] : path_l) {
            rows.push_back({Branch(ln, pn), next, 0, tn - tau/unit});
        }
        journal_state(pn, pos, s);
        chain_from(pn);
        place_down(rows, s, 0, 0.0, tn);
    };
    auto old_whole = [&](vector<Mutation_row> &rows, int xu, int xl, bool change) {
        double L = unit*(bo.upper_node->time - bo.lower_node->time);
        drop_branch(rows, bo);
        model.sample_path(xu, xl, L, model.transition_entry(L, xu, xl), path, change);
        place_on(rows, bo);
    };
    auto node_join = [&](vector<Mutation_row> &rows, unsigned mask) {
        Node *ln = bn.lower_node, *un = bn.upper_node;
        int xu = node_state(un, pos), xl = branch_state(rows, bn, ln->time, pos);
        double Lu = unit*(un->time - tn), Ll = unit*(tn - ln->time), Lq = unit*(tn - d->time), w[4];
        drop_branch(rows, bn);
        node_weights(bn, tn, d->time, xu, xl, mask, unit, w);
        int s = draw_from(w, 4);
        for (int a = 0; a < 4; a++) {
            w[a] = (mask >> a & 1) ? model.transition_entry(Lq, s, a) : 0.0;
        }
        int a = draw_from(w, 4);
        model.sample_path(xu, s, Lu, model.transition_entry(Lu, xu, s), path_u);
        model.sample_path(s, xl, Ll, model.transition_entry(Ll, s, xl), path_l);
        model.sample_path(s, a, Lq, w[a], path);
        node_place(rows, s);
    };
    auto group_column = [&](int b, const Join_case &c, bool whole, bool change_n, bool change_o) {
        vector<Mutation_row> &rows = journal_rows(pos);
        if (whole and change_o) {
            old_whole(rows, b, b, true);
        }
        if (c.kind != E1) {
            exception(b, c, change_n);
            return;
        }
        mark_states(b);
        double L[3] = {unit*(bn.upper_node->time - tn), unit*(tn - bn.lower_node->time), unit*(tn - d->time)};
        vector<pair<double, int>> *out[3] = {&path_u, &path_l, &path};
        int s = b;
        for (auto *p : out) {
            p->clear();
        }
        if (change_n) {
            double w[4];
            node_weights(bn, tn, d->time, b, b, 1u << b, unit, w);
            w[b] -= exp(model.rate(b, b)*(L[0] + L[1] + L[2]));
            s = draw_from(w, 4);
            int ends[3][2] = {{b, s}, {s, b}, {s, b}};
            double P[3], stay[3], pattern[8] = {0};
            for (int i = 0; i < 3; i++) {
                P[i] = model.transition_entry(L[i], ends[i][0], ends[i][1]);
                stay[i] = exp(model.rate(b, b)*L[i])/P[i];
            }
            for (int m = 1; s == b and m < 8; m++) {
                pattern[m] = 1;
                for (int i = 0; i < 3; i++) {
                    pattern[m] *= m >> i & 1 ? 1 - stay[i] : stay[i];
                }
            }
            int m = s == b ? draw_from(pattern, 8) : 0;
            for (int i = 0; i < 3; i++) {
                if (s != b or m >> i & 1) {
                    model.sample_path(ends[i][0], ends[i][1], L[i], P[i], *out[i], s == b);
                }
            }
        }
        node_place(rows, s);
        Node *top = flat.parents[t.top].first;
        rows.push_back({Branch(top, root.get()), b, 1, top->time});
        add_exception(b);
    };
    auto choose = [&](const Pattern_group &g, int b, long long n, long long k) {
        vector<long long> pick;
        set<long long> taken;
        for (long long j = n - k; j < n; j++) {
            long long r = min(j, (long long) (uniform_random()*(j + 1)));
            taken.insert(taken.count(r) ? j : r);
        }
        pick.assign(taken.begin(), taken.end());
        chosen.clear();
        const vector<pair<double, int>> &ex = exception_sites[bin];
        long long index = 0;
        size_t next = 0;
        for (auto &range : g.ranges) {
            for (double x = ceil(range.first); x < range.second and next < pick.size(); x += 1) {
                auto e = lower_bound(ex.begin(), ex.end(), make_pair(x, -1));
                if (e != ex.end() and e->first == x) {
                    continue;
                }
                if (!unknown) {
                    if (column_base[(size_t) x] != b) {
                        continue;
                    }
                } else {
                    if (mutation_sites.count(x) or binary_search(unassayed_sites.begin(), unassayed_sites.end(), x)) {
                        continue;
                    }
                    auto mk = upper_bound(masked.begin(), masked.end(), make_pair(x, numeric_limits<double>::infinity()));
                    if (mk != masked.begin() and prev(mk)->second > x) {
                        continue;
                    }
                }
                if (index == pick[next]) {
                    chosen.push_back(x);
                    next++;
                }
                index++;
            }
        }
    };
    auto joint_column = [&](const Joint_terms &J, const Join_case &c, int b) {
        double w[4];
        int beta = b;
        if (b < 0) {
            long double top = *max_element(J.log_c, J.log_c + 4);
            for (int a = 0; a < 4; a++) {
                w[a] = (double) (expl(J.log_c[a] - top)*(J.F[a] - J.Fc[a]));
            }
            beta = draw_from(w, 4);
        }
        bool gamma = model.continuous_gamma(), whole = J.fwd[3] > 0;
        double nu = gamma ? model.node_rate(beta, J.A[beta], J.fwd[0], J.fwd[1], J.fwd[2], J.fwd[3]) : 1.0;
        unit = unit_bin*nu;
        Join_case cv = join_case(flat, t, below, bn, tn, d, true, unit);
        double value = c.kind == E1 ? node_weights(bn, tn, d->time, beta, beta, 1u << beta, unit, w) : model.transition_entry(nu*J.fwd[0], beta, beta);
        double xn = model.node_excess(beta, nu*J.fwd[0], nu*J.fwd[1], nu*J.fwd[2], 0.0)/value, x0 = nu*J.fwd[3];
        double xo = whole ? x0*x0*model.h(beta, x0)/model.transition_entry(x0, beta, beta) : 0.0;
        double pattern[3] = {xn*(1 - xo), (1 - xn)*xo, xn*xo};
        int m = draw_from(pattern, 3);
        group_column(beta, cv, whole, m != 1, m != 0);
        if (gamma) {
            journal_rate(pos, nu, true);
        }
        unit = unit_bin;
    };
    auto joint_group = [&](const Joint_terms &J, const Join_case &c, const Pattern_group *g, int b, long long n) {
        long double log_z = joint_log_mix(J, J.F, b), log_c = joint_log_mix(J, J.Fc, b);
        long long k = binomial_random(n, max(0.0, 1 - exp((double) (log_c - log_z))));
        if (k == 0) {
            return;
        }
        if (g == nullptr) {
            chosen.assign(1, pos);
        } else {
            choose(*g, b, n, k);
        }
        for (double y : chosen) {
            pos = y;
            joint_column(J, c, b);
        }
    };
    for (bin = get_index(start); bin < get_index(end); bin++) {
        double x = coordinates[bin];
        while (nn != added.end() and nn->first <= x) {
            rn = nn->second;
            ++nn;
        }
        while (nj != joins.end() and nj->first <= x) {
            bn = nj->second;
            ++nj;
        }
        while (on != removed_branches.end() and on->first <= x) {
            ro = on->second;
            ++on;
        }
        while (oj != joining_branches.end() and oj->first <= x) {
            bo = oj->second;
            ++oj;
        }
        if (recomb_it->first == x) {
            flat.forward_update(recomb_it->second);
            ++recomb_it;
            stale = true;
        }
        if (stale) {
            t.assign(flat, root.get());
            stale = false;
        }
        unit_bin = penalty*thetas[bin]/(coordinates[bin + 1] - x);
        unit = unit_bin;
        d = rn.lower_node;
        pn = rn.upper_node;
        tn = pn->time;
        to = ro.upper_node->time;
        replaced.clear();
        listed.clear();
        while (*mut_it < coordinates[bin + 1]) {
            listed.push_back({*mut_it, -1});
            ++mut_it;
        }
        listed.insert(listed.end(), exception_sites[bin].begin(), exception_sites[bin].end());
        for (auto &[p, base] : listed) {
            pos = p;
            column_below(flat, t, pos, nullptr, below);
            unsigned mask = query_mask(d, below, pos, base);
            Join_case c = join_case(flat, t, below, bn, tn, d, mask != 0, unit_bin);
            if (!c.query and (c.kind == E2 or c.anchor < 0)) {
                continue;
            }
            Join_case co = join_case(flat, t, below, bo, to, d, mask != 0, unit_bin);
            if (base >= 0 and c.query and c.anchor >= 0 and only_replaced_changes(mutations.at(pos), base, c, co, bn, bo, pos)) {
                replaced.push_back({pos, base});
                continue;
            }
            double rate = 1;
            if (model.continuous_gamma()) {
                auto it = site_rate.find(pos);
                if (it != site_rate.end()) {
                    rate = it->second;
                } else {
                    rate = draw_rate(model.alpha, model.alpha);
                    journal_rate(pos, rate, true);
                }
            }
            unit = unit_bin*rate;
            c = join_case(flat, t, below, bn, tn, d, mask != 0, unit);
            vector<Mutation_row> &rows = journal_rows(pos);
            if (c.query and co.kind == E1 and co.anchor >= 0 and !(bo == bn)) {
                int xu = node_state(bo.upper_node, pos), xl = branch_state(rows, bo, bo.lower_node->time, pos);
                old_whole(rows, xu, xl, false);
            }
            int s = c.anchor >= 0 ? join_state(c, flat, rows, bn, tn, pos) : -1;
            double w[4];
            bool node = c.query and c.kind == E1 and c.anchor >= 0;
            if (node) {
                node_join(rows, mask);
            } else if (!c.query and c.kind == E1) {
                split(rows);
                journal_state(pn, pos, s);
            } else if (!c.query) {
                for (int a = 0; a < 4; a++) {
                    w[a] = model.transition_entry(c.length, s, a);
                }
                free_root(rows, s, draw_from(w, 4), c.length, false);
            } else if (c.anchor < 0) {
                for (int a = 0; a < 4; a++) {
                    w[a] = (mask >> a & 1) ? model.pi[a] : 0.0;
                }
                int a = draw_from(w, 4);
                for (int y = 0; y < 4; y++) {
                    w[y] = model.transition_entry(c.length, a, y);
                }
                int y = draw_from(w, 4);
                model.sample_path(y, a, c.length, model.transition_entry(c.length, y, a), path);
                join_top(rows, y);
            } else {
                for (int a = 0; a < 4; a++) {
                    w[a] = (mask >> a & 1) ? model.transition_entry(c.length, s, a) : 0.0;
                }
                int a = draw_from(w, 4);
                model.sample_path(s, a, c.length, w[a], path);
                if (c.kind == ROOT) {
                    join_root(rows, s);
                } else {
                    chain_from(flat.parents[c.anchor].first);
                    place_down(rows, s, 0, 0.0, chain[0]->time);
                }
            }
        }
        unit = unit_bin;
        pattern_groups(bin, groups);
        for (Pattern_group &g : groups) {
            column_below(flat, t, x, &g.forced, below);
            bool data = d->is_sample ? !g.forced[d->index] : below[d->flat_index] > 0;
            Join_case c = join_case(flat, t, below, bn, tn, d, data, unit);
            if (!c.query and !(c.kind == ROOT and c.anchor >= 0)) {
                continue;
            }
            double uT = unit*data_length(flat, t, below);
            if (c.query and c.anchor >= 0) {
                const Joint_terms &J = joint_terms(c, join_case(flat, t, below, bo, to, d, true, unit), bn, tn, bo, to, uT, unit);
                for (int b = unknown ? -1 : 0; b < (unknown ? 0 : 4); b++) {
                    long long n = llround(b < 0 ? g.unknown : g.count[b]);
                    if (n == 0) {
                        continue;
                    }
                    joint_group(J, c, &g, b, n);
                }
            } else if (unknown) {
                long long n = llround(g.unknown);
                if (n == 0) {
                    continue;
                }
                array<double, 4> cw;
                pattern_weights(flat, t, below, unit, cw);
                double e[4], cs[4], te = 0, tc = 0, w[4];
                for (int b = 0; b < 4; b++) {
                    e[b] = group_entry(c, b);
                    cs[b] = group_constant(c, b, uT);
                    te += cw[b]*e[b];
                    tc += cw[b]*cs[b];
                }
                long long k = binomial_random(n, max(0.0, 1 - tc/te));
                if (k == 0) {
                    continue;
                }
                choose(g, -1, n, k);
                for (double y : chosen) {
                    pos = y;
                    for (int b = 0; b < 4; b++) {
                        w[b] = cw[b]*(e[b] - cs[b]);
                    }
                    int b = draw_from(w, 4);
                    one_entry_exception(b, data, uT);
                }
            } else {
                for (int b = 0; b < 4; b++) {
                    long long n = llround(g.count[b]);
                    if (n == 0) {
                        continue;
                    }
                    double e = group_entry(c, b), cs = group_constant(c, b, uT);
                    long long k = binomial_random(n, max(0.0, 1 - cs/e));
                    if (k == 0) {
                        continue;
                    }
                    choose(g, b, n, k);
                    for (double y : chosen) {
                        pos = y;
                        one_entry_exception(b, data, uT);
                    }
                }
            }
        }
        for (auto &[p, base] : replaced) {
            pos = p;
            column_below(flat, t, pos, nullptr, below);
            Join_case c = join_case(flat, t, below, bn, tn, d, true, unit), co = join_case(flat, t, below, bo, to, d, true, unit);
            const Joint_terms &J = joint_terms(c, co, bn, tn, bo, to, unit*data_length(flat, t, below), unit);
            int b = unknown ? -1 : base;
            vector<Mutation_row> &rows = journal_rows(pos);
            rows.clear();
            mutations.erase(pos);
            mark_states(-1);
            if (model.continuous_gamma()) {
                journal_rate(pos, 0.0, false);
            }
            vector<pair<double, int>> &ex = journal_exceptions(bin);
            ex.erase(lower_bound(ex.begin(), ex.end(), make_pair(pos, -1)));
            joint_group(J, c, nullptr, b, 1);
        }
    }
}

void ARG::write_histories(string filename, const map<double, vector<Mutation_row>> &state, const function<double(Node *)> *time_of) {
    ofstream file(filename);
    file.precision(numeric_limits<double>::max_digits10);
    auto time_at = [&](Node *n) { return time_of == nullptr ? n->time : (*time_of)(n); };
    Flat_tree flat;
    Tree_order t;
    vector<int> below;
    auto r_it = recombinations.begin();
    bool stale = true;
    for (auto &[pos, rows] : state) {
        while (r_it != recombinations.end() and r_it->first <= pos) {
            flat.forward_update(r_it->second);
            ++r_it;
            stale = true;
        }
        if (stale) {
            t.assign(flat, root.get());
            stale = false;
        }
        column_below(flat, t, pos, nullptr, below);
        int total = below[t.top];
        for (const Mutation_row &r : rows) {
            Node *l = r.branch.lower_node, *u = r.branch.upper_node;
            int informed = 1;
            double time = time_at(l);
            if (!is_root_row(r)) {
                int k = l->flat_index;
                informed = below[k] > 0 and below[k] < total;
                time = time_of == nullptr ? r.time : time_at(l) + (r.time - l->time)/(u->time - l->time)*(time_at(u) - time_at(l));
            }
            file << pos << " " << l->index << " " << u->index << " " << r.state << " " << informed << " " << (time + time_offset)*Ne << "\n";
        }
    }
}

void ARG::redraw_histories(map<double, vector<Mutation_row>> &out, bool state) {
    if (!state) {
        out.clear();
    }
    Flat_tree flat;
    Tree_order t;
    Peel p;
    vector<int> st, below;
    vector<pair<double, int>> path;
    auto r_it = recombinations.begin();
    auto mut_it = mutation_sites.lower_bound(0);
    bool stale = true;
    for (int i = 0; i < bin_num; i++) {
        double x0 = coordinates[i], x1 = coordinates[i + 1];
        while (r_it != recombinations.end() and r_it->first <= x0) {
            flat.forward_update(r_it->second);
            ++r_it;
            stale = true;
        }
        if (stale) {
            t.assign(flat, root.get());
            stale = false;
        }
        double unit = penalty*thetas[i]/(x1 - x0);
        if (!state) {
            for (auto &e : exception_sites[i]) {
                out[e.first] = mutations.at(e.first);
            }
        }
        for (; *mut_it < x1; ++mut_it) {
            double pos = *mut_it;
            column_below(flat, t, pos, nullptr, below);
            if (below[t.top] == 0) {
                continue;
            }
            double u = unit*column_rate(pos);
            peel_site(flat, t, u, pos, p);
            vector<Mutation_row> &rows = out[pos];
            rows.clear();
            draw_states(flat, t, u, p, st);
            const double *P = edge_transitions(model, flat, u, p);
            for (int k = 0; k < (int) flat.parents.size(); k++) {
                if (below[k] == 0) {
                    continue;
                }
                Node *lower = flat.parents[k].first, *upper = flat.parents[k].second;
                if (k == t.top) {
                    rows.push_back({Branch(lower, upper), st[k], 1, lower->time});
                    continue;
                }
                model.sample_path(st[t.up[k]], st[k], u*flat.lengths[k], P[16*k + 4*st[t.up[k]] + st[k]], path);
                for (auto &[tau, next] : path) {
                    rows.push_back({Branch(lower, upper), next, 0, upper->time - tau/u});
                }
            }
            if (!state) {
                continue;
            }
            for (int k = 0; k < (int) flat.parents.size(); k++) {
                Node *n = flat.parents[k].first;
                if (below[k] > 0 and !n->is_sample) {
                    n->set_site(pos, 1 << st[k]);
                }
            }
        }
    }
    if (state and model.continuous_gamma()) {
        History_statistics s;
        history_statistics(s);
        redraw_rates(s);
    }
}

void ARG::history_statistics(History_statistics &s) {
    s = History_statistics();
    Flat_tree flat;
    Tree_order t;
    vector<int> below;
    vector<Pattern_group> &groups = history_groups;
    vector<tuple<int, double, int>> pieces;
    auto r_it = recombinations.begin();
    auto mut_it = mutation_sites.lower_bound(0);
    bool stale = true;
    for (int i = 0; i < bin_num; i++) {
        double x0 = coordinates[i], x1 = coordinates[i + 1];
        while (r_it != recombinations.end() and r_it->first <= x0) {
            flat.forward_update(r_it->second);
            ++r_it;
            stale = true;
        }
        if (stale) {
            t.assign(flat, root.get());
            stale = false;
        }
        double unit = penalty*thetas[i]/(x1 - x0);
        pattern_groups(i, groups);
        for (Pattern_group &g : groups) {
            column_below(flat, t, x0, &g.forced, below);
            if (below[t.top] == 0) {
                continue;
            }
            array<double, 5> &n = s.groups[unit*data_length(flat, t, below)];
            for (int b = 0; b < 4; b++) {
                n[b] += g.count[b];
            }
            n[4] += g.unknown;
        }
        auto column = [&](double pos) {
            column_below(flat, t, pos, nullptr, below);
            if (below[t.top] == 0) {
                return;
            }
            Column_statistics c;
            c.pos = pos;
            array<double, 4> credit = {0, 0, 0, 0};
            pieces.clear();
            for (const Mutation_row &r : mutations.at(pos)) {
                if (is_root_row(r)) {
                    s.root[r.state] += 1;
                } else {
                    pieces.push_back({r.branch.lower_node->flat_index, -r.time, r.state});
                }
            }
            sort(pieces.begin(), pieces.end());
            for (int k = 0; k < (int) flat.parents.size(); k++) {
                if (t.up[k] >= 0 and below[k] > 0) {
                    credit[node_state(flat.parents[k].second, pos)] += flat.lengths[k];
                }
            }
            int k = -1, x = -1;
            for (auto &[j, minus_time, next] : pieces) {
                if (j != k) {
                    k = j;
                    x = node_state(flat.parents[k].second, pos);
                }
                double below_change = -minus_time - flat.parents[k].first->time;
                credit[x] -= below_change;
                credit[next] += below_change;
                s.changes[x][next] += 1;
                x = next;
            }
            c.changes = (int) pieces.size();
            double rate = column_rate(pos);
            for (int a = 0; a < 4; a++) {
                c.time[a] = unit*credit[a];
                s.time[a] += rate*c.time[a];
            }
            s.columns.push_back(c);
        };
        while (*mut_it < x1) {
            column(*mut_it);
            ++mut_it;
        }
        for (auto &e : exception_sites[i]) {
            column(e.first);
        }
    }
}

static double group_statistics_term(const Site_model &m, const History_statistics &s, bool unknown) {
    double v = 0;
    for (auto &[uT, n] : s.groups) {
        double term[4];
        for (int b = 0; b < 4; b++) {
            term[b] = m.continuous_gamma() ? m.alpha*log(m.alpha/(m.alpha - m.rate(b, b)*uT)) : m.rate(b, b)*uT;
        }
        if (unknown) {
            double total = 0;
            for (int b = 0; b < 4; b++) {
                total += m.pi[b]*exp(term[b]);
            }
            v += n[4]*log(total);
        } else {
            for (int b = 0; b < 4; b++) {
                if (n[b] > 0) {
                    v += n[b]*(log(m.pi[b]) + term[b]);
                }
            }
        }
    }
    return v;
}

double ARG::statistics_log_density(const History_statistics &s) const {
    double v = group_statistics_term(model, s, column_base.empty());
    for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
            if (a != b and s.changes[a][b] > 0) {
                v += s.changes[a][b]*log(model.rate(a, b));
            }
        }
        v += model.rate(a, a)*s.time[a] + s.root[a]*log(model.pi[a]);
    }
    return v;
}

double ARG::alpha_log_density(const History_statistics &s) const {
    double v = group_statistics_term(model, s, column_base.empty()), al = model.alpha;
    double shared = al*log(al) - lgamma(al);
    for (const Column_statistics &c : s.columns) {
        double T = 0;
        for (int a = 0; a < 4; a++) {
            T -= model.rate(a, a)*c.time[a];
        }
        v += shared + lgamma(al + c.changes) - (al + c.changes)*log(al + T);
    }
    return v;
}

void ARG::redraw_rates(const History_statistics &s) {
    for (const Column_statistics &c : s.columns) {
        double T = 0;
        for (int a = 0; a < 4; a++) {
            T -= model.rate(a, a)*c.time[a];
        }
        site_rate[c.pos] = draw_rate(model.alpha + c.changes, model.alpha + T);
    }
}

void ARG::write_rates(string filename) {
    ofstream file(filename);
    file.precision(numeric_limits<double>::max_digits10);
    for (auto &[pos, rate] : site_rate) {
        file << pos << " " << rate << "\n";
    }
}

bool ARG::read_histories(string muts_file, string rates_file, string rescaled_file) {
    ifstream muts(muts_file), rates(rates_file);
    if (!muts.good() or (model.continuous_gamma() and !rates.good())) {
        return false;
    }
    vector<double> rescaled;
    if (!rescaled_file.empty()) {
        rescaled = Data_reader::read_numbers(rescaled_file, "rescaled node file not found");
        for (double &x : rescaled) {
            x = x/Ne - time_offset;
        }
    }
    map<double, vector<Mutation_row>> read;
    double pos, time;
    int l, u, state, informed;
    while (muts >> pos >> l >> u >> state >> informed >> time) {
        Node *lower = node_list[l], *upper = u < 0 ? root.get() : node_list[u];
        if (upper == root.get()) {
            read[pos].push_back({Branch(lower, upper), state, 1, lower->time});
            continue;
        }
        double w = time/Ne - time_offset;
        if (!rescaled.empty()) {
            w = lower->time + (w - rescaled[l])/(rescaled[u] - rescaled[l])*(upper->time - lower->time);
        }
        read[pos].push_back({Branch(lower, upper), state, 0, w});
    }
    mutations.clear();
    exception_sites.assign(bin_num, {});
    exception_journaled.assign(bin_num, 0);
    site_rate.clear();
    Flat_tree flat;
    Tree_order t;
    vector<int> below, st;
    auto r_it = recombinations.begin();
    bool stale = true;
    for (auto &[x, rows] : read) {
        while (r_it != recombinations.end() and r_it->first <= x) {
            flat.forward_update(r_it->second);
            ++r_it;
            stale = true;
        }
        if (stale) {
            t.assign(flat, root.get());
            stale = false;
        }
        column_below(flat, t, x, nullptr, below);
        int nb = (int) flat.parents.size(), leaf = -1;
        st.assign(nb, -1);
        for (const Mutation_row &r : rows) {
            if (is_root_row(r)) {
                st[t.top] = r.state;
            }
        }
        for (int k = nb - 1; k >= 0; k--) {
            if (k == t.top or below[k] == 0) {
                continue;
            }
            double lowest = numeric_limits<double>::infinity();
            st[k] = st[t.up[k]];
            for (const Mutation_row &r : rows) {
                if (!is_root_row(r) and r.branch.lower_node == flat.parents[k].first and r.time < lowest) {
                    lowest = r.time;
                    st[k] = r.state;
                }
            }
        }
        for (int k = 0; k < nb; k++) {
            Node *n = flat.parents[k].first;
            if (below[k] > 0 and !n->is_sample) {
                n->set_site(x, 1 << st[k]);
            }
            if (below[k] > 0 and n->is_sample) {
                leaf = k;
            }
        }
        if (!mutation_sites.count(x)) {
            exception_sites[get_index(x)].push_back({x, column_base.empty() ? st[leaf] : column_base[(size_t) x]});
        }
        mutations[x] = rows;
    }
    while (model.continuous_gamma() and rates >> pos >> time) {
        site_rate[pos] = time;
    }
    return true;
}
