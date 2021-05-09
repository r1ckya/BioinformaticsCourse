#include <bits/stdc++.h>

using namespace std;

size_t overlap(const string& s1, const string& s2) {
  size_t res = min(s1.size(), s2.size());
  while (res > 0) {
    if (s2.compare(0, res, s1, s1.size() - res, res) == 0) {
      return res;
    }
    --res;
  }
  return 0;
}

void merge(string& s1, const string& s2, size_t overlap) {
  copy(s2.begin() + overlap, s2.end(), back_inserter(s1));
}

typedef tuple<size_t, size_t, size_t> vertex;
typedef set<vertex, greater<vertex>> vertex_set;

void remove(const vector<string>& reads, vertex_set& q, size_t i, size_t j) {
  if (i == j) return;
  if (reads[i].empty()) return;
  if (reads[j].empty()) return;
  q.erase({overlap(reads[i], reads[j]), i, j});
}

void add(const vector<string>& reads, vertex_set& q, size_t i, size_t j) {
  if (i == j) return;
  if (reads[i].empty()) return;
  if (reads[j].empty()) return;
  q.emplace(overlap(reads[i], reads[j]), i, j);
}

signed main() {
  cin.tie(0), cout.tie(0);
  ios_base::sync_with_stdio(0);

  vector<string> reads(istream_iterator<string>(cin), {});
  reads.resize(unique(reads.begin(), reads.end()) - reads.begin());

  vertex_set q;
  for (size_t i = 0; i < reads.size(); ++i) {
    for (size_t j = 0; j < reads.size(); ++j) {
      add(reads, q, i, j);
    }
  }

  size_t last = 0;
  while (!q.empty()) {
    auto [w, i, j] = *q.begin();
    q.erase(q.begin());
    for (size_t k = 0; k < reads.size(); ++k) {
      remove(reads, q, j, k);
      remove(reads, q, k, j);
      remove(reads, q, k, i);
      remove(reads, q, i, k);
    }
    merge(reads[i], reads[j], w);
    reads[j].clear();
    for (size_t k = 0; k < reads.size(); ++k) {
      add(reads, q, k, i);
      add(reads, q, i, k);
    }
    last = i;
  }
  cout << reads[last] << endl;
  return 0;
}
