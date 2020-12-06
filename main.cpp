#include <iostream>
#include <vector>
#include <limits>
#include <math.h>

using namespace std;

const double EPS = 0.00001;
const double INF = std::numeric_limits<double>::infinity();

int gauss_jordan (vector < vector<double> > a, vector<double> & ans) {
	int n = (int) a.size();
	int m = (int) a[0].size() - 1;

	vector<int> where (m, -1);
	for (int col=0, row=0; col<m && row<n; ++col) {
		int sel = row;
		for (int i=row; i<n; ++i) {
            if (abs(a[i][col]) > abs(a[sel][col])) {
                sel = i;
            }
        }
		if (abs (a[sel][col]) < EPS) {
            continue;
        }
		for (int i=col; i<=m; ++i) {
            swap(a[sel][i], a[row][i]);
        }
		where[col] = row;

		for (int i=0; i<n; ++i)
			if (i != row) {
				double c = a[i][col] / a[row][col];
				for (int j=col; j<=m; ++j) {
                    a[i][j] -= a[row][j] * c;
                }
			}
		++row;
	}

	ans.assign (m, 0);
	for (int i=0; i<m; ++i) {
        if (where[i] != -1) {
            ans[i] = a[where[i]][m] / a[where[i]][i];
        }
    }
	for (int i=0; i<n; ++i) {
		double sum = 0;
		for (int j=0; j<m; ++j)
			sum += ans[j] * a[i][j];
		if (abs (sum - a[i][m]) > EPS) {
            return 0;
        }
	}

	for (int i=0; i<m; ++i) {
        if (where[i] == -1) {
            return INF;
        }
    }
	return 1;
}


int main() {
    vector <double> ans;
    int n, m;
    cout << "Введите количество строк и столбцов матрицы: " << endl;
    cin >> n >> m;
    vector < vector<double> > Matrix(n, vector<double>(m));
    cout << "Введите элементы расширенной матрицы A: " << endl;
    for (int i(0); i < n; ++i) {
        for (int j(0); j < m; ++j) {
            cin >> Matrix[i][j];
        }
    }
    int res = gauss_jordan(Matrix, ans);
    cout << "Количество решений системы: " << res << endl;
    if (res == 0 || res == 1) {
        cout << "Система имеет следующие решения: " << endl;
        for (int i(0); i < ans.size(); ++i) {
            cout << ans[i] << " ";
        }
    }
    return 0;
}