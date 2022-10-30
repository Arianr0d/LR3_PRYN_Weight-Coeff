#include <iostream>;
#include <sstream>;
#include <fstream>;
#include <cmath>;
#include <string>;
#include <vector>;

using namespace std;

/*
	задача теории принятия решений с весовыми коэффициентами 
	для нечетких отношений предпочтения
*/

// функция записи результата в файл
void output(pair<int, float> &result, vector<float>& mu) {

	ofstream out;
	out.open("output.txt");
	if (out.is_open()) {
		out << "Ответ:\n";
		for (int i = 0; i < mu.size(); i++) {
			out << "phi_н.д.(x" << i + 1 <<
				") = " << mu[i] << "; " << endl;
		}
		out << "\nоптимальная альтернатива - x" << result.first << ", phi_н.д.(x" << result.first <<
			") = " << result.second << endl;
	}
}

// фнкция поиска максимума и степени недоминируемости матрицы Q_s
void funDegreeDominance(int r_size, vector<vector<float>> &Q_s, vector<float> &Q_mu) {

	for (int j = 0; j < r_size; j++) {
		float max_value = Q_s[0][j];
		for (int i = 0; i < r_size; i++) {
			if (Q_s[i][j] > max_value) max_value = Q_s[i][j];
		}
		Q_mu.push_back(1 - max_value);
	}
}

// функция нечеткого отношения строгого предпочтения
void funFuzzyRelation(int r_size, vector<vector<float>> &Q, vector<vector<float>> &Q_s) {

	for (int i = 0; i < r_size; i++) {
		for (int j = 0; j < r_size; j++) {
			if (Q[i][j] - Q[j][i] < 0) Q_s[i][j] = 0;
			else Q_s[i][j] = Q[i][j] - Q[j][i];
		}
	}
}

// алгоритм поиска решения в задаче с весовыми коэффициентами
// для нечетких отношений предпочтения
void algorithm(int r_size, 
			   vector<float>& lambda, 
			   vector<vector<vector<float>>>& R, 
			   pair<int, float> &result, 
			   vector<float> &mu) {
	
	vector<vector<float>> Q1(r_size, vector<float>(r_size));
	vector<vector<float>> Q1_s(r_size, vector<float>(r_size));
	vector<float> Q1_mu;

	// поиск матрицы Q1
	for (int i = 0; i < r_size; i++) {
		for (int j = 0; j < r_size; j++) {
			float min_value = R[0][i][j];
			for (int k = 0; k < R.size(); k++) {
				if (R[k][i][j] < min_value) min_value = R[k][i][j];
			}
			Q1[i][j] = min_value;
		}
	}
	// поиск матрицы Q1s
	funFuzzyRelation(r_size, Q1, ref(Q1_s));
	// поиск максимума и степени недоминируемости матрицы Q1_s
	funDegreeDominance(r_size, Q1_s, ref(Q1_mu));

	vector<vector<float>> Q2(r_size, vector<float>(r_size));
	vector<vector<float>> Q2_s(r_size, vector<float>(r_size));
	vector<float> Q2_mu;
	// поиск матрицы Q2
	for (int i = 0; i < r_size; i++) {
		for (int j = 0; j < r_size; j++) {
			for (int k = 0; k < R.size(); k++) {
				Q2[i][j] += lambda[k] * R[k][i][j];
			}
		}
	}
	// поиск матрицы Q2s
	funFuzzyRelation(r_size, Q2, ref(Q2_s));
	// поиск максимума и степени недоминируемости матрицы Q2_s
	funDegreeDominance(r_size, Q2_s, ref(Q2_mu));

	float min_value = 0.0;
	// поиск общих степеней недоминируемости альтернатив
	for (int i = 0; i < Q1_mu.size(); i++) {
		if (Q1_mu[i] < Q2_mu[i]) min_value = Q1_mu[i];
		else min_value = Q2_mu[i];
		mu.push_back(min_value);

		if (i == 0) result = make_pair(1, min_value);
		if (min_value > result.second) result = make_pair(i + 1, min_value);
	}
}

// функция считывания матриц из файла
void importFromFile(int &lambda_count, int &r_size, int &r_count, vector<float> &lambda, vector<vector<vector<float>>> &R) {

	string str;
	vector<vector<float>> r_matrix(r_size, vector<float> (r_size));
	ifstream in("input.txt");
	if (in.is_open()) {
		int i = 0,
			k = -1;

		getline(in, str);
		r_size = stoi(str);  // считывание размера матрицы R
		getline(in, str);
		r_count = stoi(str);  // считывание количества матриц R
		getline(in, str);
		lambda_count = stoi(str);  // считывание количества весовых коэффициентов

		// считывание весовых коэффициентов
		getline(in, str);
		istringstream ist(str);
		string row;
		while (ist >> row) {
			lambda.push_back(stod(row));
		}

		// считывание нечетких отношений предпочтения
		while (getline(in, str)) {

			if (str.size() == 0) continue;

			string row;
			istringstream ist(str);
			int j = 0;
			if (i % r_size == 0) {
				k++;
				if(k >= 1) R.push_back(r_matrix);
			}

			while (ist >> row) {
				if (k >= r_count) break;
				r_matrix[i % r_size][j] = stod(row);
				j++;
			}
			i++;
		}
		R.push_back(r_matrix);
	}
	in.close();
}

int main() {
	setlocale(LC_ALL, "ru");

	int r_size = 6,
		r_count = 7,
		lambda_count = 7;
	vector<float> lambda;
	vector<vector<vector<float>>> R;
	vector<float> mu; // вектор степеней недоминируемости альтернатив

	cout << "Поиск решения задачи с весовыми коэффициентами для нечетких отношений предпочтения" << endl << endl;
	cout << "Считывание данных из файла!" << endl;

	// считывание значений из файла
	importFromFile(lambda_count, r_size, r_count, std::ref(lambda), std::ref(R));

	pair<int, float> result;
	// реализация алгоритма
	algorithm(r_size, lambda, R, ref(result), ref(mu));

	// запись результата в файл
	output(result, mu);

	cout << "Результаты успешно сохранены в файл!" << endl;
	system("pause");
}