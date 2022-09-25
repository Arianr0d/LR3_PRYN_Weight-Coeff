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
void output(pair<int, float> &result) {

	ofstream out;
	out.open("output.txt");
	if (out.is_open()) {
		out << "Ответ: оптимальная альтернатива - x" << result.first << ", phi_н.д.(x" << result.first <<
			") = " << result.second << endl;
	}
}

// алгоритм поиска решения в задаче с весовыми коэффициентами
// для нечетких отношений предпочтения
void algorithm(int r_size, vector<float>& lambda, vector<vector<vector<float>>>& R, pair<int, float> &result) {
	
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
	for (int i = 0; i < r_size; i++) {
		for (int j = 0; j < r_size; j++) {
			if (Q1[i][j] - Q1[j][i] < 0) Q1_s[i][j] = 0;
			else Q1_s[i][j] = Q1[i][j] - Q1[j][i];
		}
	}
	// поиск максимума и степени недоминируемости матрицы Q1_s
	for (int j = 0; j < r_size; j++) {
		float max_value = Q1_s[0][j];
		for (int i = 0; i < r_size; i++) {
			if (Q1_s[i][j] > max_value) max_value = Q1_s[i][j];
		}
		Q1_mu.push_back(1 - max_value);
	}

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
	for (int i = 0; i < r_size; i++) {
		for (int j = 0; j < r_size; j++) {
			if (Q2[i][j] - Q2[j][i] < 0) Q2_s[i][j] = 0;
			else Q2_s[i][j] = Q2[i][j] - Q2[j][i];
		}
	}
	// поиск максимума и степени недоминируемости матрицы Q2_s
	for (int j = 0; j < r_size; j++) {
		float max_value = Q2_s[0][j];
		for (int i = 0; i < r_size; i++) {
			if (Q2_s[i][j] > max_value) max_value = Q2_s[i][j];
		}
		Q2_mu.push_back(1 - max_value);
	}

	float min_value = 0.0;
	// поиск общих степеней недоминируемости альтернатив
	for (int i = 0; i < Q1_mu.size(); i++) {
		if (Q1_mu[i] < Q2_mu[i]) min_value = Q1_mu[i];
		else min_value = Q2_mu[i];

		if (i == 0) result = make_pair(1, min_value);
		if (min_value > result.second) result = make_pair(i + 1, min_value);
	}
}

// функция считывания матриц из консоли
void inputMatrix(int lambda_count, int r_size, int r_count, vector<float>& lambda, vector<vector<vector<float>>>& R) {

	string row;
	int j = 0;

	cin.ignore(cin.rdbuf()->in_avail());
	cout << "Введите значения весовых коэффициентов: \n";
	
	getline(cin, row);
	istringstream ist(row);
	while (ist >> row) {
		lambda.push_back(stod(row));
		j++;
	}

	cout << "Укажите матрицы нечетких отношений предпочтения \n";
	for (int k = 0; k < r_count; k++) {
		cout << "\nМатрица R" + to_string(k + 1) << endl;
		vector<vector<float>> r(r_size, vector<float>(r_size));
		
		for (int i = 0; i < r_size; i++) {
			j = 0;
			getline(cin, row);
			istringstream ist(row);

			while (ist >> row) {
				r[i][j] = stod(row);
				j++;
			}
		}
		R.push_back(r);
	}
}

// функция считывания матриц из файла
void importFromFile(int lambda_count, int r_size, int r_count, vector<float> &lambda, vector<vector<vector<float>>> &R) {

	string str;
	vector<vector<float>> r_matrix(r_size, vector<float> (r_size));
	ifstream in("input.txt");
	if (in.is_open()) {
		int i = 0,
			k = -1;

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

	string ans;
	int r_size = 6,
		r_count = 7,
		lambda_count = 7;
	vector<float> lambda;
	vector<vector<vector<float>>> R;

	cout << "Поиск решения задачи с весовыми коэффициентами для нечетких отношений предпочтения" << endl << endl;
	while (true) {
		cout << "Использовать значения по умолчанию (2 вариант): (yes / no) ";
		cin >> ans;
		if (ans == "no" || ans == "yes") break;
	}

	if (ans == "no") {
		cout << "Введите размер матриц R: ";
		cin >> r_size;
		cout << "Введите число матриц R: ";
		cin >> r_count;
		cout << "Введите число весовых коэффициентов lambda: ";
		cin >> lambda_count;
		// считывание с консоли данных
		inputMatrix(lambda_count, r_size, r_count, std::ref(lambda), std::ref(R));
	}

	// считывание значений по умолчанию
	if (ans == "yes") importFromFile(lambda_count, r_size, r_count, std::ref(lambda), std::ref(R));

	pair<int, float> result;
	// реализация алгоритма
	algorithm(r_size, lambda, R, ref(result));

	// запись результата в файл
	output(result);

	cout << "Результаты записаны в файл!";

}