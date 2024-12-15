// LAD_1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>

using namespace std;

// Структура для точки данных
struct Point {
    vector<double> x; // Независимые переменные
    double y;         // Зависимая переменная
};

// RSS(b)= ∑∣Yi−A−B⋅Xi∣


double rss(int N, double* Y, double* X, double A, double B) {
    double sum_rss = 0.0;

    for (int i = 0; i < N; ++i) {
        double Y_predicted = A + B * X[i];
        double deviation = std::abs(Y[i] - Y_predicted);
        sum_rss += deviation;
    }

    return sum_rss;
}




void printPoints(const vector<Point>& points) {
    cout << "Сгенерированные точки:" << endl;
    for (size_t i = 0; i < points.size(); ++i) {
        cout << "Точка " << i + 1 << ": ";
        cout << "x = [";
        for (size_t j = 0; j < points[i].x.size(); ++j) {
            cout << points[i].x[j];
            if (j < points[i].x.size() - 1) cout << ", ";
        }
        cout << "] ";
        cout << "y = " << points[i].y << endl;
    }
}

// Функция вычисления целевой функции
double computeObjective(const vector<Point>& points, const vector<double>& coefficients) {
    double sum = 0.0;
    for (const auto& p : points) {
        double predicted = coefficients[0]; // b_0
        for (size_t j = 0; j < p.x.size(); ++j) {
            predicted += coefficients[j + 1] * p.x[j]; // b_1, b_2, ..., b_q
        }
        sum += abs(p.y - predicted); // Сумма абсолютных отклонений
    }
    return sum;
}

// Нахождение медианы
double findWeightedMedian(vector<pair<double, double>>& values) {
    // Сортировка по значениям
    sort(values.begin(), values.end());
    double totalWeight = 0.0;
    for (const auto& v : values) {
        totalWeight += v.second;
    }

    double cumulativeWeight = 0.0;
    for (const auto& v : values) {
        cumulativeWeight += v.second;
        if (cumulativeWeight >= totalWeight / 2.0) {
            return v.first;
        }
    }
    return values.back().first; // На случай численной ошибки
}

// Алгоритм Весоловского для многомерного случая
vector<double> wesolowskyAlgorithm(const vector<Point>& points, size_t dimensions) {
    vector<double> coefficients(dimensions + 1, 0.0); // Инициализация коэффициентов b_0, b_1, ..., b_q
    double minObjective = computeObjective(points, coefficients);
    bool improved = true;

    while (improved) {
        improved = false;

        for (size_t d = 0; d <= dimensions; ++d) { // Итерация по коэффициентам
            vector<pair<double, double>> weightedValues;

            for (const auto& p : points) {
                double residual = p.y - coefficients[0]; // Начальный остаток для b_0
                for (size_t j = 0; j < dimensions; ++j) {
                    if (j != d - 1) {
                        residual -= coefficients[j + 1] * p.x[j];
                    }
                }
                double weight = (d == 0) ? 1.0 : abs(p.x[d - 1]); // Веса для медианы
                double value = (d == 0) ? residual : residual / p.x[d - 1];
                weightedValues.emplace_back(value, weight);
            }

            double newCoefficient = findWeightedMedian(weightedValues);
            double oldCoefficient = coefficients[d];
            coefficients[d] = newCoefficient;

            double newObjective = computeObjective(points, coefficients);
            if (newObjective < minObjective) {
                minObjective = newObjective;
                improved = true;
            }
            else {
                coefficients[d] = oldCoefficient; // Возврат к старому значению
            }
        }
    }

    return coefficients;
}



vector<Point> generatePoints(size_t numPoints = 10, size_t numVariables= 10, double min = 0.0, double max = 10.0 ) {
    vector<Point> points;

    // Инициализация генератора случайных чисел
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min, max);

    for (size_t i = 0; i < numPoints; ++i) {
        Point point;
        point.x.resize(numVariables);

        // Генерация независимых переменных
        for (size_t j = 0; j < numVariables; ++j) {
            point.x[j] = dist(gen);
        }

        // Генерация зависимой переменной (например, линейная зависимость + шум)
        double noise = dist(gen) * 0.1; // Малый шум
        point.y = 0.0;
        for (size_t j = 0; j < numVariables; ++j) {
            point.y += point.x[j] * (j + 1); // Пример линейной зависимости
        }
        point.y += noise; // Добавление шума

        points.push_back(point);
    }

    return points;
}

int main() {
    setlocale(LC_ALL, "Russian");
    // Пример данных
    vector<Point> points;

    int numPoints, numVariables;
    cout << "введите количество измерений:" << endl;
    cin>> numPoints;
    cout << "введите количество переменных:" << endl;
    cin >> numVariables;


    points = generatePoints(numPoints, numVariables);
    printPoints(points);


   /* vector<Point> points = {
       {{1.0, 2.0}, 3.0},
       {{2.0, 1.0}, 4.0},
       {{3.0, 3.0}, 6.0},
       {{4.0, 5.0}, 8.0},
       {{5.0, 4.0}, 9.0}
    };*/



    size_t dimensions = numVariables; // Количество независимых переменных
    auto coefficients = wesolowskyAlgorithm(points, dimensions);

    // Вывод результатов
    cout << "Результаты регрессии:" << endl;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        cout << "b" << i << " = " << coefficients[i] << endl;
    }

    return 0;
}
