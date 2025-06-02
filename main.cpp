#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <cfloat> // Для DBL_MAX, DBL_MIN

using namespace std;

// Структура для возврата результата и погрешности
struct ArctanResult {
    double sum;
    double error;
};

double ulp(double value) {
    if (isnan(value))
        return numeric_limits<double>::quiet_NaN();
    if (isinf(value))
        return numeric_limits<double>::infinity();
    if (abs(value) == DBL_MAX)
        return pow(2.0, 971.0);
    if (abs(value) == 0.0)
        return numeric_limits<double>::denorm_min();
    if (abs(value) < DBL_MIN)
        return numeric_limits<double>::denorm_min();

    unsigned long long value_as_int;
    memcpy(&value_as_int, &value, sizeof(double));
    unsigned long long next_value_as_int = value_as_int + 1U;
    double next_value;
    memcpy(&next_value, &next_value_as_int, sizeof(double));
    return abs(next_value - value);
}

const double U = numeric_limits<double>::epsilon() / 2.0;

double gamma(int n) {
    return n * U / (1 - n * U);
}

ArctanResult arctan_and_errors(double x, int nMax) {
    double* terms = new double[nMax];
    double sum = 0.0;

    double x_power = x; // x^(2k+1), начнём с k=0 -> x^1 = x

    for (int k = 0; k < nMax; ++k) {
        int denominator = 2 * k + 1;
        double sign = (k % 2 == 0) ? 1.0 : -1.0;

        double term = sign * x_power / denominator;
        terms[k] = term;
        sum += term;

        // Подготовка x_power для следующего члена: умножаем на x^2
        x_power *= x * x;
    }

    //ошибка представления
    double error_presentation = sum / (1.0 + U);

    // Оценка погрешности обрыва ряда (остаточного члена ряда) - ok тех чл ряда которые не играют роли откидываются не влазят в дабл 
    double first_error = 1 / ((nMax + 1.0) * (nMax + 1));
    first_error += ulp(first_error);


    // Оценка погрешности округления каждого члена ряда - ok 
    double second_error = 0.0;
    for (int i = 0; i < nMax; ++i) {
        second_error += abs(terms[i]);
    }

    // Погрешность суммирования с ростом оценки - в main
    double third_error = gamma(nMax - 1) * second_error;
    second_error *= U;

    //total error
    double total_error = first_error + second_error + third_error + error_presentation;

    delete[] terms;
    return { sum, total_error };
}

int main() {
    double arctg5 = 1.0 / 5.0;
    double arctg239 = 1.0 / 239.0;

    ArctanResult atg1_5 = arctan_and_errors(arctg5, 11);
    ArctanResult atg1_239 = arctan_and_errors(arctg239, 4);

    double result1 = 16 * atg1_5.sum;
    double result2 = 4 * atg1_239.sum;

    double prodError = atg1_5.sum * U + atg1_239.sum * U; //ошибка при умножении

    double calculated_pi = 4.0 * (4.0 * atg1_5.sum - atg1_239.sum);

    double divError = abs(calculated_pi) * U; //ошибка деления

    double absolute = atg1_5.error + atg1_239.error + prodError + divError;

    cout << setprecision(16) << calculated_pi << endl;
    cout << "absolute error: " << fixed << setprecision(16) << absolute << endl;
    cout << "3.141592653589793" << endl;

    return 0;
}