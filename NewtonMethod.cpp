#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

void gauss(vector<vector<double>>& matrix, int N, int M);
void swapLines(int line1, int line2, vector<vector<double>>& a);
void getAnswerGauss(vector<vector<double>>& matrix, int N, int M, vector<double>& B);

double getD1(const vector<double>& x);
double getD2(const vector<double>& x, const vector<double>& xKplus1);
double FirstFunction(const double& x01, const double& x02);
double SecondFunction(const double& x01, const double& x02);
double jacobianFirstFuncX1(const double& x01, const double& x02);
double jacobianFirstFuncX2(const double& x01, const double& x02);
double jacobianSecondFuncX1(const double& x01, const double& x02);
double jacobianSecondFuncX2(const double& x01, const double& x02);
void jacobian(vector<vector<double>>& matrix, const vector<double>& x, const double& M);

double getD1(const vector<double>& x)
{
    double f1 = FirstFunction(x[0], x[1]);
    double f2 = SecondFunction(x[0], x[1]);
    return (abs(f1) > abs(f2) ? abs(f1) : abs(f2));
}

double getD2(const vector<double>& x, const vector<double>& xKplus1)
{
    double max = 0;
    for (int i = 0; i < x.size(); i++)
    {
        if (abs(xKplus1[i]) < 1)
        {
            if (abs(xKplus1[i] - x[i]) > max)
                max = abs(xKplus1[i] - x[i]);
        }
        else
        {
            if (abs((xKplus1[i] - x[i]) / xKplus1[i]) > max)
                max = abs((xKplus1[i] - x[i]) / xKplus1[i]);
        }
    }
    return max;
}

double FirstFunction(const double& x01, const double& x02)
{
    return x01 * x01 * x02 * x02 - 3 * x01 * x01 - 6 * x02 * x02 * x02 + 8;
}

double SecondFunction(const double& x01, const double& x02)
{
    return x01 * x01 * x01 * x01 - 9 * x02 + 2;
}

double jacobianFirstFuncX1(const double& x01, const double& x02)
{
    return 2 * x01 * x02 * x02 - 6 * x01;
}

double jacobianFirstFuncX2(const double& x01, const double& x02)
{
    return 2 * x02 * x01 * x01 - 18 * x02 * x02;
}

double jacobianSecondFuncX1(const double& x01, const double& x02)
{
    return 4 * x01 * x01 * x01;
}

double jacobianSecondFuncX2(const double& x01, const double& x02)
{
    return -9;
}


void jacobian(vector<vector<double>>& matrix, const vector<double>& x, const double& M)
{
    matrix[0][0] = (FirstFunction(x[0] + M * x[0], x[1]) - FirstFunction(x[0], x[1])) / M * x[0];
    matrix[0][1] = (FirstFunction(x[0], x[1] + M * x[1]) - FirstFunction(x[0], x[1])) / M * x[1];
    matrix[1][0] = (SecondFunction(x[0] + M * x[0], x[1]) - SecondFunction(x[0], x[1])) / M * x[0];
    matrix[1][1] = (SecondFunction(x[0], x[1] + M * x[1]) - SecondFunction(x[0], x[1])) / M * x[1];
}

void compareNewtonMethod(const double &mParameter)
{
    const double eps = 1e-9;
    double d1 = 1, d2 = 1;
    vector<double> x = { -1.5, 1.5 };
    vector<double> xkPlus1(x);
    int k = 1;
    int NIT = 100;
    int N = 2;
    int M = 3;
    vector<vector<double>> matrix(N, vector<double>(M));
    vector<double> ans(N);

    while ((d1 > eps || d2 > eps))
    {
        if (k >= NIT)
        {
            cout << "IER=2" << endl;
            break;
        }

        else
        {
            jacobian(matrix, x, mParameter);

            matrix[0][2] = -FirstFunction(x[0], x[1]);
            matrix[1][2] = -SecondFunction(x[0], x[1]);

            getAnswerGauss(matrix, N, M, ans);
            xkPlus1[0] += ans[0];
            xkPlus1[1] += ans[1];

            d1 = getD1(x);
            d2 = getD2(x, xkPlus1);
            cout << "Iteration:" << k << " d1: " << d1 << " d2: " << d2 << endl;
            k++;
            x = xkPlus1;
        }
    }
    cout << "Solution: x = " << setprecision(20) << x[0] << setw(8) << " y = " << x[1] << endl << "Number of iterations: " << k << endl;
}

int main()
{
    vector<double> mParameter = {0.01, 0.05, 0.1 };
    for (const double& M : mParameter)
    {
        cout << " Parameter M = " << M << endl;
        compareNewtonMethod(M);
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
    }

    return 0;
}