#pragma once
#include <vector>
using namespace std;
void gauss(vector<vector<double>>& matrix, int N, int M);
void swapLines(int line1, int line2, vector<vector<double>>& a);
void getAnswerGauss(vector<vector<double>>& matrix, int N, int M, vector<double>& B);

//Matrisin ust ucgen haline getirilmesi.
void gauss(vector<vector<double>>& matrix, int N, int M)
{
    for (int k = 0; k < N - 1; k++)
    {
        double d = 0;
        if (matrix[k][k] == 0)
        {
            for (int i = k + 1; i < N; i++)
            {
                if (matrix[i][k] != 0)
                {
                    swapLines(k, i, matrix);
                    break;
                }
            }
        }
        for (int i = k + 1; i < N; i++)
        {
            d = matrix[i][k] / matrix[k][k];

            for (int j = k; j < M; j++)
            {
                matrix[i][j] = matrix[i][j] - matrix[k][j] * d;
            }
        }
    }
    for (int i = 0; i < N; i++) {
        double diagElement = matrix[i][i];
        for (int j = i; j < M; j++) {
            matrix[i][j] /= diagElement;
        }
    }
}
//Tersten gidilip sonuclarin hesaplanmasi.
void getAnswerGauss(vector<vector<double>>& matrix, int N, int M, vector<double>& ans)
{
    gauss(matrix, N, M);
    ans[N - 1] = matrix[N - 1][M - 1];
    for (int i = N - 2; i > -1; i--) {
        double sum = 0;
        for (int k = i + 1; k < N; k++) {
            sum += matrix[i][k] * ans[k];
        }
        ans[i] = matrix[i][N] - sum;
    }
}

void swapLines(int line1, int line2, vector<vector<double>>& a)
{
    vector<double> tmp = a[line1];
    a[line1] = a[line2];
    a[line2] = tmp;
}
