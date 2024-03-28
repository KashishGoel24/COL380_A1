#include<bits/stdc++.h>
#include <random>
#include<cmath>
#include<chrono>
using namespace std;
using namespace std::chrono;

double **A,**L,**U;
int* p;
int n;

void LUdecomp()
{
    for(int k=0;k<n;k++)
    {
        double pvt=0;
        int kmax=k;
        for(int i=k;i<n;i++)
        {
            if(fabs(A[i][k])>pvt)
            {
                pvt=fabs(A[i][k]);
                kmax=i;
            }
        }
        if(pvt==0)
        {
            cout<<"Singular matrix\n";
            return;
        }
        if(kmax!=k)
        {
            swap(A[k],A[kmax]);
            swap(p[k],p[kmax]);

            for(int i=0;i<k;i++)
            {
                swap(L[k][i],L[kmax][i]);
            }
        }
        U[k][k]=A[k][k];

        for(int i=k+1;i<n;i++)
        {
            L[i][k]=A[i][k]/U[k][k];
            U[k][i]=A[k][i];
        }

        for(int i=k+1;i<n;i++)
        {
            for(int j=k+1;j<n;j++)
            {
                A[i][j]=A[i][j]-L[i][k]*U[k][j];
            }
        }
    }
}

double checker(double** A)
{
    double **LU=new double*[n];
    for(int i=0;i<n;i++)
    {
        LU[i]=new double[n];
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            for(int k=0;k<n;k++)
            {
                LU[i][j]+=L[i][k]*U[k][j];
            }
        }
    }

    // cout<<"LU\n";
    // for(int i=0;i<n;i++){
    //     for(int j=0;j<n;j++){
    //         cout<<LU[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<"PA\n";
    // for(int i=0;i<n;i++){
    //     for(int j=0;j<n;j++){
    //         cout<<A[p[i]][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    double sum=0;
    for(int i=0;i<n;i++)
    {
        double tmp=0;
        for(int j=0;j<n;j++)
        {
            LU[j][i] = LU[j][i] - A[p[j]][i];
            tmp += LU[j][i]*LU[j][i];
        }
        sum+=sqrt(tmp);
    }
    return sum;
}

int main()
{
    cin>>n;
    A=new double*[n];
    L=new double*[n];
    U=new double*[n];
    p=new int[n];
    for(int i=0;i<n;i++)
    {
        A[i]=new double[n];
        L[i]=new double[n];
        U[i]=new double[n];
    }

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0);
    // cout<<"Matrix A:\n";
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = dis(gen);
            U[i][j] = 0;
            L[i][j] = 0;
            // cout<<A[i][j]<<" ";
        }
        // cout<<endl;
    }
    // cout<<endl;
    cout<<"Done Initialisation\n";
    double **A1=new double*[n];
    for(int i=0;i<n;i++)
    {
        A1[i]=new double[n];
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            A1[i][j]=A[i][j];
        }
    }
    for(int i=0;i<n;i++)
    {
        p[i]=i;
        L[i][i]=1.0;
    }
    auto start_time = high_resolution_clock::now(); 
    LUdecomp();
    auto end_time = high_resolution_clock::now();
    auto duration = (float) duration_cast<milliseconds>(end_time - start_time).count();
    cout << "Time taken for LU decomposition serial: " << duration/1000 << " seconds" << endl;
    cout<<checker(A1)<<endl;
    return 0;
}
