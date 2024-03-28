#include<bits/stdc++.h>
#include <random>
#include<cmath>
#include<chrono>
#include<omp.h>
using namespace std::chrono;
using namespace std;

double **A,**L,**U;
int* p;
int n,num_threads;

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

            #pragma omp parallel for num_threads(num_threads)
            for(int i=0;i<k;i++)
            {
                swap(L[k][i],L[kmax][i]);
            }
        }

        U[k][k]=A[k][k];
        double ukk=U[k][k];

        #pragma omp parallel for num_threads(num_threads)
        for(int i=k+1;i<n;i++)
        {
            L[i][k]=A[i][k]/ukk;
            U[k][i]=A[k][i];
        }

        #pragma omp parallel for num_threads(num_threads)
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

    #pragma omp parallel for num_threads(num_threads)
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
    // for(int i=0;i<n;i++)
    // {
    //     for(int j=0;j<n;j++)
    //     {
    //         cout<<LU[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // cout<<"PA\n";
    // for(int i=0;i<n;i++)
    // {
    //     for(int j=0;j<n;j++)
    //     {
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
            LU[j][i]=LU[j][i]-A[p[j]][i];
            tmp+=LU[j][i]*LU[j][i];
        }
        sum+=sqrt(tmp);
    }
    return sum;
}

int main()
{
    cin>>n;
    cin>>num_threads;
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

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        srand48(time(nullptr) + tid);

        #pragma omp for
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = drand48();
                L[i][j] = 0;
                U[i][j] = 0;
            }
        }
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
    double start_time = omp_get_wtime();
    LUdecomp();    
    double end_time = omp_get_wtime();
    cout << "Time taken for LU decomposition: OMP " << (end_time - start_time) << " seconds\n";
    
    cout << "Error: " << checker(A1) << endl;
    return 0;
}