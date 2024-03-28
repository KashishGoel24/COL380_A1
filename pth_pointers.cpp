#include<bits/stdc++.h>
#include <random>
#include<cmath>
#include<chrono>
#include<pthread.h>

using namespace std;
using namespace std::chrono;

double **A,**L,**U;
int* p;
int n,thread_count;
long int seed;

void print_matrix(double** A)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}  

void* matrixInitialisation(void* rank)
{
    int my_rank = *(int*)rank ;
    int start = my_rank*(n/thread_count);
    int end = (my_rank + 1)*(n/thread_count);
    srand48(time(nullptr) + my_rank);
    if(my_rank==thread_count-1)
    {
        end=n;
    }

    for (int row = start ; row < end ; row++)
    {
        for (int col = 0 ; col < n ; col++)
        {
            A[row][col] = drand48();
            L[row][col] = 0;
            U[row][col] = 0;
        }
    }
    return NULL;
}

void* update(void* rank)
{
    pair<int,int> rank1=*(pair<int,int>*)rank;
    int k=rank1.second;
    int start = rank1.first * ((n-k-1)/ thread_count);
    int end = (rank1.first + 1) * ((n-k-1)/ thread_count);
    if(rank1.first==thread_count-1)
    {
        end=n-k-1;
    }
    for(int i=start;i<end;i++)
    {
        for(int j=k+1;j<n;j++)
        {
            A[k+i+1][j]=A[k+i+1][j]-L[i+k+1][k]*U[k][j];
        }
    }
    return NULL;   
}

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
        double ukk=U[k][k];
        for(int i=k+1;i<n;i++)
        {
            L[i][k]=A[i][k]/ukk;
            U[k][i]=A[k][i];
        }

        pthread_t threads[thread_count];
        for(int i=0;i<thread_count;i++)
        {
            pthread_create(&threads[i],NULL,update,new pair<int,int>(i,k));
        }

        for(int i=0;i<thread_count;i++)
        {
            pthread_join(threads[i],NULL);
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
    // print_matrix(LU);
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
    cin>>thread_count;
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

    pthread_t threads[thread_count];
    int thread_ranks[thread_count];
    for (int i = 0; i < thread_count; i++) 
    {
        thread_ranks[i] = i;
        pthread_create(&threads[i], NULL, matrixInitialisation, (void*)&thread_ranks[i]);
    }
    for(int i=0;i<thread_count;i++)
    {
        pthread_join(threads[i],NULL);
    }


    for (int i = 0; i < n; i++)
    {
        p[i] = i; 
        L[i][i]=1.0;
    }
    // print_matrix(A);
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

    auto start_time = high_resolution_clock::now();
    LUdecomp();
    auto end_time = high_resolution_clock::now();
    auto duration = (float) duration_cast<milliseconds>(end_time - start_time).count();
    cout << "Time taken for LU decomposition(PTHREAD): " << duration/1000 << " seconds" << endl;

    //printing U and L
    // print_matrix(U);
    // cout<<endl;
    // print_matrix(L);

    cout<<"Error:"<<checker(A1)<<endl;

    return 0;
}