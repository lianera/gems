#include <iostream>
#include <limits>
#include <vector>

using std::vector;
using std::numeric_limits;
using std::cout;
using std::endl;

template<typename T>
static constexpr T inf()
{
	return numeric_limits<T>::max();
}

template<typename T>
static void swap(T& a, T& b)
{
	T t = a;
	a = b;
	b = t;
}

template<typename T>
void insertion_sort(T A[], int n)
{
	for(int i = 1; i < n; i++){
		T key = A[i];
		int j = i-1;
		while(j >= 0 && A[j] > key){
			A[j] = A[j+1];
			j--;			
		}
		A[j+1] = key;
	}
}

template<typename T>
void merge(T A[], int p, int q, int r)
{
	int n1 = q-p+1;
	int n2 = r-q;
	T* L = new T[n1+1];
	T* R = new T[n2+1];
	for(int i = 0; i < n1; i++)
		L[i] = A[p+i];
	for(int j = 0; j < n2; j++)
		R[j] = A[q+1+j];
	L[n1] = inf<T>();
	R[n2] = inf<T>();
	int i = 0, j = 0;
	for(int k = p; k <= r; k++){
		if(L[i] <= R[j]){
			A[k] = L[i];
			i++;
		}else{
			A[k] = R[j];
			j++;
		}
	}
	delete[] R;
	delete[] L;
}

template<typename T>
void merge_sort(T A[], int p, int r)
{
	if(p < r){
		int q = (p+r)/2;
		merge_sort(A, p, q);
		merge_sort(A, q+1, r);
		merge(A, p, q, r);
	}
}

template<typename T>
void bubble_sort(T A[], int n)
{
	for(int i = 0; i < n-1; i++){
		for(int j = n-1; j > i; j--)
			if(A[j] < A[j-1])
				swap(A[j], A[j-1]);
	}
}

template<typename T>
static void validate(vector<T> A)
{
	for(int i = 0; i < A.size()-1; i++){
		if(A[i] > A[i+1]){
			cout << "invalid" << endl;
			break;
		}
	}
}

int main()
{
	int N = 10000;
	vector<int> A(N);
	for(int& a:A)
		a = rand();
	auto B = A;
	insertion_sort(B.data(), N);
	validate(B);
	B = A;
	merge_sort(B.data(), 0, N-1);
	validate(B);
	B = A;
	bubble_sort(B.data(), N);
	validate(B);	
	cout << "done" << endl;
	return 0;
}