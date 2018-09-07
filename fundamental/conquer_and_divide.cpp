#include <limits>
#include <iostream>
#include <vector>
#include <tuple>

using namespace std;

template<typename T>
static constexpr T neginf()
{
	return numeric_limits<T>::min();
}

template<typename T>
tuple<int, int, T> find_max_crossing_subarray(T A[], int low, int mid, int high)
{
	int max_left = mid;
	T left_sum = neginf<T>();
	T sum = T();
	for(int i = mid; i >= low; i--){
		sum += A[i];
		if(sum > left_sum){
			left_sum = sum;
			max_left = i; 
		}
	}
	int max_right = mid+1;
	T right_sum = neginf<T>();
	sum = T();
	for(int j = mid+1; j <= high; j++){
		sum += A[j];
		if(sum > right_sum){
			right_sum = sum;
			max_right = j;
		}
	}

	return make_tuple(max_left, max_right, left_sum+right_sum);
}

template<typename T>
tuple<int, int, T> find_maximun_subarray(T A[], int low, int high)
{
	if(low == high){
		return make_tuple(low, high, A[low]);
	}else{
		int mid = (low+high)/2;
		int left_low, left_high;
		T left_sum;
		tie(left_low, left_high, left_sum) = 
			find_maximun_subarray(A, low, mid);
		int right_low, right_high;
		T right_sum;
		tie(right_low, right_high, right_sum) = 
			find_maximun_subarray(A, mid+1, high);
		int cross_low, cross_high;
		T cross_sum;
		tie(cross_low, cross_high, cross_sum) = 
			find_max_crossing_subarray(A, low, mid, high);
		if(left_sum >= right_sum && left_sum >= cross_sum)
			return make_tuple(left_low, left_high, left_sum);
		else if(right_sum >= left_sum && right_sum >= cross_sum)
			return make_tuple(right_low, right_high, right_sum);
		else
			return make_tuple(cross_low, cross_high, cross_sum);
	}
}

int main()
{
	vector<int> A = {13, -3, -25, 20, -3, -16, -23, 18, 20, -7, 12, -5, -22, 15, -4, 7};
	int low, high, sum;
	tie(low, high, sum) = find_maximun_subarray(A.data(), 0, A.size()-1);
	cout << low << " " << high << " " << sum << endl;
}