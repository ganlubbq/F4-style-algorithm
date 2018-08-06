#pragma once

#include "Poly.h"
#include <iostream>
//windows
#ifdef _MSC_VER
#include <intrin.h>
//linux
#else
#include <x86intrin.h>
#endif

using namespace std;

// 256bitにucharが32個入る
static const size_t single_size = 32;
static const size_t ALIGNMENT = 32;

//macro
#define TYPE_AVX __m256i
// u : unsigned char(8bit) を256bitに32個敷き詰める
#define SET_VALUE_AVX(u) _mm256_set_epi8(u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u,u)
// 演算子
#define XOR_AVX(a,b) _mm256_xor_si256(a,b)
#define OR_AVX(a,b) _mm256_or_si256(a,b)
#define AND_AVX(a,b) _mm256_and_si256(a,b)
#define RIGHT_AVX(a,cnt) _mm256_srli_epi16(a,cnt)
#define LEFT_AVX(a,cnt) _mm256_slli_epi16(a,cnt)
#define ADD_AVX(a,b) _mm256_add_epi8(a,b)
#define SUB_AVX(a,b) _mm256_sub_epi8(a,b)
// ロード・セーブ
#define LOAD_AVX(arr) _mm256_load_si256((__m256i*)(&arr))
#define STORE_AVX(arr,v) _mm256_store_si256((__m256i*)(&arr),v)

/*
CONST_0xYY はYYは16進数での定数。
マクロ定義でも良かったが、一応static constで定義。
*/
static const unsigned char CONST_0x01 = 1;
static const unsigned char CONST_0x03 = 3;
static const unsigned char CONST_0x07 = 7;
static const unsigned char CONST_0x0f = 15;
static const unsigned char CONST_0x10 = 16;
static const unsigned char CONST_0x18 = 24;
static const unsigned char CONST_0x1c = 28;
static const unsigned char CONST_0x1e = 30;
static const unsigned char CONST_0x1f = 31;
static const unsigned char CONST_0x20 = 32;

/*
AVX演算に用いる定数を事前に用意しておく。
static constで定義することで、呼び出すだけで用いることができる。

構造は、AVX_CONST_0xYYは256bitの変数になっており
8bitのunsigned char YY が32個連続して格納されている.
*/
static const auto AVX_CONST_0x01 = SET_VALUE_AVX(CONST_0x01);
static const auto AVX_CONST_0x03 = SET_VALUE_AVX(CONST_0x03);
static const auto AVX_CONST_0x07 = SET_VALUE_AVX(CONST_0x07);
static const auto AVX_CONST_0x0f = SET_VALUE_AVX(CONST_0x0f);
static const auto AVX_CONST_0x10 = SET_VALUE_AVX(CONST_0x10);
static const auto AVX_CONST_0x18 = SET_VALUE_AVX(CONST_0x18);
static const auto AVX_CONST_0x1c = SET_VALUE_AVX(CONST_0x1c);
static const auto AVX_CONST_0x1e = SET_VALUE_AVX(CONST_0x1e);
static const auto AVX_CONST_0x1f = SET_VALUE_AVX(CONST_0x1f);
static const auto AVX_CONST_0x20 = SET_VALUE_AVX(CONST_0x20);

class Poly_GF31_simd : public Poly
{	
public:
	Poly_GF31_simd(vector<unsigned char> &coeff) : Poly(coeff) 
	{
		_Coeff_size = coeff.size();
		_Div_single_size = _Coeff_size / single_size;
	};

	//variable
	int _Coeff_size;//coeffのsize
	int _Div_single_size;//coeff_size をsingle_sizeで割った値　simdまわす数

	//function
	//unsigned char MOD31(unsigned char a);
	void MOD31(TYPE_AVX& vec);
	void MOD31_First(TYPE_AVX& vec);
	void mul0();
	void mul1();
	void mul2();
	void mul3();
	void mul4();
	void mul5();
	void mul6();
	void mul7();
	void mul8();
	void mul9();
	void mul10();
	void mul11();
	void mul12();
	void mul13();
	void mul14();
	void mul15();
	void mul16();
	void mul17();
	void mul18();
	void mul19();
	void mul20();
	void mul21();
	void mul22();
	void mul23();
	void mul24();
	void mul25();
	void mul26();
	void mul27();
	void mul28();
	void mul29();
	void mul30();
	
	//operator
	void operator+(Poly &poly);
	void operator*(unsigned char &n);
	void operator*(vector<unsigned char> &monomial_deg);
};