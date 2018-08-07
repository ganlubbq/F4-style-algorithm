#include "GF31.h"

/*
inline void GF31::MOD31(TYPE_AVX& vec) {
	vec = ADD_AVX(vec, AVX_CONST_0x01);
	vec = SUB_AVX(AND_AVX(vec, AVX_CONST_0x1f), RIGHT_AVX(XOR_AVX(AND_AVX(vec, AVX_CONST_0x20), AVX_CONST_0x20), 0x5));
}

inline void GF31::MOD31_First(TYPE_AVX& vec) {
	vec = ADD_AVX(AND_AVX(vec, AVX_CONST_0x1f), RIGHT_AVX(AND_AVX(vec, AVX_CONST_0x20), 5));
}

inline void GF31::mul0() {
	cout << "mul0 occur error" << endl;
	system("pause");
	return;
}

inline void GF31::mul1() {
	return;
}

// vec *= 2*vec2の計算
inline void GF31::mul2() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		//vec += ((vec2 & 0x10) >> 4) + ((vec2 & 0x0f) << 1);
		vx[i] = ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 0x4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1));
		MOD31(vx[i]);
	}

	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (2 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul3() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		//vec += 2*vec2 + vec2;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 0x4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)), vx[i]);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (3 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul4() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += ((vec2 & 0x18) >> 3) + ((vec2 & 0x7) << 2);
		vx[i] = ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 0x3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 0x2));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (4 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul5() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += 4*vec2 + vec2;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 0x3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 0x2)), vx[i]);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (5 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul6() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += 4*vec + 2*vec;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (6 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul7() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += 8*vec + vec^0x1f;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)), XOR_AVX(vx[i], AVX_CONST_0x1f));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (7 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul8() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += (vec2 & 0x1c) >> 2 + (vec2 & 0x3) << 3;
		vx[i] = ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (8 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul9() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += 8*vec2 + vec2;
		vx[i] = ADD_AVX(vx[i], ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (9 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul10() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += 8*vec2 + 2*vec2;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (10 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul11() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (int i = 0; i < _Div_single_size; ++i) {
		// vec += (4 * vec2 + 16 * vec2) ^ 0x1f;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4)));
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (11 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul12() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		//vec += (vec2 * 4 + vec2 * 8);
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (12 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul13() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += (vec2 * 2 + vec2 * 16) ^ 0x1f;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4)));
		MOD31_First(vx[i]);
		vx[i] = XOR_AVX(vx[i], AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (13 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul14() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += (vec2 * 17) ^ 0x1f;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4)), vx[i]);
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (14 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul15() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += (vec2 * 16) ^ 0x1f;
		vx[i] = XOR_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4)), AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (15 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul16() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += ((vec2 & 0x1e) >> 1) + ((vec2 & 0x1) << 4);
		vx[i] = ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (16 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul17() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += 16 * vec2 + vec2;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4)), vx[i]);
		//MOD31_First(temp);
		//vx[i] = ADD_AVX(temp, vx[i]);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (17 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul18() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec += 16 * vec + 2 * vec;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4)));
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (18 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul19() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (8 * vec + 4 * vec) ^ 0x1f;
		temp = LOAD_AVX(_Coeff[i]);
		//cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << _Div_single_size l;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)));
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (19 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul20() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = 16 * vec + 4 * vec;
		vx[i] = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1e), 1), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x01), 4)));;
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (20 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul21() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (8 * vec + 2 * vec) ^ 0x1f;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)));
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (21 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul22() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (9 * vec) ^ 0x1f;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)), vx[i]);
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (22 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul23() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (8 * vec) ^ 0x1f;
		vx[i] = XOR_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)), AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (23 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul24() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (8 * vec) ^ 0x1f + vec;
		// vec = 23 * vec + vec ;
		vx[i] = ADD_AVX(XOR_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x1c), 2), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x03), 3)), AVX_CONST_0x1f), vx[i]);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (24 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul25() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (2 * vec + 4 * vec) ^ 0x1f;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)), ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)));
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (25 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul26() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (5 * vec) ^ 0x1f;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)), vx[i]);
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (26 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul27() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (4 * vec) ^ 0x1f;
		vx[i] = XOR_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x18), 3), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x07), 2)), AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (27 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul28() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	TYPE_AVX temp;
	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (3 * vec) ^ 0x1f;
		temp = ADD_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)), vx[i]);
		MOD31_First(temp);
		vx[i] = XOR_AVX(temp, AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (28 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul29() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = (2 * vec) ^ 0x1f;
		vx[i] = XOR_AVX(ADD_AVX(RIGHT_AVX(AND_AVX(vx[i], AVX_CONST_0x10), 4), LEFT_AVX(AND_AVX(vx[i], AVX_CONST_0x0f), 1)), AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (29 * _Coeff[i]) % 31;
	}
}

inline void GF31::mul30() {

	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		// vec = vec ^ 0x1f;
		vx[i] = XOR_AVX(vx[i], AVX_CONST_0x1f);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i <_Coeff_size; ++i) {
		_Coeff[i] = (30 * _Coeff[i]) % 31;
	}
}

inline void GF31::operator+(GF31 poly)
{
	//AVX 専用の型にデータをロードする
	TYPE_AVX *vx = (TYPE_AVX *) &(_Coeff[0]);
	TYPE_AVX *vy = (TYPE_AVX *) &(poly._Coeff[0]);

	for (size_t i = 0; i < _Div_single_size; ++i) {
		//vec += vec2;
		vx[i] = ADD_AVX(vx[i], vy[i]);
		MOD31(vx[i]);
	}

	// SIMD計算で残った部分
	for (size_t i = _Div_single_size * single_size; i < _Coeff_size; ++i) {
		_Coeff[i] = (_Coeff[i] + poly._Coeff[i]) % 31;
	}
}

//coeff倍
inline void GF31::operator*(unsigned char &n)
{
	switch (n)
	{
	case 1:
		mul1();
		break;
	case 2:
		mul2();
		break;
	case 3:
		mul3();
		break;
	case 4:
		mul4();
		break;
	case 5:
		mul5();
		break;
	case 6:
		mul6();
		break;
	case 7:
		mul7();
		break;
	case 8:
		mul8();
		break;
	case 9:
		mul9();
		break;
	case 10:
		mul10();
		break;
	case 11:
		mul11();
		break;
	case 12:
		mul12();
		break;
	case 13:
		mul13();
		break;
	case 14:
		mul14();
		break;
	case 15:
		mul15();
		break;
	case 16:
		mul16();
		break;
	case 17:
		mul17();
		break;
	case 18:
		mul18();
		break;
	case 19:
		mul19();
		break;
	case 20:
		mul20();
		break;
	case 21:
		mul21();
		break;
	case 22:
		mul22();
		break;
	case 23:
		mul23();
		break;
	case 24:
		mul24();
		break;
	case 25:
		mul25();
		break;
	case 26:
		mul26();
		break;
	case 27:
		mul27();
		break;
	case 28:
		mul28();
		break;
	case 29:
		mul29();
		break;
	case 30:
		mul30();
		break;
	default:
		cout << "operator* error" << endl;
		system("pause");
	}
}

//monomialの係数は1が前提　Coeffのスライドをする LMDeg,LMdeg_index,Coeff_size,div_single_size,Max_degreeを更新
inline void GF31::operator*(vector<unsigned char> &monomial_deg)
{
	//もろもろ更新
	_LMdeg = _Degree.vec_add(_LMdeg, monomial_deg);
	_LMdeg_index = _Degree.degree_to_index(_LMdeg);

	//_Coeffのresizeと_Coeff_sizeとdiv_single_size更新
	while (_Coeff_size < _LMdeg_index + 1)
	{
		_Coeff_size << 1;
	}
	_Coeff.resize(_Coeff_size);
	_Div_single_size = _Coeff_size / single_size;
	if (_Max_degree < _Degree.calc_total_deg(_LMdeg))
	{
		_Max_degree = _Degree.calc_total_deg(_LMdeg);
	}
	//更新ここまで

	//スライド処理
	vector<unsigned char> coeff_temp(_Coeff_size);
	for (int i = 0; i < _Coeff_size; i++)
	{
		if (_Coeff[i] != 0)
		{
			coeff_temp[_Degree.degree_to_index(_Degree.vec_add(_Degree.index_to_degree(i), monomial_deg))] = _Coeff[i];
		}
	}
	_Coeff = coeff_temp;
}
*/